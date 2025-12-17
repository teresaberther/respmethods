%% Data preparation script: respiratory phase extraction, surrogate generation & hitrate binning
%
% This script extracts respiratory phase from raw respiration data using
% two-point interpolation, generates surrogate phases via the Iterative
% Amplitude-Adjusted Fourier Transform (IAAFT), and bins behavioral data
% (e.g., hit rates) into empirical and surrogate respiratory phase bins.
% A detailed description of the preparation pipeline can be found in
% [ref].
%
% Input:
%   - raw respiration time series vector
%   - events structure with stimulus onset and outcome field for each subject
%     .onset: contains timestamp of stimulus presentation for each trial (in samples from start of recording)
%     .outcome: contains outcome value for each trial (here: hit vs miss)
%
% Output (here: hit rate):
%   binhr.mat   = mean empirical outcome for each respiratory phase bin, for each subject
%                 (subjects x phase bins matrix)
%   sbinhr.mat  = mean surrogate outcome for each respiratory phase bin and surrogate, for each subject
%                 (subjects x phase bins x surrogate iterations matrix) 
%
% Required Helper Functions:
%   - two_point_interp: extracts phase of a signal by two-point interpolation
%   - generate_surrogate_iaaft: creates a phase-scrambled surrogate time series
%                               of a signal
%   - PLV: computes phase-locking value of two phase time series
%
%
% Copyright (C) 2025, Daniel Kluger & Teresa Berther, University of Münster, Germany

clearvars
clc
close all

addpath(genpath('~/respmethods/'));
datapath = '~/respmethods/_data/';                                          % contains raw respiration traces and events table
outpath = '~/respmethods/_out/';


% subject ids
ids = {'002' '003' '004' '005' '006' '007' '009' '010' '011' '012' '013' '014' '015' '016' '019' '020' '022' '023'};

% config
orgfs = 1000;                                                               % original sampling frequency of raw data
fs = 100;                                                                   % final sampling frequency after downsampling

niter = 5000;                                                               % number of surrogate iterations
plvcrit = 0.1;                                                              % threshold for IAAFT surrogate generation, max PLV of surrogate with data
nbin = 60;                                                                  % number of phase bins
phw = 2*pi/10;                                                              % width of each phase bin
pb = linspace(-pi,pi,nbin);                                                 % vector containing centre of each phase bin (= phase bin vector)

load(fullfile(datapath, 'events.mat'));                                     % load stimulus onsets ("onsets") and outcome (here: "corr" for hit vs miss trials)

%% Compute true and surrogate respiration phase values at stimulus onset 
for isub = 1:numel(ids)

    disp(['Computing phase and surrogates for subject #' num2str(isub) '/' num2str(length(ids))]);

    pattern = sprintf('*_%s.mat', ids{isub});
    respfiles = dir(fullfile(datapath, pattern));
    load([datapath respfiles.name]);

    x = downsample(zscore(data(1,:)), orgfs/fs);                            % downsample & zscore raw respiration signal
    x = movmean(x, 0.4*fs);                                                 % some smoothing

    % extract phase vector of empirical data using two-point interpolation
    [pv, first, last] = two_point_interp(x);                                % first/last indices saved for NaN exclusion later

    %%% generate surrogate phase vectors using IAAFT
    allsresp = nan(niter,length(x));
    k = 1;
    while k <= niter

        disp(['Computing surrogate resp trace #' num2str(k) '/' num2str(niter)]);

        sresp = generate_surrogate_iaaft(x,'verbose',false);                % single mock solution for respiration time series
        [spv, ftmp, ltmp] = two_point_interp(sresp);                        % extract surrogate phase via two-point interpolation
        spvtmp = spv(max(first, ftmp):min(last, ltmp));                     % temporarily get rid of NaNs in phase vectors as PLV function cannot handle NaNs
        pvtmp = pv(max(first, ftmp):min(last, ltmp));
        plv = PLV(pvtmp', spvtmp);                                          % compute phase-locking value
        if plv < plvcrit                                                    % only keep solution if plv is smaller than current constraint
            allsresp(k,:) = spv;                                            % save full surrogate phase vector
            k = k+1;
        else
            continue                                                        % otherwise start over
        end
    end
    save(fullfile(outpath, ['surrogates_iaaft_' num2str(isub) '.mat']),'allsresp', '-v7.3');

    % get stimulus onsets in the respiration time frame
    % (= time from recording start in ms)
    disp(['Computing empirical and surrogate onsets for subject #' num2str(isub) '/' num2str(length(ids))]);

    onsets = events{isub}.onsets;                                           % the onsets field contains the sample time stamps for the stimuli presentation
    onsets = round(onsets/(orgfs/fs));                                      % downsample the sample time stamps to match new fs
    empphase = pv(onsets);                                                  % empirical phase at onset of stimulus for each trial
    surrphases = [];
    for k = 1:niter
        surrphases(:,k) = allsresp(k, onsets);                              % surrogate phase at onset of stimulus for each trial
    end

    % save
    save(fullfile(outpath, ['empiricalphase_' num2str(isub) '.mat']),'empphase');
    save(fullfile(outpath, ['surrogatephases_' num2str(isub) '.mat']),'surrphases');

end


%% Compute outcome (hit rates) ~ respiration phase
% example outcome here is hit rate, but can be any outcome variable
binhr = nan(numel(ids), nbin);                                              
sbinhr = nan(numel(ids), nbin, niter);

for isub = 1:numel(ids)

    disp(['Computing empirical and surrogate binned hit rates for subject #' num2str(isub) '/' num2str(length(ids))]);

    % empirical binHR
    load(fullfile(outpath, ['empiricalphase_' num2str(isub) '.mat']),'empphase');
    for ibin = 1:nbin
        phsel = find((empphase>pb(ibin)-phw) & (empphase<pb(ibin)+phw) | ...     % find all trials whose phase falls within the current phase bin
            (empphase-2*pi>pb(ibin)-phw) & (empphase-2*pi<pb(ibin)+phw)| ...
            (empphase+2*pi>pb(ibin)-phw) & (empphase+2*pi<pb(ibin)+phw));
        binhr(isub,ibin) = sum(events{isub}.outcome(phsel))/numel(phsel);        % based on the outcome for these trials, compute hit rate for this phase bin
    end

    % surrogate binHR
    load(fullfile(outpath, ['surrogatephases_' num2str(isub) '.mat']),'surrphases');
    for k = 1:niter
        sphase = surrphases(:, k);                                              
        for ibin = 1:nbin
            phsel = find((sphase>pb(ibin)-phw) & (sphase<pb(ibin)+phw) | ...     % again find trials for each phase bin...
                (sphase-2*pi>pb(ibin)-phw) & (sphase-2*pi<pb(ibin)+phw)| ...
                (sphase+2*pi>pb(ibin)-phw) & (sphase+2*pi<pb(ibin)+phw));
            sbinhr(isub,ibin,k) = sum(events{isub}.outcome(phsel))/numel(phsel); % ...and compute bin HR based on the outcome
        end
    end
end

% save mean hit rates ~ respiration phase bins for empirical & surrogate data
% this is the data we run the circular clustering on
save(fullfile(outpath, 'binhr.mat'), 'binhr');                              % subjects x nbins array of empirical phase-binned outcome values
save(fullfile(outpath, 'surrbinhr.mat'), 'sbinhr');                         % subjects x nbins x niter array of surrogate phase-binned values
save(fullfile(outpath, 'phasebinvect.mat'), 'pb');                          % also save the vector with the phase bin centres
