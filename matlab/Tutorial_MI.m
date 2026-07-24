%% Example script: Phase-amplitude coupling using Modulation Index
%
% This script demonstrates how to calculate phase-amplitude coupling
% between respiration and global field power in neural resting state
% data using the Modulation Index introduced by Tort (2010)
% as described in [ref].
%
% Input:
%   - z-scored respiration time series vector
%   - preprocessed continuous neural data (nsignals x timepoints)
%   -> usually both inputs are at same sampling frequency (concurrent
%      recordings), already downsampled brain data provided here for storage
%      space reasons
%
% Output:
%   mis.mat     = empirical MIs, normalized on surrogate MI distribution and
%                 expressed in units of SD (nsubj x nfreqs)
%   surrmis.mat = surrogate MIs (nsubjs x nfreqs x niter)
%
% Required Helper Functions:
%   - two_point_interp: extracts phase of a signal by two-point interpolation
%   - generate_surrogate_iaaft: creates a phase-scrambled surrogate time series
%                               of a signal
%   - PLV: computes phase-locking value of two phase time series
%
% Copyright (C) 2026, Teresa Berther University of Münster, Germany

clearvars
clc
close all

addpath(genpath('~/respmethods/'));
datapath = '~/respmethods/_exampledata/';                                   % contains raw respiration traces and neural data
outpath = '~/respmethods/_out/';

% subject ids
ids = {'004' '005' '006' '007' '009' '010' '011' '012' '013' '014' '015' '016' '019' '020' '022' '023'};

% config
orgfs = 300;                                                                % original sampling frequency of raw data
fs = 100;                                                                   % final sampling frequency after downsampling

niter = 5000;                                                               % number of surrogate iterations
plvcrit = 0.1;                                                              % threshold for IAAFT surrogate generation, max PLV of surrogate with data
nbin = 60;                                                                  % number of bins for MI calculation
pb = linspace(-pi,pi,nbin);                                                 % vector containing centre of each phase bin (= phase bin vector)
nfreqs = 64;                                                                % number of frequencies for MI calculation, derived from wavelet filter bank

% set up returned structures
mis = nan(length(ids), nfreqs);
surrmis = nan(length(ids), nfreqs, niter);

% iterate across subjects
for isub = 1:length(ids)

    %% Extract respiratory phase and generate surrogates
    disp(['Computing phase and surrogates for subject #' num2str(isub) '/' num2str(length(ids))]);

    pattern = sprintf('*resp_%s.mat', ids{isub});
    respfiles = dir(fullfile(datapath, pattern));
    load([datapath respfiles.name]);

    x = downsample(data, orgfs/fs);                                         % downsample respiration signal
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

    %% Prepare neural data
    pattern = sprintf('*meg_%s.mat', ids{isub});
    neurofiles = dir(fullfile(datapath, pattern));
    load([datapath neurofiles.name]);

    % resbrain = downsample(brain.', orgfs/fs).';                           % downsample neural data from original sampling freq, here directly loading resampled data  

    gfp = zeros(nfreqs, size(resbrain, 2));                                 % struct for global field power, size frequencies x nsignals
    fb = cwtfilterbank('SignalLength', size(resbrain, 2), ...               % create wavelet filter bank to extract amplitude envelopes
        'SamplingFrequency', fs, 'FrequencyLimits',[0.5 40], ...            % the parameters FrequencyLimits and WaveletParams used here
        'WaveletParameters',[3,20]);                                        % determine the frequency resolution (here: 64 freqs between 0.5 - 40Hz)

    % calculate global field power for each channel/sensor/parcel/voxel
    for ichan = 1:size(resbrain, 1)
        [Envc,f] = wt(fb, resbrain(ichan,:));                               % apply wavelet transform
        y = abs(Envc);
        gfp = gfp + abs(y);                                                 % calculate global field power as sum of amplitude envelopes
    end

    %% Calculate empirical and surrogate Modulation Index
    % remove NaNs first: all phase extraction methods introduce NaNs
    % in beginning & end of phase-vector, we find the common valid
    % timepoints across surrogate and empirical phase time series
    nancols = any(isnan(allsresp), 1);
    first_surr = find(~nancols, 1, 'first');                                % find last column with NaNs at beginning across surrogates
    last_surr  = find(~nancols, 1, 'last');                                 % find first column with NaNs at end across surrogates

    first = max(first, first_surr);                                         % match against empirical NaN indices
    last  = min(last, last_surr);

    pv = pv(first:last);                                                    % NAs need to be removed for binning in mi code
    allsresp = allsresp(:, first:last);                                     % remove the same indices from the surrogate data
    gfp = gfp(:, first:last)';                                              % match indices of neural data

    %%% calculate empirical MI
    lnbin = log(nbin);
    edg = eqpop(pv, nbin);                                                  % work on the whole array
    [~,bin] = histc(pv, edg);
    binamp = zeros(size(gfp, 2), nbin);                                     % binned amplitude

    nip = length(find(bin==1));                                             % samples per bin

    for k = 1:nbin
        binamp(:,k) = mean(gfp(bin==k,:),1);                                % calculate mean amplitude for each bin
    end
    binampn = binamp./repmat(sum(binamp,2),[1 nbin]);
    MI = (lnbin-(-sum(binampn.*log(binampn),2)))/lnbin;                     % this is our empirical MI value

    %%% calculate surrogate MIs
    MI2 = zeros(niter, nfreqs);
    for k = 1:niter
        tmp = allsresp(k,:)';
        edg = eqpop(tmp,nbin);
        [~,bin] = histc(tmp,edg);
        binamp = zeros(size(gfp,2),nbin);

        nip = length(find(bin==1));

        for k2 = 1:nbin
            binamp(:,k2) = mean(gfp(bin==k2,:),1);
        end
        binampn = binamp./repmat(sum(binamp,2),[1 nbin]);
        MI2(k,:)=(lnbin-(-sum(binampn.*log(binampn),2)))/lnbin;             % this is the MI for the surrogate signal
    end

    bm = mean(MI2,1);
    bs = std(MI2,0,1);
    mi = (MI'-bm)./bs;                                                      % normalise empirical MI values with surrogate distribution

    %%% save everything
    mis(isub, :) = mi;
    surrmis(isub, :, :) = MI2';

end

save([outpath 'misemp.mat'], "mis");                  % MIs for all subjects
save([outpath 'surrmis.mat'], "surrmis", "-v7.3");    % surrogate MIs for all subjects
