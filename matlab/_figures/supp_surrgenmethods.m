%% Example script: Surrogate generation methods 
% 
% This script demonstrates the generation of surrogate phase time series 
% on one example time seires using the alternate procedures described in 
% Berther, T., Balestrieri, E., Saltafossi, M., Paulsen, L. B., Andersen, L. M., Kluger, D. S. (2025). 
% Robust circular cluster-based statistics for respiration-brain coupling.
% PsyArXiv, 2025-11.
%
% Input:
%   - z-scored, smoothed respiration time series vector
%
% Output:
%   - surrogate phase time series, k = number of iterations chosen 
%
% Required Helper Functions:
%   - two_point_interp: extracts phase of a signal by two-point interpolation
%   - generate_surrogate_iaaft: creates a phase-scrambled surrogate time series
%                               of a signal
%   - PLV: computes phase-locking value of two phase time series
%
% Copyright (C) 2026, Teresa Berther University of Münster, Germany 

% setup
addpath(genpath('~respmethods/matlab/'));
datapath = '~respmethods/matlab/_figures/_data4plotting/';

niter = 1000;                                                               % number of surrogate sets to create
plvcrit = 0.1;

% prepare the original respiration time series
load([datapath 'zresp_001.mat']);
[tmp, first, last] = two_point_interp(x); 
pv = tmp(first:last);                                                       % remove NaNs introduced by phase extraction so they are not shuffled into the surrogates

%% 1) Circular shifting
% circularly shift the phase time series in time to generate a surrogate
% time series, by a random shift with a set minimum 
allsresp = nan(niter,length(pv));
shifts = randi([30*fs, length(pv)-30*fs], niter, 1);                        % create niter random shifts, min. +/- 30s
for iter = 1:niter
    sresp = circshift(pv, shifts(iter));
    allsresp(iter, :) = sresp;
end

%% 2) Segment shuffling
% randomly shuffle segments of the phase vector to generate a surrogate
% phase time series 
allsresp = nan(niter,length(tmppv));
binl = 50;                                                                  % segment length in samples
nl = length(tmppv);
nseg = floor(nl/binl);
segments = reshape(tmppv(1:nseg*binl), binl, nseg);

for iiter = 1:niter
    order = randperm(nseg);
    shuffled = segments(:, order);
    shuffled = shuffled(:)';
    allsresp(iiter, 1:nseg*binl) = shuffled;

    if nseg * binl < nl
        leftover = tmppv(nseg*binl + 1:end);
        allsresp(iiter, nseg*binl + 1:end) = leftover;
    end
end

%% 3) Random shuffling 
% randomly shuffle all datapoints of the phase time series to generate a
% surrogate phase time series 
allsresp = nan(niter,length(pv));
for iiter = 1:niter
    allsresp(iiter, :) = pv(randperm(length(pv)));
end

%%
% Replace the lines for surrogate generation in Tutorial_MI
% or Tutorial_DataPrep to switch to the alternate surrogate generation 
% methods demonstrated here. 