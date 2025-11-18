function [phasevec, first, last] = two_point_interp(ts)

% Two_point_interp extracts the phase of a oscillatory time series 
% (e.g. respiration) using the two-point interpolation method.
% 
%
% Use as:
%   [phasevec, first, last] = two_point_interp(ts)
%
% 
% Input arguments: 
%    ts = vector of zscored, possibly smoothed oscillatory time series (e.g. respiration trace) 
% 
% Output:  
%   phasevector = vector of phase value time series corresponding to the input timeseries
%   first       = first real phase value (method introduces NaNs at start of phase vector)
%   last        = last real phase value (method introduces NaNs at end of phase vector)
%
% Copyright (C) 2025, Daniel Kluger & Teresa Berther, University of Münster, Germany 

[~, peaks] = findpeaks(ts, 'MinPeakProminence',1); % use a criterion of z = 1 for peak prominence
troughs = [];
for k = 2:numel(peaks)                  % start with peak #2
    tmp = ts(peaks(k-1):peaks(k));      % get respiration course between peaks #1 and #2
    inds = peaks(k-1):peaks(k);         % get indices of the respiration course
    troughinds = inds(tmp == min(tmp)); % minimum of respiration between peaks is the trough
    troughs(k-1) = troughinds(1);       % for the rare but annoying case there is a peak longer than 1 sample
end

phasevec = NaN(size(ts));               % initiate phase vector with NaNs
phasevec(peaks) = 0;                    % inspiration maximum = phase zero
phasevec(troughs) = pi;                 % inspiration minimum = phase +/- pi

for k = 1:numel(peaks)-1
    tmp = phasevec(peaks(k)+1:troughs(k));              % phase angles between peak #1  and trough #1
    tmp2 = phasevec(troughs(k)+1:peaks(k+1));           % phase angles between trough #1 and peak #2
    sub = linspace(0+pi/numel(tmp),pi,numel(tmp));      % linear interpolation peak2trough
    sub2 = linspace(-pi+pi/numel(tmp2),0,numel(tmp2));  % same for trough2peak
    phasevec(peaks(k)+1:troughs(k)) = sub;              % substitute with phase angles (peak2trough)
    phasevec(troughs(k)+1:peaks(k+1)) = sub2;           % same for trough2peak
end

first = find(~isnan(phasevec), 1, 'first');   % get index of first phase value so we can exclude NaNs later
last = find(~isnan(phasevec), 1, 'last');     % get index of last phase value so we can exclude NaNs later

end
