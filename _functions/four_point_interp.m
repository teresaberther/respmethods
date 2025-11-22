function [phasevec, first, last] = four_point_interp(zresp)

% Four_point_interp extracts the phase of a oscillatory time series 
% (e.g. respiration) using the four-point interpolation method.
% 
%
% Use as:
%   [phasevec, first, last] = four_point_interp(ts)
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

[~, peaks] = findpeaks(zresp, 'MinPeakProminence',1); % the last argument uses a criterion of z = 1 for peak prominence
troughs = [];
inflinsp = [];
inflexp = [];

for k = 2:numel(peaks)                           % start with peak #2
    tmp = zresp(peaks(k-1):peaks(k));            % get respiration course between peaks #1 and #2
    inflpt1 = find(diff(tmp) == max(diff(tmp))); % inflection point 1 (= strongest increase in this cycle)
    inflpt2 = find(diff(tmp) == min(diff(tmp))); % inflection point 2 (= strongest decrease in this cycle)
    inds = peaks(k-1):peaks(k);                  % get indices of the respiration course
    troughinds = inds(tmp == min(tmp));          % minimum of respiration between peaks is the trough
    troughs(k-1) = troughinds(1);                % for the rare but annoying case there is a peak longer than 1 sample
    inflinsp(k-1) = inds(inflpt1);
    inflexp(k-1) = inds(inflpt2);
end

phasevec = NaN(size(zresp));            % initiate phase vector with NaNs
phasevec(peaks) = 0;                    % resp maximum = phase zero
phasevec(troughs) = pi;                 % resp minimum = phase +/- pi
phasevec(inflinsp) = -pi/2;             % inspiratory infliction point = phase -pi/2
phasevec(inflexp) = pi/2;               % expiratory infliction point = phase pi/2

for k = 1:numel(peaks)-1
    t1 = phasevec(peaks(k)+1:inflexp(k));               % peak to expiratory infliction point
    s1 = linspace(0+pi/numel(t1),pi/2,numel(t1));       % linear interpolation
    phasevec(peaks(k)+1:inflexp(k)) = s1;               % substitute with phase angles (peak2inflexp)

    t2 = phasevec(inflexp(k):troughs(k));               % expiratory infliction point to trough
    s2 = linspace(pi/2+pi/numel(t2),pi,numel(t2));
    phasevec(inflexp(k):troughs(k)) = s2;

    t3 = phasevec(troughs(k):inflinsp(k));              % trough to inspiratory infliction point
    s3 = linspace(-pi+pi/numel(t3),-pi/2,numel(t3));
    phasevec(troughs(k):inflinsp(k)) = s3;

    t4 = phasevec(inflinsp(k):peaks(k+1));              % inspiratory infliction point to peak
    s4 = linspace(-pi/2+pi/numel(t4),0,numel(t4));
    phasevec(inflinsp(k):peaks(k+1)) = s4;
end

first = find(~isnan(phasevec), 1, 'first');   % get index of first phase value so we can exclude NAs later
last = find(~isnan(phasevec), 1, 'last');     % get index of last phase value so we can exclude NAs later

end