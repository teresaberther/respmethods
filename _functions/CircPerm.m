function [summary, bounds] = CircPerm(EMPdat, PERMdat, theta, alternative)

% CircPerm performs a cluster permutation test on circular data.
%
%
% Use as:
%   [summary, bounds] = CircPerm(EMPdat, PERMdat, theta, 'two_sided')
%
%
% Input arguments:
%   EMPdat      = vector of empirical data across phase bins to be tested 
%   PERMdat     = phase bins x permutation iterations matrix of surrogate data
%   theta       = a vector of angles (in radians) from -pi to pi/0 to 2*pi, representing the phase bins of the data (phase bins x 1)
%   alternative = string specifying alternative hypothesis, "two_sided", "greater" or "less" (default = "two_sided")
%
% Output:
%   summary            = structure containing the following fields:
%     .idxs            = indices of the phase bins for each cluster, each row represents a separate cluster 
%     .clustmass_stat  = clustering statistic (i.e. sum of z-scores) for each cluster
%     .p_ecdf          = empirical cumulative distribution value for each cluster statistic (proportion of surrogate statistics <= observed value)
%     .p               = permutation p-value (proportion of absolute surrogate stats values >= observed values)
%   
%   bounds             = empirical boundaries used for clustering in each bin, based on ecdf distribution
%
% Copyright (C) 2025, Elio Balestrieri & Teresa Berther, University of Münster, Germany 

% check input
if size(theta, 1) == 1
    theta = theta';
    warning('Given angle/phase bin vector did not match required dimensions, was transposed automatically - check input.')
end 
if any(isnan(EMPdat), 'all') || any(isnan(PERMdat), 'all')
    error('Data contains NaN values. Remove NaNs before performing the cluster permutation test.');
end
if any(theta < -pi) || any(theta > 2*pi)
    error('Phase values are outside the expected range (-π to π or π to 2π). Use phase values in radians and scale to expected range.');
end
if length(EMPdat) ~= length(theta)
    error('Dimension mismatch in empirical data: expected input is a vector of empirical data across the given phase bins.');
end
if size(PERMdat, 1) ~= length(theta)
    error('Dimension mismatch in surrogate data: expected input is a matrix of size phase bins x perm iterations.');
end

% set boundaries for cluster definition according to alternative hypothesis
if nargin < 4 || isempty(alternative) || strcmp(alternative, "two_sided")
    low_p = 0.025;
    up_p = 0.975;
elseif strcmp(alternative, "greater")
    low_p = 0;
    up_p = 0.95;
elseif strcmp(alternative, "lesser")
    low_p = 0.05;
    up_p = 1;
end 


% config
PERMEMPdat = [PERMdat, EMPdat]; 
nbins = length(theta); 
nperms = size(PERMdat, 2);

% find empirical boundaries
[up_bounds, low_bounds] = deal(nan(nbins, 1));

for ibin = 1:nbins

    this_bin = PERMEMPdat(ibin, :);
    [ECDF, val] = ecdf(this_bin);
    tmp_low = val(ECDF <= low_p);
    low_bounds(ibin) = tmp_low(end);
    tmp_up = val(ECDF >= up_p);
    up_bounds(ibin) = tmp_up(1);

end

%%% perform the permutations using surrogate data
surrogatedist = nan(nperms, 1); % store maximum clustering statistics from each permutation
for iperm = 1:nperms
    
    % standardize the permuted data by computing the z-score
    thisperm = PERMEMPdat(:, iperm);
    
    % cluster the permuted data
    tmp_out = CircClust(theta, thisperm, low_bounds, up_bounds, 'empirical');
    
    % combine the clustering statistics from both low and high clusters
    all_cluststats = [tmp_out.low_cluststats, tmp_out.up_cluststats];

    % select the maximum clustering statistic (by absolute value)
    if ~isempty(all_cluststats)
        [~, maxabs_idx] = max(abs(all_cluststats));
        surr_val = all_cluststats(maxabs_idx);
    else
        % in case there are no valid clusters, use the max absolute z-score
        [~, maxabs_idx] = max(abs(thisperm));
        surr_val = thisperm(maxabs_idx);
    end

    % store the surrogate statistic
    surrogatedist(iperm) = surr_val;
end


%%% clustering statistic for actual data
actual_avg = EMPdat;
out = CircClust(theta, actual_avg, low_bounds, up_bounds, 'empirical');

% combine the clustering statistics from both low and high clusters
all_cluststats = [out.low_cluststats, out.up_cluststats];
all_clustidxs = [out.low_clusts_idxs, out.up_clusts_idxs];


%%% compute p-values for each cluster by comparing with surrogate statistics
summary = struct('idxs', {}, 'clustmass_stat', [], 'p_ecdf', [], 'p', []);
acc = 1;
for imassstat = all_cluststats

    % combine surrogate & actual clustering statistic +inf to always have at least one value above the cluster stat  
    tmpstats = [surrogatedist; imassstat; inf]; 
    % compute the p-value as the proportion of surrogate statistics less than or equal to the actual statistic
    p_ecdf = mean(tmpstats <= imassstat);
    p = mean(abs(surrogatedist) >= abs(imassstat)); 

    summary(acc).idxs = all_clustidxs(acc);
    summary(acc).clustmass_stat = imassstat;
    summary(acc).p_ecdf = p_ecdf;
    summary(acc).p = p;
    
    acc = acc + 1;
end

bounds = [low_bounds, up_bounds];

end
