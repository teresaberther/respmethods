function out = CircClust(theta, target, low_p, up_p, threshtype)

% CircClust performs clustering on circular data using the specified 
% thresholds.
%
% Use as:
%   [out] = CircClust(theta, target, low_p, up_p, threshtype)
%
%
% Input arguments:
%   theta      = a vector of angles (in radians), representing the phase of the target values
%   target     = a vector of target values corresponding to the theta values (data points to be clustered)
%   low_p      = lower p-value threshold (default =.025, two-sided testing)
%   up_p       = upper p-value threshold (default =.975, two-sided testing)
%   threshtype = thresholding type: 'percentile', 'empirical' or 'z' (default = 'empirical')
%
% Output:
%   out                 = structure containing the following fields:
%    .up_clusts_idxs    = indices of clusters above the upper threshold
%    .low_clusts_idxs   = indices of clusters below the lower threshold
%    .up_cluststats     = sum of target values for each cluster above the upper threshold
%    .low_cluststats    = sum of target values for each cluster below the lower threshold
%
% Copyright (C) 2025, Elio Balestrieri & Teresa Berther, University of Münster, Germany 

% check input
if any(isnan(theta)) || any(isnan(target))
    error('Data contains NaN values. Remove NaNs before clustering.');
end
if any(theta < -pi) || any(theta > 2*pi)
    error('Phase values for clustering are outside the expected range (-π to π or 0 to 2π). Use phase values in radians.');
end


if nargin < 5 || isempty(threshtype)
    threshtype = 'empirical'; 
end
if nargin < 4 || isempty(up_p)
    up_p = .975;
end
if nargin < 3 || isempty(low_p)
    low_p = .025;
end


% threshold based on the specified method
switch lower(threshtype)
    case 'percentile' 
        bounds = prctile(target, [low_p, up_p] * 100);
    case 'z'
        bounds = [norminv(min(low_p)), norminv(max(up_p))];
    case 'empirical'
        bounds = [low_p; up_p];    
    otherwise
        error('Unsupported threshold type. Use "empirical", "percentile" or "z".');
end

% classify points below and above the thresholds
below_ptile = target < min(bounds);
low_clusts = localClust(below_ptile, theta);

above_ptile = target > max(bounds);
up_clusts = localClust(above_ptile, theta);

% compute statistics for clusters above the upper threshold
up_cluststats = nan(size(up_clusts)); 
acc_clust = 1;
for iclust = up_clusts
    up_cluststats(acc_clust) = sum(target(iclust{:}));
    acc_clust = acc_clust + 1;
end

% compute statistics for clusters below the lower threshold
low_cluststats = nan(size(low_clusts)); 
acc_clust = 1;
for iclust = low_clusts
    low_cluststats(acc_clust) = sum(target(iclust{:}));
    acc_clust = acc_clust + 1;
end

% output structure with clustering results
out.up_clusts_idxs = up_clusts;
out.low_clusts_idxs = low_clusts;
out.up_cluststats = up_cluststats;
out.low_cluststats = low_cluststats;

end




function clusts = localClust(masklength, anglevect)

% localClust performs local clustering based on angular distances 
% between points, given a mask and an angle vector. It finds clusters 
% of contiguous points with similar angular distances.
%
% Use as:
%   [clusts] = localClust(masklength, anglevect)
%
%
% Input arguments:
%   masklength = binary vector or logical mask where 'true' indicates points to be considered for clustering
%   anglevect  = vector of angular data corresponding to the mask
%
% Output:
%   clusts     = cell array of clusters, each cell contains the indices of the points belonging to this cluster
%
% Copyright (C) 2025, Elio Balestrieri & Teresa Berther, University of Münster, Germany 


if ~any(masklength)
    clusts = {}; % return empty if no points are selected
else
    % convert the mask into a complex representation using angles
    cmplx_repr = masklength .* exp(anglevect * 1i);
    
    % compute pairwise distances between all points
    [X, Y] = meshgrid(cmplx_repr);
    dists = round(abs(X - Y), 5);  % calculate distance and round to 5 decimals
    minsigdist = unique(dists);    % find unique distances
    
    % identify the smallest non-zero distance to form clusters
    shortlenarc = dists == minsigdist(2);
    [row, col] = find(shortlenarc);  % find rows and columns of minimum distance
    
    % initialize cluster accumulation
    acc = 1; 
    clusts = {};
    
    while ~isempty(col)
        % accumulate points that belong to the same cluster
        idxs_acc = [];
        seed_entry = col(1); 
        while true
            idxs_acc = [idxs_acc; seed_entry]; 
            tmp_mask = ismember(col, seed_entry);
            
            if ~any(tmp_mask)
                break
            else
                seed_entry = unique(row(tmp_mask));
                row(tmp_mask) = [];
                col(tmp_mask) = [];
            end
        end
    
        clusts{acc} = unique(idxs_acc);
        acc = acc + 1;
    end
end

end
