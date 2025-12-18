%% Example script: Permutation test of circular data
% 
% This script demonstrates how to run clustering and permutation
% testing of real and simulated circular data as described in [ref].
% 
% An example dataset of behavioural hit rates in a spatial attention 
% paradigm across respiration phase bins is available at [github_link].
%
% Copyright (C) 2025, Teresa Berther & Elio Balestrieri, University of Münster, Germany 

%% Cluster permutation test on example dataset
clearvars
clc
close all

addpath(genpath('~/respmethods/'));
outpath = '~/respmethods/_out/';

load(fullfile(outpath, 'binhr.mat'));         % subjects x phase bins matrix of mean empirical hit rates across phase bins
load(fullfile(outpath, 'surrbinhr.mat'));     % subjects x phase bins x iterations matrix of surrogate hit rates across phase bins
load(fullfile(outpath, 'phasebinvect.mat'));  % vector of phase bins corresponding to data phase bins (phase bins x 1), range -pi to pi

% config
nbins = size(binhr, 2);
nsubj = size(binhr, 1); 
nsurr = size(sbinhr, 3);

% some pretty colors for plots (adapted from Wes Anderson Zissou1 palette, (C) Karthik Ram)
negcolors = [85, 152, 175;
             133, 182, 195;
             65, 130, 160] / 255;

poscolors = [230, 204, 79;
             218, 176, 59;
             249, 233, 148] / 255;

% z-score values
permemp_mat = cat(3, sbinhr, binhr);                                        
bsl_ = repmat(mean(permemp_mat, 3), 1, 1, nsurr+1);                        
scale_ = repmat(std(permemp_mat, [], 3), 1, 1, nsurr+1);                   
tmp_z = (permemp_mat-bsl_)./scale_;                                         

% calculate t-values 
permemp_t = squeeze(sqrt(nsubj)*mean(tmp_z, 1)./std(tmp_z, [], 1));        
emp_t = permemp_t(:, end);                                                        % empirical t-values
perm_t = permemp_t(:, 1:end-1);                                                   % permuted t-values


%%% run the cluster permutation test 
[summary, bounds] = CircPerm(emp_t, perm_t, pb, 'two_sided');       % two sided testing used here, use 'lesser'/'greater' for one sided tests 

%%%% plot significant clusters (linear plot)                                           
figure; 
hold on;
plot(pb, perm_t, 'Color', [0.85 0.85 0.85]);                        % permutation values
plot(pb, bounds(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);     % lower empirical cluster threshold
plot(pb, bounds(:, 2), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);     % upper empirical cluster threshold
plot(pb, emp_t, 'Color', 'k', 'LineWidth', 2);                      % empirical values

poscolor_idx = 1;
negcolor_idx = 1;

for iclust = 1:length(summary)
    if summary(iclust).p <= 0.05
        this_idxs = summary(iclust).idxs{:};
        if summary(iclust).p_ecdf <= 0.025      % negative clusters
            color = negcolors(min(negcolor_idx, size(negcolors, 1)), :);
            negcolor_idx = negcolor_idx + 1;
        elseif summary(iclust).p_ecdf >= 0.975  % positive clusters
            color = poscolors(min(poscolor_idx, size(poscolors, 1)), :);
            poscolor_idx = poscolor_idx + 1;
        else
            continue
        end

        scatter(pb(this_idxs), emp_t(this_idxs), 90, 'filled', ...
            'MarkerFaceColor', color);                                            % color points belonging to cluster

        if ismember(1, this_idxs) && ismember(nsteps, this_idxs)                  % handle wraparound clusters
            start = this_idxs(this_idxs >= nbins/2);
            last  = this_idxs(this_idxs <  nbins/2);
            plot(pb(start), repmat(6, size(start)), 'LineWidth', 8, 'Color', color);
            plot(pb(last),  repmat(6, size(last)),  'LineWidth', 8, 'Color', color);
        else
            plot(pb(this_idxs), repmat(6, size(this_idxs)), 'LineWidth', 8, 'Color', color);
        end
    end
end
xlim([-pi pi]);
xticks([-pi 0 pi]);
xticklabels({'-\pi', '0', '\pi'});
xlabel(['phase (' num2str(length(pb)) ' bins)'], 'FontSize', 15);
ylabel('tvals', 'FontSize', 15);
ylim([-6.5 6.5]);
yticks([-5 -2.5 0 2.5 5]);
set(gca, 'FontSize', 15); 
title('Significant clusters in hit rate ~ respiration phase', 'FontSize', 20);

%%%% plot significant clusters as a polar plot 
% wrap around data for polar plot
pb2 = [pb; pb(1)];
emp_t2 = [emp_t; emp_t(1)];
perm_t2 = [perm_t; perm_t(1, :)];
bounds2 = [bounds; bounds(1, :)];

% plot
figure;
ax = polaraxes; 
hold on;
polarscatter(pb2', perm_t2, 8, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9]);
polarplot(pb2, bounds2(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);    % lower empirical cluster threshold
polarplot(pb2, bounds2(:, 2), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);    % upper empirical cluster threshold
polarscatter(pb2', emp_t2, 100, 'filled', 'MarkerFaceColor', 'k');

poscolor_idx = 1;
negcolor_idx = 1;

for iclust = 1:length(summary)

    if summary(iclust).p <= 0.05
        this_idxs = summary(iclust).idxs{:};
        if summary(iclust).p_ecdf <= 0.025      % negative clusters
            color = negcolors(min(negcolor_idx, size(negcolors, 1)), :);
            negcolor_idx = negcolor_idx + 1;
        elseif summary(iclust).p_ecdf >= 0.975  % positive clusters
            color = poscolors(min(poscolor_idx, size(poscolors, 1)), :);
            poscolor_idx = poscolor_idx + 1;
        else
            continue
        end

        polarscatter(pb2(this_idxs), emp_t2(this_idxs), 100, 'filled', ...
            'MarkerFaceColor', color);                                                      % color points belonging to cluster

        if ismember(1, this_idxs) && ismember(nbins, this_idxs)                             % handle wraparound clusters
            start = this_idxs(this_idxs >= nbins/2);
            last  = this_idxs(this_idxs <  nbins/2);
            polarplot(pb2(start), repmat(4.8, size(start)), 'LineWidth', 8, 'Color', color);
            polarplot(pb2(last),  repmat(4.8, size(last)),  'LineWidth', 8, 'Color', color);
        else
            polarplot(pb2(this_idxs), repmat(4.8, size(this_idxs)), 'LineWidth', 8, 'Color', color);
        end
    end
end

ax.ThetaZeroLocation = 'left'; 
ax.ThetaTick = [0 180];
ax.ThetaTickLabels = {'0', '-\pi/\pi'};
ax.RTick = [];
rlim([-6.5 5]); 
set(ax, 'FontSize', 15); 
title('Significant clusters in hit rate ~ respiration phase', 'FontSize', 20);



%% Cluster permutation test on some simulated data 
% config
nsteps = 60;        % continuous phase segmented into 60 phase bins 
nsubjs = 20;        % 20 subjects
nperms = 10000;     % 10000 permutations

% some pretty colors for plots (adapted from Wes Anderson Zissou1 palette, (C) Karthik Ram)
negcolors = [85, 152, 175;
             133, 182, 195;
             65, 130, 160] / 255;

poscolors = [230, 204, 79;
             218, 176, 59;
             249, 233, 148] / 255;

% simulate some random empirical data
anglevect = linspace(-pi, pi, nsteps)';
anglemat = repmat(anglevect, 1, nsubjs);
phasediff_mat = repmat(randn(nsteps, 1)*pi/6, 1, nsubjs);
freqdiff_mat = repmat(randn(nsteps, 1)*.1+1, 1, nsubjs);
simemp = (cos(anglemat.*freqdiff_mat + phasediff_mat) + randn(nsteps, nsubjs)*.3)';

% circularly shift empirical data for each subject to simulate surrogate data
simperm = nan(nsubjs, nsteps, nperms);
for iperm = 1:nperms
    surr = nan(nsubjs, nsteps);
    for isub = 1:nsubjs
        shift = randi(nsteps);  % random circular shift
        surr(isub, :) = circshift(simemp(isub, :), shift);
    end
    simperm(:, :, iperm) = surr;
end

% calculate t-vals
permemp_mat = cat(3, simperm, simemp);                                         
permemp_t = squeeze(sqrt(nsubjs)*mean(permemp_mat, 1)./std(permemp_mat, [], 1));        
emp_t = permemp_t(:, end);                           
perm_t = permemp_t(:, 1:end-1);  

%%% run the cluster permutation test
[summary, bounds] = CircPerm(emp_t, perm_t, anglevect, 'two_sided');       % two sided testing used here, use 'lesser'/'greater' for one sided tests 

%%%% plot significant clusters as linear plot                                     
figure; 
hold on;
plot(anglevect, perm_t, 'Color', [0.85 0.85 0.85]);                        % permutation values
plot(anglevect, bounds(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);     % lower empirical cluster threshold
plot(anglevect, bounds(:, 2), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);     % upper empirical cluster threshold
plot(anglevect, emp_t, 'Color', 'k', 'LineWidth', 2);                      % empirical values

poscolor_idx = 1;
negcolor_idx = 1;

for iclust = 1:length(summary)
    if summary(iclust).p <= 0.05
        this_idxs = summary(iclust).idxs{:};
        if summary(iclust).p_ecdf <= 0.025      % negative clusters
            color = negcolors(min(negcolor_idx, size(negcolors, 1)), :);
            negcolor_idx = negcolor_idx + 1;
        elseif summary(iclust).p_ecdf >= 0.975  % positive clusters
            color = poscolors(min(poscolor_idx, size(poscolors, 1)), :);
            poscolor_idx = poscolor_idx + 1;
        else
            continue
        end

        scatter(anglevect(this_idxs), emp_t(this_idxs), 90, 'filled', ...
            'MarkerFaceColor', color);                                            % color points belonging to cluster

        if ismember(1, this_idxs) && ismember(nsteps, this_idxs)                  % handle wraparound clusters
            start = this_idxs(this_idxs >= nsteps/2);
            last  = this_idxs(this_idxs <  nsteps/2);
            plot(anglevect(start), repmat(29, size(start)), 'LineWidth', 8, 'Color', color);
            plot(anglevect(last),  repmat(29, size(last)),  'LineWidth', 8, 'Color', color);
        else
            plot(anglevect(this_idxs), repmat(29, size(this_idxs)), 'LineWidth', 8, 'Color', color);
        end
    end
end
xlim([-pi pi]);
xticks([-pi 0 pi]);
xticklabels({'-\pi', '0', '\pi'});
xlabel(['phase (' num2str(length(anglevect)) ' bins)'], 'FontSize', 15);
ylabel('tvals', 'FontSize', 15);
ylim([-30 30]);
set(gca, 'FontSize', 15); 
title('Significant clusters in simulated data', 'FontSize', 20);

%%%% plot significant clusters as a polar plot 
% wrap around data for polar plot
anglevect2 = [anglevect; anglevect(1)];
emp_t2 = [emp_t; emp_t(1)];
perm_t2 = [perm_t; perm_t(1, :)];
bounds2 = [bounds; bounds(1, :)];

% plot
figure;
ax = polaraxes; 
hold on;
polarscatter(anglevect2', perm_t2, 8, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9]);
polarplot(anglevect2, bounds2(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);      % lower empirical cluster threshold
polarplot(anglevect2, bounds2(:, 2), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);      % upper empirical cluster threshold
polarscatter(anglevect2', emp_t2, 100, 'filled', 'MarkerFaceColor', 'k');
polarplot(anglevect2', emp_t2, 'Color', 'k', 'LineWidth', 2);

poscolor_idx = 1;
negcolor_idx = 1;

for iclust = 1:length(summary)
    if summary(iclust).p <= 0.05
        this_idxs = summary(iclust).idxs{:};
        if summary(iclust).p_ecdf <= 0.025      % negative clusters
            color = negcolors(min(negcolor_idx, size(negcolors, 1)), :);
            negcolor_idx = negcolor_idx + 1;
        elseif summary(iclust).p_ecdf >= 0.975  % positive clusters
            color = poscolors(min(poscolor_idx, size(poscolors, 1)), :);
            poscolor_idx = poscolor_idx + 1;
        else
            continue
        end

        polarscatter(anglevect2(this_idxs), emp_t2(this_idxs), 100, 'filled', ...
            'MarkerFaceColor', color);                                                      % color points belonging to cluster

        if ismember(1, this_idxs) && ismember(nsteps, this_idxs)                            % handle wraparound clusters
            start = this_idxs(this_idxs >= nsteps/2);
            last  = this_idxs(this_idxs <  nsteps/2);
            polarplot(anglevect2(start), repmat(29, size(start)), 'LineWidth', 8, 'Color', color);
            polarplot(anglevect2(last),  repmat(29, size(last)),  'LineWidth', 8, 'Color', color);
        else
            polarplot(anglevect2(this_idxs), repmat(29, size(this_idxs)), 'LineWidth', 8, 'Color', color);
        end
    end
end

ax.ThetaZeroLocation = 'left'; 
ax.ThetaTick = [0 180];
ax.ThetaTickLabels = {'0', '-\pi/\pi'};
ax.RTick = [];
rlim([-30 30]);
set(ax, 'FontSize', 15); 
title('Significant clusters in simulated data', 'FontSize', 20);

