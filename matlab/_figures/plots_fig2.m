%% Plotting script: Figure 2 
% From raw respiratory recordings to quasi-continuous binned phase
% 
% This script creates the subplots shown in Figure 2 of
% Berther, T., Balestrieri, E., Saltafossi, M., Paulsen, L. B., Andersen, L. M., Kluger, D. S. (2025). 
% Robust circular cluster-based statistics for respiration-brain coupling. PsyArXiv, 2025-11.
% 
% Input:
%  - raw respiration time series vector
%  - multiple versions of binned outcome variable (here: hit rate) based on
%    different binning parameters (cf. TutorialDataPrep)
%  
% Required Helper Functions/Scripts:
%  - TutorialDataPrep
%  - all interpolation functions from this toolbox 
%
% Copyright (C) 2026, Teresa Berther, University of Münster, Germany 

% setup
addpath(genpath('~respmethods/matlab/'));
datapath = '~respmethods/matlab/_figures/_data4plotting/';

colors = [217/255 216/255 134/255;
          159/255 193/255 120/255;
          135/255 198/255 194/255;
          103/255 150/255 189/255;
          30/255 83/255 128/255;
          20/255 46/255 84/255];

%% a) Phase extraction methods
% comparison of phase vectors extracted by different methods 
load([datapath 'zresp_001.mat']);                                          % z-scored, smoothed & downsampled (to 100Hz) respiration time series z_1331
fs = 100;                                                                   % sampling frequency is 100 Hz

% extract phase with different methods
hilb = angle(hilbert(x));
[tpt, f_tpt, l_tpt] = two_point_interp(x);
[fpt, f_fpt, l_fpt] = four_point_interp(x);
[pp, f_pp, l_pp] = protophase_interp(x);
[tpz, f_tpz, l_tpz] = trapez_interp(x, fs);

% deal with NaNs
firsts = [f_tpt, f_fpt, f_pp, f_tpz];
lasts = [l_tpt, l_fpt, l_pp, l_tpz];

f = max(firsts);
l = min(lasts);

% a) different interpolation methods plot
figure;
hold on;
plot(x(f:l), 'Color', 'k', 'LineWidth', 2); 
scatter(1:numel(x(f:l)), hilb(f:l), 8,'filled','MarkerFaceColor', colors(1,:, :),'MarkerFaceAlpha',0.75);
scatter(1:numel(x(f:l)), tpt(f:l), 8,'filled','MarkerFaceColor', colors(2,:, :),'MarkerFaceAlpha',0.75);
scatter(1:numel(x(f:l)), fpt(f:l), 8,'filled','MarkerFaceColor', colors(3,:, :),'MarkerFaceAlpha',0.75);
scatter(1:numel(x(f:l)), pp(f:l), 8,'filled','MarkerFaceColor', colors(4,:, :),'MarkerFaceAlpha',0.75);
scatter(1:numel(x(f:l)), tpz(f:l), 8,'filled','MarkerFaceColor', colors(5,:, :),'MarkerFaceAlpha',0.75);
xlim([3924 5069]);
axis off;

%% b) Potential artefacts
% show potential artefacts caused by hilbert phase extraction
load([datapath('zresp_002.mat')]);                                          % z-scored, smoothed & downsampled (to 100Hz) respiration time series z_1960

% extract phase
hilb = angle(hilbert(x));
[tpt, ~] = two_point_interp(x);

% plot
figure;
hold on;
plot(movmean(x, 0.4*fs), 'Color', 'k', 'LineWidth', 2);                % bit more smoothing for aesthetics
scatter(1:numel(x), hilb, 8,'filled','MarkerFaceColor', colors(1,:, :),'MarkerFaceAlpha',0.75);
scatter(1:numel(x), tpt, 8,'filled','MarkerFaceColor', colors(2,:, :),'MarkerFaceAlpha',0.75);
fill([88197 88708 88708  88197], [-pi-0.1 -pi-0.1 pi+0.1 pi+0.1], [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlim([87290 89593]);
axis off;

%% c) Phase binning 
% created in Power Point

%% d) Binning parameters 
% show effect of number of bins for two different phase bin widths 
nbins = [2, 4, 8, 15, 30, 60];

% phw = 2*pi/10
load([datapath 'binhrs_phw10.mat']);
load([datapath 'pbs_phw10.mat']);

figure;
hold on;
for ibin = 1:length(nbins)
    binhr = binhr_10{ibin};
    pb = pb_10{ibin};
    plot(pb, mean(binhr, 1), 'Color', colors(ibin, :, :), 'LineWidth', 4);
end
xlim([-pi pi]);
xticks([]);
ylim([0.52 0.64]);
yticks([0.52 0.64]);
% legend created manually in power point

% phw = 2*pi/5
load([datapath 'binhrs_phw5.mat']);
load([datapath 'pbs_phw5.mat']);

figure;
hold on;
for ibin = 1:length(nbins)
    binhr = binhr_5{ibin};
    pb = pb_5{ibin};
    plot(pb, mean(binhr, 1), 'Color', colors(ibin, :, :), 'LineWidth', 4);
end
xlim([-pi pi]);
xticks([]);
ylim([0.52 0.64]);
yticks([0.52 0.64]);
% legend created manually in power point
