%% Plotting script: Figure 3
% Implications of different surrogate procedures 
% 
% This script creates the subplots shown in Figure 3 of
% Berther, T., Balestrieri, E., Saltafossi, M., Paulsen, L. B., Andersen, L. M., Kluger, D. S. (2025). 
% Robust circular cluster-based statistics for respiration-brain coupling. PsyArXiv, 2025-11.
% 
% Input:
%  - z-scored, smoothed & downsampled respiration timeseries
%  - k = 5000 surrogate phase time series for each surrogate generation
%    method (see supplemental_surrgen_methods) based on the resp time series 
%  - empirical & surrogate Modulation Index spectra (see Tutorial_MI)
%  
% Required Helper Functions/Scripts:
%  - PLV 
%  - supplemental script for different surrogate generation methods

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

niter = 5000;
methods = {'iaaft', 'circshift', 'segment', 'shuffle'};
nbin = 60;
pb = linspace(-pi, pi, nbin);                                               % include -pi and pi for nicer plotting

%% a) Temporal Autocorrelation 
% compare autocorrelations for different surrogate generation methods
load([datapath 'zresp_003.mat']);
[ph,~] = two_point_interp(x);

% calculate AC for original signal
[tmpacf,~] = autocorr(x, NumLags = 600);                                    % max shift 6s                
empac = tmpacf(1:10:end);                                                   % downsampling for faster plotting

% calculate AC for all surrogate generation methods 
surrac = {};
for im = 1:length(methods)
    load([datapath 'surrogates_' methods{im} '_001.mat']);                   % k = 5000 surrogate phase time series 
    for iiter = 1:niter
        sresp = allsresp(iiter, :);
        [tmpacf,~] = autocorr(sresp,NumLags = 600);
        surrac{im}(:,iiter) = tmpacf(1:10:end);
    end 
end 

% plot
figure; 
hold on;
for im = 1:length(methods)
    plot(surrac{im}, 'Color', [colors(im, :, :) 0.05], 'LineWidth', 2);
end
plot(empac, 'Color', 'k', 'LineWidth',2);
xlim([0 60]);
ylim([-1 1]);
yticks(-1:0.5:1);
xticks(0:10:60);
xticklabels(string(0:6));                                                   % x-axis is lag in seconds
% legend created manually in power point 

%% b) Surrogate phase independence 
% surrogate generaton and phase-locking value calculation for surrogates
% see supp_phasegenmethods
% raincloud plots generated in R, see plots_fig3_rainclouds.R

%% c) Event-based statistics 
% hit rate ~ resp phase calculation for different surrogate generation
% procedures see supp_phasegenmethods
load([datapath 'surrts.mat']);

% main plot
figure;
hold on;
for im = 1:length(methods)
    plot(pb, surrts{im}, 'Color', colors(im, :, :), 'LineWidth', 4);
end
xlim([-pi pi]);
xticks([-pi 0 pi]);
xticklabels({'-\pi' '0' '\pi'});
yticks([-5 0 3]);

% negative peak inset
xlim([-0.8 -0.15]);
axis off;

%% d) Continuous statistics 
% Modulation Index spectrum calculation for empirical data and different 
% surrogate generation procedures see TutorialMI and supp_phasegenmethods
load([datapath 'freqvec.mat']);
load([datapath 'empmis.mat']);
load([datapath 'surrmis_.mat']);

% normalize the empirical MI values on the different surrogate distributions 
normmis = nan(20, length(f), length(methods));
for im = 1:length(methods)
    surr = surrmis{im};
    for isub = 1:20
        MI = empmis(isub, :);
        MI2 = squeeze(surr(isub, :, :))';
        bm = mean(MI2,1);
        bs = std(MI2,0,1);
        mi = (MI-bm)./bs;
        normmis(isub, :, im) = mi;
    end
end

% during MI calculation, spectra get sorted along descending frequencies as
% indicated in frequencies vector -> values are flipped manually to ensure
% correct plotting 
fflip = fliplr(f);
for im = 1:length(methods)
    flipmis(im, :) = fliplr(mean(squeeze(normmis(:, :, im)), 1));
end 

% plot for iaaft, circular shifting and segment shuffling
figure;
hold on;
for im = 1:length(methods)
    plot(fflip, flipmis(im, :), 'LineWidth', 4, 'Color', colors(im, :, :));
end
ylim([-1 1.5]);
yticks([-1 0 1]);
xlim([2 40]);
set(gca, 'FontSize', 27);

% random shuffling 
ylim([20 80]);
yticks([60]);
xticks([]);

