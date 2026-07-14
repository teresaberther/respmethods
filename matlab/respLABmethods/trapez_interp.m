function [phasevec, first, last] = trapez_interp(ts, fs)

% Trapez_interp extracts the phase of a oscillatory time series 
% (e.g. respiration) using the trapezium-area interpolation method.
% 
%
% Use as:
%   [phasevec, first, last] = trapez_interp(ts, fs)
%
% 
% Input arguments: 
%    ts = vector of zscored, possibly smoothed oscillatory time series (e.g. respiration trace) 
%    fs = sampling frequency of the time series
% 
% Output:  
%   phasevector = vector of phase value time series corresponding to the input timeseries
%   first       = first real phase value (method introduces NaNs at start of phase vector)
%   last        = last real phase value (method introduces NaNs at end of phase vector)
%
% Copyright (C) 2025, Teresa Berther, University of Münster, Germany 
% code adapted from workshop tutorial by Tahnée Engelen & Daniel Kluger

firstder = diff(ts);                                                        % get first derivative of time series
smfirstder = movmean(firstder, round(0.4 * fs));                            % smooth it 
[~, exp_begin_sample] = findpeaks(ts, 'MinPeakProminence',1);               % use a criterion of z = 1 for peak prominence

insp_begin_sample = [];
insp_begin_sample(1) = NaN;

for i_exps = 2:length(exp_begin_sample)

    cycle_length = exp_begin_sample(i_exps) - exp_begin_sample(i_exps-1);   % length of resp cycle in samples
    window_begin = round(exp_begin_sample(i_exps) - 0.7*cycle_length);      % start of search window (70% of length prior to exp begin)
    
    [~,xm] = max(smfirstder(window_begin:exp_begin_sample(i_exps)));        % find index of local max first derivative
    xm = xm + window_begin;                                                 % putting back the offset
    ym = ts(xm);                                                            % value of resp signal at point from where we start searching backwards
    xr = window_begin;                                                      % x val at which our search window begins
    yr = ts(xr);

    window_length = xm - xr;                                                % window length is used to determine nr of samples we fit trapeziums on

    xi = [];
    yi = [];
    trap_area = [];
    
    % iteratively calculating trapezium area while moving trap
    % towards start of our search window
    for isearch = 1:window_length-2
        xi(isearch) = xm - isearch;
        yi(isearch) = ts(xi(isearch));
        if yi(isearch) > ym
            break
        else
            trap_area(isearch) = abs(0.5.*(ym-yi(isearch)).*(2.*xr - xi(isearch) - xm));
        end

    end

    % noisy data catch
    if isempty(trap_area)
        insp_begin_sample(i_exps) = xr;
    else
        [~,max_area_id] = max(trap_area);
        insp_begin_sample(i_exps) = xm - max_area_id;
    end

end

% save cycles 
resp_cycles = [];
resp_cycle_counter = 1;
for irespcycles = 1:(length(insp_begin_sample)-1)
    cur_insp_sample = insp_begin_sample(irespcycles);
    next_insp_sample = insp_begin_sample(irespcycles+1);
    cur_exp_sample = exp_begin_sample(find(exp_begin_sample>cur_insp_sample & exp_begin_sample<next_insp_sample));
    % we only fill this structure with cycles that don't contain any NaNs 
    if numel(cur_exp_sample) == 1
        if ~isnan(cur_insp_sample) && ~isnan(next_insp_sample) && ~isnan(cur_exp_sample) 
            % filling the cycle sample information
            resp_cycles(resp_cycle_counter).beginSample = cur_insp_sample; 
            resp_cycles(resp_cycle_counter).expSample = cur_exp_sample;     % expiration start within this sample
            resp_cycles(resp_cycle_counter).endSample = next_insp_sample-1; % subtracting 1 sample so it's not the same value as start of next cycle
    
            resp_cycle_counter = resp_cycle_counter + 1;
        end
    end  
end


%%% phase assignment
phasevec = nan(1,length(ts));
resp_cycle_starts = [resp_cycles.beginSample];
resp_exp_starts   = [resp_cycles.expSample];
resp_cycle_ends   = [resp_cycles.endSample];

for ncycles_resp = 1:length(resp_cycles)

    Tar_bp = resp_cycle_starts(ncycles_resp);
    Tbr_bp = resp_exp_starts(ncycles_resp);
    Tcr_bp = resp_cycle_ends(ncycles_resp);
    t = Tar_bp:Tbr_bp;
    phasevec(Tar_bp:Tbr_bp) = -pi + pi*(t-Tar_bp)/(Tbr_bp-Tar_bp);
    t = Tbr_bp+1:Tcr_bp;
    phasevec(Tbr_bp+1:Tcr_bp) = pi*(t-Tbr_bp)/(Tcr_bp-Tbr_bp);

end

%%% find first and last valid phase
valid = find(~isnan(phasevec));
first = valid(1);
last  = valid(end);

end
