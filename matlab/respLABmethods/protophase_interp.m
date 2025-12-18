function [phasevec, first, last] = protophase_interp(ts)

% Protophase_interp extracts the phase of a oscillatory time series 
% (e.g. respiration) using the protophase interpolation method described by 
% Rosenblum, M. & Pikovsky, A. in e.g. "Inferring connectivity of an oscillatory 
% network via the phase dynamics reconstruction". Front. Netw. Physiol. 3, 1298228 (2023).
%
% Use as:
%   [phasevec, first, last] = protophase_interp(ts)
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
% This code is a simpified version of the co_mmzproto and co_fbtrT functions of the DAMOCO
% toolbox provided by Kralemann, B., Rosenblum, M., and Pikovsky, A., 2014:
% http://www.stat.physik.uni-potsdam.de/~mros/damoco2.html
%
% Please cite one of these papers when the provided function is used:
%
% Rosenblum, M. & Pikovsky, A. Inferring connectivity of an oscillatory 
% network via the phase dynamics reconstruction. Front. Netw. Physiol. 3, 1298228 (2023).
%
% B. Kralemann, M. Frühwirth, A. Pikovsky, M. Rosenblum,  T. Kenner,  J. Schaefer, and M. Moser
% In vivo cardiac phase response curve elucidates human respiratory heart rate variability, 
% Nature Communications, 4, 2418 (2013).

% Adapted by Daniel Kluger & Teresa Berther, University of Münster, Germany (2025)

%%% use findpeaks to get peak/trough indices
ts = -ts;                                                                   % flip zscored resp trace
[~, IN1] = findpeaks(ts, 'MinPeakProminence',1);                            % use a criterion of z = 1 for peak prominence
[~, IN3] = findpeaks(-ts, 'MinPeakProminence',1);                           % use a criterion of z = 1 for peak prominence
IN2 = find(ts(1:end-1)>0 & ts(2:end)<0);                                    % indices of the zero crossings: positive to negative
IN4 = find(ts(1:end-1)<0 & ts(2:end)>0);                                    % indices of the zero crossings: negative to positive

for n = 1:length(IN1)-1
    D(n) = IN1(n+1)-IN1(n);                                                 % length of each cycle
    R2(n) = (IN2(n)-IN1(n)) / D(n);                                         % position of the 1st zero crossing relative to the cycle length
    R3(n) = (IN3(n)-IN1(n)) / D(n);                                         % position of the minimum relatively to the cycle length
    R4(n) = (IN4(n)-IN1(n)) / D(n);                                         % position of the 2nd zero crossing relative to the cycle length
end

R2 = mean(R2);                                                              % average position of the 1st zero crossing
R3 = mean(R3);                                                              % average position of the minima
R4 = mean(R4);                                                              % average position of the 2nd zero crossings

IN = [];
Pin = [];
for n = 1: length(IN1)-1
    IN = [IN IN1(n) IN2(n) IN3(n) IN4(n)];                                  % sort the positions of markers for interpolation
    Pin = [Pin 2*pi*(n-1)  2*pi*(n-1+R2) 2*pi*(n-1+R3) 2*pi*(n-1+R4)];      % sort the values of the protophase for corresponding markers
end 

IN = [IN IN1(end)];
Pin = [Pin 2*pi*(length(IN1)-1)];
protophase = interp1(IN, Pin, (IN1(1)+1:1:IN1(end)),'linear');              % compute the protophase using linear interpolation between markers 

first = IN1(1)+1;                                                           % get index of first phase value so we can exclude NaNs later
last = IN1(end);                                                            % get index of last phase value so we can exclude NaNs later

protophase = mod(protophase,2*pi);

nfft = 100;                                                                 % maximal number of Fourier coefficients
Spl = zeros(nfft,1);                                                        % fourier coefficients 1,...,nfft
Hl = zeros(nfft,1);                                                         % Tenreiro function to be minimized

IN = find(diff(mod(protophase,2*pi))<0);                                    % correction for short time series:
npt = length(protophase(IN(1) : IN(end)));                                  % only full periods are used

S = 0; 
c = double(npt+1)/double(npt-1);
for k = 1:nfft                                                              % computing Fourier coefficients
    Spl(k) = sum(exp(-1i*k*protophase(IN(1):IN(end))))/npt;
    S = S+Spl(k)*conj(Spl(k))-1./double(npt);
    Hl(k) = k/npt-c*S;                                                      % Tenreiro function
end
[~,indopt] = min(Hl);

phasevec = protophase;                                                      % transformation protophase --> phasevec
for k = 1:indopt
    phasevec = phasevec+2*imag(Spl(k)*(exp(1i*k*protophase)-1)/k);
end

end 
