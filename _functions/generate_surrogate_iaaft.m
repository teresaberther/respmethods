function surrogateData = generate_surrogate_iaaft(originalData, varargin)
% GENERATE_SURROGATE_IAAFT Generates surrogates of a time series.
%
%	surrogateData =  GENERATE_SURROGATE_IAAFT(originalData, ...)
%	Generates, by means of the Iterative Amplitude Adjusted Fourier Transform
%	(IAAFT) algorithm, surrogate time series corresponding to the original
%	time series originalData.
%
%	Options, passed as ('name',value) pairs:
%
%	'M'		Number of surrogate time series to generate.
%			By default, M = 1. Each of the M surrogate time series
%			is a column in the returned array, which therefore
%			has a size N x M, where N is the original time series
%			length.
%
%	'detrend'	Specifies whether the time series has to be detrended
%			prior to surrogate generation. Default is false.
%
%	'verbose'	Sets the verbosity of the function. If true (default),
%			all messages are printed on the command line. If
%			false, only critical errors are displayed.
%
%	NSE Laboratory, Dept. of Physics, University of Trento
%
%	REFERENCE:
%	The IAAFT (Iterative Amplitude-Adjusted Fourier Transform) algorithm
%	for surrogate generation was originally proposed by T. Schreiber
%	and A. Schmitz in Phys. Rev. Lett. 77 (1996), 635,
%	<a href="matlab:web('https://doi.org/10.1103/PhysRevLett.77.635')">doi:10.1103/PhysRevLett.77.635</a>

	% --- Input parsing & validation
	ip = inputParser;
	addRequired(ip, 'originalData', @isvector);
	addParameter(ip, 'M', 1, @isnumeric);
	addParameter(ip, 'detrend', true, @islogical);
	addParameter(ip, 'verbose', true, @islogical);
	ip.KeepUnmatched = false;
	parse(ip, originalData, varargin{:});
	if (~iscolumn(originalData))
		originalData = originalData';
	end
	N = length(originalData);
	M = fix(ip.Results.M);
	if (M < 1)
		error('Function argument "M" must be a positive integer.')
	end
	beVerbose = ip.Results.verbose;
	doDetrend = ip.Results.detrend;

	% --- Provide some information
	if (beVerbose)
		fprintf('\n### Starting IAAFT routine ###\n');
	end

	% --- Detrend, if requested
	if (doDetrend)
		x0 = originalData(1);
		xN = originalData(N);
		if (beVerbose)
			fprintf('\n### Detrending time series ###\n');
		end
		timeIdx = [0:1:N-1]';
		originalData = originalData - (xN-x0)*timeIdx/(N-1.0);
	end

	% --- Assess spectrum & distribution of original sequence
	spectrumMagnitudeOriginal = abs(fft(originalData));
	distributionOriginal = sort(originalData);

	% --- Initialize output array (each of the M generated surrogates is a column)
	surrogateData = [];

	% --- Generate surrogates: iterate M times
	for iter=1:M
		if(beVerbose)
			fprintf('# Surrogate number %d out of %d.\n', iter, M);
		end

		% --- Starting conditions
		rng('shuffle', 'twister');
		runIterations = true;
		thisSurrogate(:,1) = originalData(randperm(N));
		lastIterSurrogate = thisSurrogate;

		% --- Iterative algorithm (switched depending on target function to be matched exactly.
		nIter = 0;
		while (runIterations)

			nIter = nIter + 1;

			spectrumSurrogate = fft(thisSurrogate);
			phasesSurrogate = angle(spectrumSurrogate);
			spectrumSurrogate = spectrumMagnitudeOriginal .* (cos(phasesSurrogate) + 1i.*sin(phasesSurrogate));
			thisSurrogate = ifft(spectrumSurrogate, 'symmetric');

			[~, iSorted] = sort(thisSurrogate);
			thisSurrogate(iSorted) = distributionOriginal;

			if (iterationConverged())
				break;
			end

			lastIterSurrogate = thisSurrogate;
		end

		if (beVerbose)
			fprintf('Process converged after %d iterations.\n\n', nIter);
		end

		% --- Append sequence to the set of those already generated
		surrogateData = [surrogateData, thisSurrogate];
	end

	% --- Define check of convergence (sequence did not change with respect to last iteration)
	function convergence = iterationConverged()
		z = sum((lastIterSurrogate - thisSurrogate).^2);
		I = sum(thisSurrogate.^2);
		if (z / I <= 1e-6)
			convergence = true;
		else
			convergence = false;
		end
	end

end