%% Generate Gaussuan Spectrum (generateGaussianSpectrum)
% Generates normalized Gaussain spectral functions  in the frequency domain
% when spreads of wavelengths and central wavelengths are given.

%% Inputs
% lambda: 1xq vector, each element is a central wavelengths (in nm).
% sigma : 1xq vector, each element is a standard deviations (in nm).
% p: size of each of the q output function.

%% Outputs
% freq: px1 vector of frequencies at which the spectral function
% varargout: q vectors of amplitude at the respective frequncy. Each vector
% has size px1.

%% Procedure 

function [freq,varargout] = generateGaussianSpectrum(lambda,sigma,p)

c = 3e8;
% Speed of light. Please do not change!
nsigmas = 2;
% How many widths of the Gaussian are captured in the generated output.

lambda = lambda*10^(-9);
sigmaLambda = sigma*10^(-9);
% Convert to SI Units

omega = 2*pi*c./lambda;
sigmaFreq = 2*pi*c*sigmaLambda./lambda.^2;
% Mean and standard deviation of the frequency spectral function.

freq = linspace(min(omega)-nsigmas*max(sigmaFreq),max(omega)+...
    nsigmas*max(sigmaFreq),p);
% Frequncy range is uniformly spaced between highest and lowest values.

varargout = cell(1,length(lambda));
for i=1:length(lambda)
    varargout{i} = GaussianFunction(freq,omega(i),sigmaFreq(i));
end
% Generate Gaussian spectrum in the frequency

end

%% Gaussian Function
% Finds value at w of Gaussian function with mean w0 and std. s
function amp = GaussianFunction(w,w0,s)
amp = exp(-(w-w0).^2/(4*s^2))/sqrt(s*sqrt(2*pi));
end
