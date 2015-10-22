%% Amplitude Estimation
% Compute the amplitudes of the transition matrix. Uses single photon 
% statistics to compute the amplitudes of the transition matrix relative
% to the elements of the first row and column 

%% Inputs.
% singlePhotonClicks: m x m matrix where the (i,j) entry is the number of 
% clicks in output channel j when single photons are inpinged at input 
% channel i.

%% Outputs
% relativeAmplitudes: m x m matrix which records the amplitude of the
%           transition matrix relative to the elements of the first row and
%           column. The entries of the first row and column are 1.
% relativeAmplitudeErrors: m x m matrix which holds the estimated
%           shot noise error in each entry. The first row and column are 0.

%% Procedure

function [relativeAmplitudes,relativeAmplitudeErrors] = ...
		AmplitudeEstimation(singlePhotonClicks)

% size of interferometer
m = length(singlePhotonClicks);


relativeAmplitudes = ones(m,m);
relativeAmplitudeErrors = zeros(m,m);


for g=2:m
		for h=2:m
				% the relative amplitudes are given by the ratio of clicks
				relativeAmplitudes(g,h) = sqrt(singlePhotonClicks(1,1)*...
						singlePhotonClicks(g,h)/(singlePhotonClicks(1,h)*...
						singlePhotonClicks(g,1)));
				% the errors are given by the adding the variance
				relativeAmplitudeErrors(g,h) = 0.5*...
						(1/sqrt(singlePhotonClicks(1,1)) + ...
						1/sqrt(singlePhotonClicks(g,h)) + ...
						1/sqrt(singlePhotonClicks(1,h)) + ...
						1/sqrt(singlePhotonClicks(g,1)));
		end
end

end

