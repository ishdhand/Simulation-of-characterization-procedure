%% Calculate coincidence probability 
% This function computes the expected coincidence probability for a given 
% transition matrix, multiphoton rate, filter functions and relative delay
% values. SPDC type of light source is assumed which either emits 1 or 2 
% pairs of photons in two channels. The setup is assumed to a heraled one
% for now.
%
%% Inputs:
%       amplitudeMatrix: 2x2 matrix containing the amplitudes of the transition matrix.
%       phaseMatrix: 2x2 matrix containting the phases of the transition matrix.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%       eta: positive real number between 0 and 1 inclusive. Probability of multiphoton events.
%       modeMatch: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%       tauvals: 1xn values of relative delay for which the rate is computed
%
%% Outputs
%       prob: 1xn values of expected coincidence rates.
%
%% Procedure

function prob = calculateCoincidenceProbability(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,eta,modeMatch,tauvals)
% find scaling parameter
sp = scalingParameter(amplitudeMatrix,freq,phiv1,phiv2,fv1,fv2);

% scale spectra
fv1 = fv1*sp;
fv2 = fv2*sp;
phiv1 = phiv1*sp;
phiv2 = phiv2*sp;

% find coincidence probability due to a single pair
prob11 = calculateCoincidenceProbability112(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)/sp^4;

% if there are multiphotons
if eta>0
    % compute lowest order multiphoton probabilities
    prob21 = calculateCoincidenceProbability21(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)/sp^6;
    prob12 = calculateCoincidenceProbability12(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)/sp^6;
    prob = prob11+eta^2*(prob21+prob12);
else
    prob = prob11;
end



end