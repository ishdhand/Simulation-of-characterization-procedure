% Calculate experimental coincidence rate
% Computes the expected coincidence
% rate assuming a given transition matrix, multiphoton rate, filter
% functions and relative delay values. SPDC type of light source is assumed
% which either emits 1 or 2 pairs of photons in two channels. The setup is
% assumed to a heraled one for now. The coincidence probability is shifted and scaled.

%% Inputs:
%       exptauvals: 1xn values of relative delay for which the rate is computed
%       expTrans: 1x3 vector. scale and shift of experimental data.
%       amplitudeMatrix: 2x2 matrix containing the amplitudes of the transition matrix.
%       phaseMatrix: 2x2 matrix containting the phases of the transition matrix.
%       modeMatch: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%       eta: positive real number between 0 and 1 inclusive. Probability of multiphoton events.
%
%% Outputs
%       exprate: 1xn values of expected coincidence rates in an experiment

%% Procedure


function exprate = calculateExperimentalCoincidenceRate(exptauvals,expTrans,amplitudeMatrix,phaseMatrix,modeMatch,...
    freq,phiv1,phiv2,fv1,fv2,eta)

c = 3e8; % speed of light

% shift and scale tau values. There is an implicit factor of the speed of
% light.
tauvals = (exptauvals - expTrans(1))/expTrans(2)/c;
% calculate the coincidence probabily with these values.
prob = calculateCoincidenceProbability(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,eta,modeMatch,tauvals);
% scale and shift the prob to get the experimental rate.
exprate = prob*expTrans(3);

end

