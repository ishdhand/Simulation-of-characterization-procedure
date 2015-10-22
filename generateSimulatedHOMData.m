% Generate Simulated HOM data
% Generates simulated hom data
%% Inputs:
%       tauvals: 1xn vector of delay values.
%       invExpTrans: 1x3 vector. They first two scale and shift the time
%           while the last scales the rate.
%       amplitudeMatrix: 2x2 matrix that has the amplitudes of the
%           interferometer.
%       phaseMatrix: 2x2 matrix containting the phases of the transition
%           matrix.
%       modeMatch: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%       eta: positive real number between 0 and 1 inclusive. Probability of
%           multiphoton events.
%% Outputs:
%       exptauvals: 1xn vector of scaled and shifted delay values.
%       exprate: 1xn vector of scaled exprate.
%% Procedure

function [exptauvals,exprate] = generateSimulatedHOMData(tauvals,...
    invExpTrans,amplitudeMatrix,phaseMatrix,modeMatch,...
    freq,phiv1,phiv2,fv1,fv2,eta)

% compute the coincidence probability
prob = calculateCoincidenceProbability(amplitudeMatrix,phaseMatrix,freq,...
    phiv1,phiv2,fv1,fv2,eta,modeMatch,tauvals);

% the delays and prob are transformed to measured quantities in experiment
exptauvals = invExpTrans(3)*tauvals + invExpTrans(2);

% scale the rate so that number of clicks at wings are equal to scale
% parameter
% scaledrate = invExpTrans(3)*prob;
S = invExpTrans(1)/prob(1);
scaledrate = S*prob;

% add shot noise
exprate = poissrnd(scaledrate);
exprate(exprate < 0) = 0;

end

