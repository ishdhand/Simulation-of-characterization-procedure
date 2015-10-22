% Phase estimation two channel
% Estimates the magnitude of one unknown phase of
% the transition matrix of a 2 channel interferometer using two photon
% statistics.
%% Inputs
%       rateEstimatingInterferometer: 1xn vector of the experimental
%           coincidence rates correspoding to the delays values in tauvals
%       exptauvals: 1xn vector of delay values from experiment.
%       beta0: 1x4 vector of guesses for fitting parameters. First three
%           are scale and shift of the experimental data, the fourth the
%           phaseMatrix(2,2) guess.
%       amplitudeMatrix: 2x2 matrix that has the amplitudes of the
%           interferometer. The amplitudes may only be correct up to losses
%           at the inputs and outputs of the interferometer.
%       phaseMatrix: 2x2 matrix containting the phases of the transition
%           matrix. The last entry is to be estimated and it's value is not
%           used.
%       modeMatch: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       sp: scaler integral scaling paramater.
%       scaleVec: 1x3 vector of optimization scaling parameters.
%       tauScaling: scaling of the delay values.
%% Outputs
%       phaseMag: the magnitude of the phase.
%       expTrans: 1x2 vector. scale and shift of experimental data.
%% Procedure

function [phaseMag, expTrans] = PhaseEstimation2Channel(rateEstimatingInterferometer,exptauvals,beta0,...
    amplitudeMatrix,phaseMatrix,modeMatch,landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVec,tauScaling)

% scale the parameters for optimization sake
beta0 = beta0./scaleVec(1:3);

% function to calculate the shifted and scaled coincidence rate given a
% phase guess. 
    function rate = modelfun(beta,taus)
        beta(4) = tauScaling;
        beta = beta.*scaleVec;
        phaseMatrix(2,2) = beta(1);
        rate = calculateCoincidenceRateTransformed(taus,beta(2:4),amplitudeMatrix,phaseMatrix,modeMatch,landscapeCosFit,landscapeSinFit,background1,background2,sp);
    end

% weights to account for shot noise
weights = 1./abs(rateEstimatingInterferometer);
        
% optimize
% options = optimset('Display','off');
%fittedparams = nlinfit(exptauvals,rateEstimatingInterferometer,@modelfun,beta0,options,'Weights',weights);
fittedparams = nlinfit(exptauvals,rateEstimatingInterferometer,@modelfun,beta0,'Weights',weights);

% scale params back
fittedparams = fittedparams.*scaleVec(1:3);

% the magnitude of the phase
phaseMag = fittedparams(1);

% scale and shift params
expTrans = fittedparams(2:3);

end
