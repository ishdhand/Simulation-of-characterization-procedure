%% Calibration
% Estimates the mode match and the scaling and shifting of
% experimental data that matches th the theoretical curve for a given
% reflectivity beamsplitter.
% inputs:
%       rateCalibBS: 1xn vector of the experimental coincidence
%           rates correspoding to the delays values in tauvals
%       exptauvals: 1xn vector of delay values from experiment.
%       beta0: 1x4 vector of guesses for fitting parameters. The guesses
%           are modeMatch, rate scale, tau shift, tau scale.
%       reflectivityCalibBS: reflectivity of calibrating beamsplitter.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       sp: scaler integral scaling paramater.
%       scaleVec: 1x4 vector of optimization scaling parameters.
% outputs:
%       modeMatch: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       expTrans: 1x3 vector. The guesses are modeMatch, rate scale, tau
%           shift, tau scale.


function [modeMatch,expTrans] = Calibration(rateCalibBS,exptauvals,...
    beta0,reflectivityCalibBS,landscapeCosFit,landscapeSinFit,...
    background1,background2,sp,scaleVec)

% scale the parameters for optimization sake
beta0 = beta0./scaleVec;

theta = asin(sqrt(reflectivityCalibBS));
amplitudeMatrix = [cos(theta) sin(theta);sin(theta) cos(theta)];
phaseMatrix = [0 pi; 0 0];

% function to calculate the shifted and scaled coincidence rate given a
% modeMatch guess.
    function rate = modelfun(beta,taus)
        beta = beta.*scaleVec;
        modeMatchGuess = beta(1);
        rate = calculateCoincidenceRateTransformed(taus,beta(2:4),...
            amplitudeMatrix,phaseMatrix,modeMatchGuess,landscapeCosFit,...
            landscapeSinFit,background1,background2,sp);
        
    end

% weights to account for shot noise
weights = 1./abs(rateCalibBS);

options = statset('Display','off');
fittedparams = nlinfit(exptauvals,rateCalibBS,@modelfun,beta0,options,...
    'Weights',weights);

% scale params back
fittedparams = fittedparams.*scaleVec;

% outputs
modeMatch = fittedparams(1);
expTrans = fittedparams(2:4);

end

