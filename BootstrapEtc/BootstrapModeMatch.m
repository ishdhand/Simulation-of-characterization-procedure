function [modeMatchPrecision] = BootstrapModeMatch(rateCalibBS,exptauvalsCalibBS,expTransCalib,modeMatchEstimated,...
    reflectivityCalibBS,landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVecCalibration,bootstrapSamples)
%BOOTSTRAPMODEMATCH bootstraps the mode match and reports its precision.
% inputs:
%       rateCalibBS: 1xn vector of the experimental coincidence
%           rates correspoding to the delays values in tauvals
%       exptauvalsCalibBS: 1xn vector of delay values from experiment.
%       expTransCalib: 1x3 vector of fitting parameters. The parameters
%           are rate scale, tau shift, tau scale.
%       modeMatchEstimated: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       reflectivityCalibBS: reflectivity of calibrating beamsplitter.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       sp: scaler integral scaling paramater.
%       scaleVecCalibration: 1x4 vector of optimization scaling parameters.
%       bootstrapSamples: number of times to bootstrap
% outputs:
%       modeMatchPrecision: precision of mode match estimate.


nCalibBS = length(exptauvalsCalibBS);

%% estimate precision of modeMatch

% generate beamsplitter matrix.
theta = asin(sqrt(reflectivityCalibBS));
amplitudeMatrixCalibBS = [cos(theta) sin(theta);sin(theta) cos(theta)];
phaseMatrixCalibBS = [0 pi; 0 0];

% refit and find residuals
rateCalibBSEstimated = calculateCoincidenceRateTransformed(exptauvalsCalibBS,expTransCalib,amplitudeMatrixCalibBS,...
    phaseMatrixCalibBS,modeMatchEstimated,landscapeCosFit,landscapeSinFit,background1,background2,sp);
 
residualsCalibBS = rateCalibBS-rateCalibBSEstimated;
normalresidualsCalibBS = residualsCalibBS./rateCalibBSEstimated;

% bootstrap modeMatch
modeMatchEstimatedBootstrap = zeros(1,bootstrapSamples);
beta0 = [modeMatchEstimated expTransCalib];
for counter=1:bootstrapSamples
    redrawnnormalresidualsCalibBS = randsample(normalresidualsCalibBS,nCalibBS,'true');
    redrawnresidualsCalibBS = redrawnnormalresidualsCalibBS.*rateCalibBSEstimated;
    rateCalibBSBootrap = rateCalibBSEstimated + redrawnresidualsCalibBS;
    % find mode match parameter
    modeMatchEstimatedBootstrap(counter) = Calibration(rateCalibBSBootrap,exptauvalsCalibBS,beta0,reflectivityCalibBS,...
        landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVecCalibration);

end
% estimate modeMatch precision
modeMatchPrecision = std(modeMatchEstimatedBootstrap);


end

