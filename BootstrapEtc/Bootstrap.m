function [phasesMatrixErrors,modeMatchError] = Bootstrap(rateCalibBS,exptauvalsCalibBS,modeMatchEstimated,expTransCalib,...
    ratesInterferometer,exptauvalsInterferometer,relativeAmplitudesMatrix,reflectivityCalibBS,...
    freq,phiv1,phiv2,fv1,fv2,bootstrapSamples)
%BOOTSTRAP 

m = length(amplitudesTransitionMatrix);
nCalibBS = length(exptauvalsCalibBS);

%% estimate precision of modeMatch

% generate beamsplitter matrix.
theta = asin(sqrt(reflectivityCalibBS));
amplitudeMatrixCalibBS = [cos(theta) sin(theta);sin(theta) cos(theta)];
phaseMatrixCalibBS = [0 pi; 0 0];

% refit and find residuals
rateCalibBSEstimated = calculateCoincidenceRateTransformed(exptauvalsCalibBS,expTransCalib,amplitudeMatrixCalibBS,phaseMatrixCalibBS,modeMatchEstimated,landscapeCosFit,landscapeSinFit,background1,background2,sp);
 
residualsCalibBS = rateCalibBS-rateCalibBSEstimated;
normalresidualsCalibBS = residualsCalibBS./rateCalibBSEstimated;

% bootstrap modeMatch
modeMatchEstimatedBootstrap = zeros(1,boostrapSamples);
for counter=1:bootstrapSamples
    redrawnnormalresidualsCalibBS = randsample(normalresidualsCalibBS,nCalibBS,'true');
    redrawnresidualsCalibBS = redrawnnormalresidualsCalibBS.*rateCalibBSEstimated;
    rateCalibBSBootrap = rateCalibBSEstimated + redrawnresidualsCalibBS;
    % find mode match parameter
    modeMatchEstimatedBootstrap(counter) = Calibration(rateCalibBSBootrap,exptauvalsCalibBS,expTransCalib,reflectivityCalibBS,landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVecCalibration);

end
% estimate modeMatch precision
modeMatchError = std(modeMatchEstimatedBootstrap);

%% Estimate precision of phases

% find phases and residuals
amplitudeMatrixgh = ones(2,2);
phaseMatrixgh = zeros(2,2);
coincidenceRatesEstimated = cell(m,m);
residualsArray = cell(m,m);
normalresidualsArray = cell(m,m);

% refit and find residuals for all phases
for g = 2:m
    for h = 2:m
        
        amplitudeMatrixgh(2,2) = relativeAmplitudesMatrix(g,h);
        phaseMatrixgh(2,2) = phasesMatrix(g,h);
        
        coincidenceRatesEstimated{g,h} = calculateCoincidenceRateTransformed(exptauvalsInterferometer{g,h},expTransArray{g,h},amplitudeMatrixgh,phaseMatrixgh,modeMatchEstimated,landscapeCosFit,landscapeSinFit,background1,background2,sp);
        residualsArray{g,h} = ratesInterferometer{g,h} - coincidenceRatesEstimated{g,h};
        normalresidualsArray{g,h} = residualsArray(g,h)./coincidenceRatesEstimated{g,h};
    end
end


% bootstrap phases
ratesInterferometerBootrap = cell(m,m);
phasesMatrixBootstrap = cell(bootstrapSamples,m,m);
for counter=1:bootstrapSamples
    % redraw residuals
    for g=2:m
        for h=2:m
            redrawnnormalresiduals = randsample(normalresidualsArray{g,h},length(exptauvalsInterferometer{g,h}),'true');
            redrawnresiduals = redrawnnormalresiduals.*coincidenceRatesEstimated{g,h};
            ratesInterferometerBootrap{g,h} = coincidenceRatesEstimated{g,h} + redrawnresiduals;

        end
    end
    % run phase estimation on entire interferometer.
    phasesMatrixBootstrap{counter,g,h} = PhaseEstimation2(ratesInterferometerBootrap,exptauvalsInterferometer,relativeAmplitudesMatrix,modeMatchEstimated,freq,phiv1,phiv2,fv1,fv2,landscapeCosFit,landscapeSinFit,background1,background2,scaleVecPhase,tauScaling);

end

phasesMatrixErrors = zeros(m,m);
phaseghBootstrap = zeros(1,bootstrapSamples);
for g=2:m
    for h=2:m
        % std of the estimated phases is the error
        for counter=1:bootstrapSamples
            phaseghBootstrap(counter) = phasesMatrixBootstrap{counter,g,h};
            phasesMatrixErrors(g,h) = std(phaseghBootstrap);
        end
    end
end


end

