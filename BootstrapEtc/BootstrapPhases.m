function phasesMatrixPrecision = BootstrapPhases(ratesInterferometer,exptauvalsInterferometer,expTransArray,modeMatchEstimated,...
    relativeAmplitudesMatrix,phasesMatrix,freq,phiv1,phiv2,fv1,fv2,...
    landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVecPhase,tauScaling,bootstrapSamples)
%BOOTSTRAPPHASES Bootstraps the phases and estimates their precision
% inputs:
%       ratesInterferometer:mxm cell array with each cell holding
%           the coincidence rates between input channels 1 and h and output
%           channels 1 and g. The first row and column are empty and not
%           used.
%       exptauvalsInterferometer: mxm cell array with each cell holding the
%           a 1xn vector of delay values for the corresponding HOM dip.
%       expTransArray: mxm cell array of estimated curve fitting
%           parameters. The first row and col are empty and not used.
%       modeMatchEstimated: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%       relativeAmplitudesMatrix: mxm matrix having the amplitudes of the
%           transition matrix relative to the entries of the first row and
%           column. First row and col are 1.
%       phasesMatrix: mxm matrix with the phases of the
%           transition matrix. The first row and column are 0.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%       eta: positive real number between 0 and 1 inclusive. Probability of
%           multiphoton events.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       sp: scaler integral scaling paramater.
%       scaleVecPhase: 1x4 vector of optimization scaling parameters.
%       tauScaling: parameter that scales tau.
%       bootstrapSamples: number of times to bootstrap
% outputs:
%       phasesMatrixPrecision: mxm matrix which has the precision of the
%       phases of the reconstructed matrix. The first row and col are 0.

m = length(relativeAmplitudesMatrix);

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
        expTrans = [expTransArray{g,h} tauScaling.*scaleVecPhase(4)];
        coincidenceRatesEstimated{g,h} = calculateCoincidenceRateTransformed(exptauvalsInterferometer{g,h,1,1},expTrans,...
            amplitudeMatrixgh,phaseMatrixgh,modeMatchEstimated,landscapeCosFit,landscapeSinFit,background1,background2,sp);
        residualsArray{g,h} = ratesInterferometer{g,h,1,1} - coincidenceRatesEstimated{g,h};
        normalresidualsArray{g,h} = residualsArray{g,h}./coincidenceRatesEstimated{g,h};
    end
end


% bootstrap phases
ratesInterferometerBootrap = cell(m,m);
phasesMatrixBootstrap = cell(1,bootstrapSamples);
for counter=1:bootstrapSamples
    % redraw residuals
    for g=2:m
        for h=2:m
            redrawnnormalresiduals = randsample(normalresidualsArray{g,h},length(exptauvalsInterferometer{g,h,1,1}),'true');
            redrawnresiduals = redrawnnormalresiduals.*coincidenceRatesEstimated{g,h};
            ratesInterferometerBootrap{g,h} = coincidenceRatesEstimated{g,h} + redrawnresiduals;

        end
    end
    % run phase estimation on entire interferometer.
    phasesMatrixBootstrap{counter} = PhaseEstimation(ratesInterferometerBootrap,exptauvalsInterferometer,...
        relativeAmplitudesMatrix,modeMatchEstimated,freq,phiv1,phiv2,fv1,fv2,...
        landscapeCosFit,landscapeSinFit,background1,background2,scaleVecPhase,tauScaling);

end

% std of the estimated phases is the error
phasesMatrixPrecision = zeros(m,m);
phaseghBootstrap = zeros(1,bootstrapSamples);
for g=2:m
    for h=2:m
        % extract the bootstrap phases corresponding to g,h
        for counter=1:bootstrapSamples
            phaseMatrixBootstrap = phasesMatrixBootstrap{counter};
            phaseghBootstrap(counter) = phaseMatrixBootstrap(g,h);
        end
        % the precision
        phasesMatrixPrecision(g,h) = std(phaseghBootstrap);
    end
end



end

