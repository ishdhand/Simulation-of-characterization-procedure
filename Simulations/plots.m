
% figure;
% hold on;
k=1; %input
j=1; %output
for g=2:msize
    for h=2:msize
        plot(simulatedtauvals{g,h,j,k},...
            simulatedRatesInterferometer{g,h,j,k},'.');
        hold on;
        amplitudeMatrix = relativeAmplitudesMatrix([j g],[k h]);
        phaseMatrix = [0 0; 0 phasesMatrix(g,h)];
        r = calculateCoincidenceRateTransformed(...
            simulatedtauvals{g,h,j,k},[expTransArray{g,h,j,k}...
            expTransCalib(3)],amplitudeMatrix,phaseMatrix,...
            modeMatchEstimated,landscapeCosFit,landscapeSinFit,...
            background1,background2,sp);
        plot(simulatedtauvals{g,h,j,k},r)
    end
end