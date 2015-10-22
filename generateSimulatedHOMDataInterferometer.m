% Generate Simulated HOM Data Interferometer 
% Generates simulated hom data for (m-1)^2 different experiments, where the
% first input and output channel is fixed.
%% Inputs:
%       tauvals: 1xn vector of delay values.
%       invExpTrans: 1x3 vector. They first two scale and shift the time
%           while the last scales the rate.
%       interferometerMatrix: mxm unitary matrix.
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
%       simulatedtauvals: mxm cell array with each cell holding the a 1xn vector
%           of delay values for the corresponding HOM dip.
%       simulatedRatesInterferometer: mxm cell array with each cell holding
%           the coincidence rates between input channels 1 and h and output
%           channels 1 and g. The first row and column are empty and not
%           used.
%% Procedure
function [simulatedtauvals,simulatedRatesInterferometer] = generateSimulatedHOMDataInterferometer(tauvals,invExpTrans,...
    interferometerMatrix,modeMatch,freq,phiv1,phiv2,fv1,fv2,eta)

m = length(interferometerMatrix);

% seperate the interferometer into amplitudes and phases
amplitudesMatrixInterferometer = abs(interferometerMatrix);
phasesMatrixInterferometer = angle(interferometerMatrix);

simulatedRatesInterferometer = cell(m,m,m,m);
simulatedtauvals = cell(m,m,m,m);

% simulate HOM data assuming that one input and one output is always in
% channel 1. j=1,k=1
k=1; %input
j=1; %output
for g=2:m
    for h=2:m
        amplitudeMatrix = amplitudesMatrixInterferometer([j g],[k h]);
        phaseMatrix = phasesMatrixInterferometer([j g],[k h]);
        [simulatedtauvals{g,h,j,k},simulatedRatesInterferometer{g,h,j,k}] = generateSimulatedHOMData(tauvals,invExpTrans,amplitudeMatrix,phaseMatrix,modeMatch,freq,phiv1,phiv2,fv1,fv2,eta);
    end
end

% simulate HOM data assuming that the inputs are in channel 1 and 2 and one
% output is in channel 2.  
k=1; %input
h=2; %input
j=2; %output
for g=3:m
    amplitudeMatrix = amplitudesMatrixInterferometer([j g],[k h]);
    phaseMatrix = phasesMatrixInterferometer([j g],[k h]);
    [simulatedtauvals{g,h,j,k},simulatedRatesInterferometer{g,h,j,k}] = generateSimulatedHOMData(tauvals,invExpTrans,amplitudeMatrix,phaseMatrix,modeMatch,freq,phiv1,phiv2,fv1,fv2,eta);
end


% simulate HOM data assuming that one input is in channel 2 and both
% outputs are in channel 1 and 2.
k=2; %input
j=1; %output
g=2; %output
for h=2:m
    amplitudeMatrix = amplitudesMatrixInterferometer([j g],[k h]);
    phaseMatrix = phasesMatrixInterferometer([j g],[k h]);
    [simulatedtauvals{g,h,j,k},simulatedRatesInterferometer{g,h,j,k}] = generateSimulatedHOMData(tauvals,invExpTrans,amplitudeMatrix,phaseMatrix,modeMatch,freq,phiv1,phiv2,fv1,fv2,eta);
end

% simulate HOM data assuming that one input is in channel 2 and one output
% is in channel 2.
k=2; %input
j=2; %output
for g=3:m
    for h=3:m
        amplitudeMatrix = amplitudesMatrixInterferometer([j g],[k h]);
        phaseMatrix = phasesMatrixInterferometer([j g],[k h]);
        [simulatedtauvals{g,h,j,k},simulatedRatesInterferometer{g,h,j,k}] = generateSimulatedHOMData(tauvals,invExpTrans,amplitudeMatrix,phaseMatrix,modeMatch,freq,phiv1,phiv2,fv1,fv2,eta);
    end
end

end






