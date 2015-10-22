%% Phase estimation 
% Computes the phases of the transition matrix assuming
% known coincidence rates between input/output port 1 and all other input
% output port combinations, and the amplitudes of the transition matrix.
%% Inputs:
%       ratesInterferometer:mxm cell array with each cell holding
%           the coincidence rates between input channels 1 and h and output
%           channels 1 and g. The first row and column are empty and not
%           used.
%       exptauvals: mxm cell array with each cell holding the a 1xn vector
%           of delay values for the corresponding HOM dip.
%       relativeAmplitudesMatrix: mxm matrix having the amplitudes of the
%           transition matrix relative to the entries of the first row and
%           column. First row and col are 1.
%       modeMatch: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       scaleVecPhase: 1x4 vector of optimization scaling parameters.
%       tauScaling: parameter that scales tau.
%% Outputs:
%       phasesMatrix: mxm matrix with the phases of the
%           transition matrix. The first row and column are 0.
%       expTransArray: mxm cell array of estimated curve fitting parameters. The
%       first row and col are empty and not used.
%% Procedure

function [phasesMatrix,expTransArray] = PhaseEstimation(ratesInterferometer,exptauvals,relativeAmplitudesMatrix,modeMatch,...
    freq,phiv1,phiv2,fv1,fv2,landscapeCosFit,landscapeSinFit,background1,background2,scaleVecPhase,tauScaling)

% size of interferometer
m = length(relativeAmplitudesMatrix);
% phasesMatrixMags = zeros(m,m);
% phasesMatrixSigns = ones(m,m);
expTransArray = cell(m,m);
betaAbs = zeros(m,m);
phasesMatrix = zeros(m,m);

phaseMatrix = zeros(2,2);

    function M = PhaseEstimation(M,gVals,hVals,jInd,kInd)
        
        for gInd = gVals
            for hInd = hVals
                        rategh = ratesInterferometer{gInd,hInd,jInd,kInd};
                        exptauvalsgh = exptauvals{gInd,hInd,jInd,kInd};
                        amplitudeMatrix = relativeAmplitudesMatrix([jInd gInd],[kInd hInd]);
                        sp = scalingParameter(amplitudeMatrix,freq,phiv1,phiv2,fv1,fv2);
                        % get beta
                        beta0 = HOMFittingGuess(exptauvalsgh,rategh,amplitudeMatrix,'phase',phaseMatrix,modeMatch);
                        [M(gInd,hInd), expTransArray{gInd,hInd}] = PhaseEstimation2Channel(rategh,exptauvalsgh,beta0,amplitudeMatrix,...
                            phaseMatrix,modeMatch,landscapeCosFit,landscapeSinFit,background1,background2,sp,scaleVecPhase,tauScaling);
            end
        end
    end


% compute the magnitude of each of the phases.
% k=1; %input
% j=1; %output
phasesMatrix = PhaseEstimation(phasesMatrix,2:m,2:m,1,1);

% only the magnitude of the recovered values is relevant.
phasesMatrix = abs(phasesMatrix);

% k=1; %input
% h=2; %input
% j=2; %output
betaAbs = PhaseEstimation(betaAbs,3:m,2,2,1);

% k=2; %input
% j=1; %output
% g=2; %output
betaAbs = PhaseEstimation(betaAbs,2,3:m,1,2);


% k=2; %input
% j=2; %output
betaAbs = PhaseEstimation(betaAbs,3:m,3:m,2,2);

% only the absolute value of beta matters
betaAbs = abs(betaAbs);

% 2,2 has positive sign. find signs of other entries

% find signs of second column
k=1; %input
h=2; %input
j=2; %output
for g = 3:m
    phasesMatrix(g,h) = computeSigns(phasesMatrix([j g],[k h]),betaAbs(g,h))*phasesMatrix(g,h);
end

%find signs of second row
k=2; %input
j=1; %output
g=2; %output
for h = 3:m
    phasesMatrix(g,h) = computeSigns(phasesMatrix([j g],[k h]),betaAbs(g,h))*phasesMatrix(g,h);
end

% find signs of all other entries
k=2;
j=2;
for g = 3:m
    for h = 3:m
        phasesMatrix(g,h) = computeSigns(phasesMatrix([j g],[k h]),betaAbs(g,h))*phasesMatrix(g,h);
    end
end

end

function S = computeSigns(phaseMatrix,beta)
% phase matrix has entries
% [alpha(j,k) alpha(j,h);
%  alpha(g,k) alpha(g,h)];

% beta =  abs(phaseMatrix(1,1)-phaseMatrix(1,2)-phaseMatrix(2,1)+phaseMatrix(2,2));
betaplus =  fixAngleRange(abs(phaseMatrix(1,1)-phaseMatrix(1,2)-phaseMatrix(2,1)+abs(phaseMatrix(2,2))));
betaminus =  fixAngleRange(abs(phaseMatrix(1,1)-phaseMatrix(1,2)-phaseMatrix(2,1)-abs(phaseMatrix(2,2))));

% disp([beta,betaplus,betaminus]*180/pi);

S = sign(abs(beta-betaminus) - abs(beta-betaplus));

end

function ang = fixAngleRange(ang)
% makes angle between 0 and pi such that cos(ang) remains the same

if ang>2*pi
    ang = ang-2*pi;
end

if ang>pi
    ang = 2*pi-ang;
end
end

