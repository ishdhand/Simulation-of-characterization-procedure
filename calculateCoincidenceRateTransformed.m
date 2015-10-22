%% Calculate coincidence rate transformed
% Computes the expected coincidence rate using integrals over spectra of 
% source field.

%% Inputs:
%       exptauvals: 1xn values of relative delay for which the rate is computed
%       expTrans: 1x3 vector of transformations: rate scale, tau
%           shift, tau scale.
%       amplitudeMatrix: 2x2 matrix containing the amplitudes of the transition matrix.
%       phaseMatrix: 2x2 matrix containting the phases of the transition matrix.
%       modeMatch: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%       sp: scaler integral scaling paramater.

%% outputs
%       exprate: 1xn values of expected coincidence rates in an experiment

%% Procedure

function exprate = calculateCoincidenceRateTransformed(exptauvals,...
    expTrans,amplitudeMatrix,phaseMatrix,modeMatch,landscapeCosFit,...
    landscapeSinFit,background1,background2,sp)

% shift and scale tau
tauvals = (exptauvals - expTrans(2))/expTrans(3);

% background coefficients
a1 = amplitudeMatrix(1,1)^2*amplitudeMatrix(2,2)^2;
a2 = amplitudeMatrix(1,2)^2*amplitudeMatrix(2,1)^2;

% cross coefficient
cc = 2*modeMatch*amplitudeMatrix(1,1)*amplitudeMatrix(1,2)*...
    amplitudeMatrix(2,1)*amplitudeMatrix(2,2);

% phase sum term
phaseSum = phaseMatrix(1,1)+phaseMatrix(2,2)-phaseMatrix(1,2)-...
    phaseMatrix(2,1);

% the prob
prob = a1*background1+a2*background2 + cc*(cos(phaseSum)*...
    landscapeCosFit(tauvals) + sin(phaseSum)*landscapeSinFit(tauvals));

% scale prob by intergration scale
probs = sp^4*prob;

% scale rate by experimental scale
% exprate = expTrans(1)*probs';
S = expTrans(1)/probs(1);
exprate = S*probs';

end

