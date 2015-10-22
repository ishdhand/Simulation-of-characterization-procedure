%% HOM Fitting guess
% Guesses the values of the parameters that fit the
% theoretical curve to the experimental rate.
%% Inputs:
%       exptauvals: 1xn vector of delay values from experiment.
%       exprate: 1xn vector of the experimental
%           coincidence rates correspoding to the delays values in tauvals.
%       amplitudeMatrix: 2x2 matrix that has the amplitudes of the
%           interferometer. The amplitudes may only be correct up to losses
%           at the inputs and outputs of the interferometer.
% if the fourth input is 'mode match' then the subsequent inputs are
%       reflectivityCalibBS: reflectivity of calibrating beamsplitter
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
% if the fourth input is 'phase' then the subsequent inputs are
%       phaseMatrix: 2x2 matrix containting the phases of the transition
%           matrix. The last entry is to be estimated and it's value is not
%           used.
%       modeMatch: scaler. 1 means complete mode overlap, 0 means
%                   no mode overlap.
%% Outputs:
%       expTransGuess: for fourth argument 'mode match' 1x4 vector of
%           guesses. The guesses are modeMatch, rate scale, tau shift, tau
%           scale. For fourth argument 'phase' 1x3 vector of guesses. The
%           guesses are phase, rate scale, tau shift.
%% Procedure

function expTransGuess = HOMFittingGuess(exptauvals,exprate,amplitudeMatrix,varargin)

expTransGuess = zeros(1,3);
c = 3e8;

L = length(exprate);
l = round(0.1*L);

% various exp rates
minRate = min(exprate);
maxRate = max(exprate);
wingRate = mean([exprate(1:l) exprate(end-(l-1):end)]);

% various theory rates
baseRate = amplitudeMatrix(1,1)^2*amplitudeMatrix(2,2)^2 + amplitudeMatrix(1,2)^2*amplitudeMatrix(2,1)^2;

% shifting of tau values so stationary point is at zero
if abs(wingRate-minRate) > abs(wingRate -maxRate) % normal HOM curve
    expTransGuess(3) = mean(exptauvals(exprate == minRate));
    zeroRate = minRate;
else %inverted HOM curve
    expTransGuess(3) = mean(exptauvals(exprate == maxRate));
    zeroRate = maxRate;
end

% scaling of rate
% expTransGuess(2) = wingRate/baseRate;
expTransGuess(2) = wingRate;


% the experimental visibility is needed for mode match or phase
expvisibility = (wingRate - zeroRate)/wingRate;

% estimate mode match, needs reflectivity
if strcmp(varargin{1},'mode match')
    
    % scaling of delay.
    freq = varargin{3};
    phiv1 = varargin{4};
    phiv2 = varargin{5};
    fv1 = varargin{6};
    fv2 = varargin{7};
    
    % width of spectra
    sf1 = var(freq,fv1)/2;
    sf2 = var(freq,fv2)/2;
    sp1 = var(freq,phiv1)/2;
    sp2 = var(freq,phiv2)/2;
    
    % width of HOM curve can be ascertained via these
    sig1 = (sf1*sp1*sp2)/(2*sp1*sp2 + sf1*(sp1 + sp2));
    sig2 = (sf2*sp1*sp2)/(2*sp1*sp2 + sf2*(sp1 + sp2));
    sigT = 1/(sig1+sig2);
    sigE = var(exptauvals,exprate)/2;
    expTransGuess(4) = sqrt(sigE/sigT);
    
    % beamsplitter angle from reflectivity
    theta = asin(sqrt(varargin{2}));
    
    expectedvisibility = 2*cos(theta)^2*sin(theta)^2/(cos(theta)^4 + sin(theta)^4);
    % the mode match is the ratio of the experimental visibility to the
    % expected visibility
    expTransGuess(1) = expvisibility/expectedvisibility;
    
elseif strcmp(varargin{1},'phase')
    % estimate phase, needs 3 other phases
    phaseMatrix = varargin{2};
    modeMatch = varargin{3};
    b = 2*prod(prod(amplitudeMatrix))*modeMatch;
    d = -baseRate/b*expvisibility;
    % ensure d remains bounded
    if d<-1
        d = -1;
    elseif d>1
        d = 1;
    end
    expTransGuess(1) = -(acos(d) - (phaseMatrix(1,1)-phaseMatrix(1,2)-phaseMatrix(2,1)));
end

end

