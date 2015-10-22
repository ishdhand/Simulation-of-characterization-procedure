%% Calculate Coincidence Probability 
% Calculate coincidence probability when single photons with known spectra
% are incident at a known interferometer.

%% Inputs:
%       amplitudeMatrix: 2x2 matrix containing the amplitudes of the transition matrix.
%       phaseMatrix: 2x2 matrix containting the phases of the transition matrix.
%       freq: nx1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: nx1 vector of the first photon's spectrum.
%       phiv2: nx1 vector of the second photon's spectrum.
%       fv1: nx1 vector of the first filter function.
%       fv2: nx1 vector of the second filter function.
%       modeMatch: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       tauvals: 1xm vector of delay values
%% Output:
%     prob: a vector, the same size as tauvals, of the coincidence
%           probability correspoding to the delays values in tauvals
%% Procedure
function prob = calculateCoincidenceProbability112(amplitudeMatrix,...
    phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)

polyn = 20; % degree to which to fit the spectrum and filters.

% fit the spectrum
[pp1,~,mup1] = polyfit(freq,phiv1,polyn);
[pp2,~,mup2] = polyfit(freq,phiv2,polyn);

% squared spectrum
[pp1s,~,mup1s] = polyfit(freq,phiv1.^2,polyn);
[pp2s,~,mup2s] = polyfit(freq,phiv2.^2,polyn);

% fit the filter functions
[pf1,~,muf1] = polyfit(freq,abs(fv1).^2,polyn);
[pf2,~,muf2] = polyfit(freq,abs(fv2).^2,polyn);

backfun11 = @(x) polyval(pp1s,x,[],mup1s).*polyval(pf1,x,[],muf1);
backfun12 = @(x) polyval(pp2s,x,[],mup2s).*polyval(pf2,x,[],muf2);
backfun21 = @(x) polyval(pp1s,x,[],mup1s).*polyval(pf2,x,[],muf2);
backfun22 = @(x) polyval(pp2s,x,[],mup2s).*polyval(pf1,x,[],muf1);

% some constants
crossCoefficient = 2*modeMatch*amplitudeMatrix(1,1)*amplitudeMatrix(1,2)*...
    amplitudeMatrix(2,1)*amplitudeMatrix(2,2);
phaseSum = phaseMatrix(1,1)+phaseMatrix(2,2)-phaseMatrix(1,2)-phaseMatrix(2,1);

% compute background term
background = amplitudeMatrix(1,1)^2*amplitudeMatrix(2,2)^2*integral...
    (backfun11,freq(1),freq(end))*integral(backfun12,freq(1),freq(end))+...
    amplitudeMatrix(1,2)^2*amplitudeMatrix(2,1)^2*integral(backfun21,...
    freq(1),freq(end))*integral(backfun22,freq(1),freq(end));

% compute landscape term
landscape1c = zeros(1,length(tauvals));
landscape2c = zeros(1,length(tauvals));
landscape1s = zeros(1,length(tauvals));
landscape2s = zeros(1,length(tauvals));

for i=1:length(tauvals)
    
    landscapefun1c = @(x) polyval(pf1,x,[],muf1).*polyval(pp1,x,[],mup1).*...
        polyval(pp2,x,[],mup2).*cos(x*tauvals(i));
    landscapefun2c = @(x) polyval(pf2,x,[],muf2).*polyval(pp1,x,[],mup1).*...
        polyval(pp2,x,[],mup2).*cos(x*tauvals(i));
    landscapefun1s = @(x) polyval(pf1,x,[],muf1).*polyval(pp1,x,[],mup1).*...
        polyval(pp2,x,[],mup2).*sin(x*tauvals(i));
    landscapefun2s = @(x) polyval(pf2,x,[],muf2).*polyval(pp1,x,[],mup1).*...
        polyval(pp2,x,[],mup2).*sin(x*tauvals(i));
    
    landscape1c(i) = integral(landscapefun1c,freq(1),freq(end));
    landscape2c(i) = integral(landscapefun2c,freq(1),freq(end));
    landscape1s(i) = integral(landscapefun1s,freq(1),freq(end));
    landscape2s(i) = integral(landscapefun2s,freq(1),freq(end));
    
    
    
end

% the final coincidence prob
% prob = background+crossCoefficient*landscape;
prob = background + crossCoefficient*cos(phaseSum)*(landscape2c.*landscape1c...
    + landscape2s.*landscape1s)- crossCoefficient*sin(phaseSum)*...
    (landscape2s.*landscape1c - landscape2c.*landscape1s);


end
