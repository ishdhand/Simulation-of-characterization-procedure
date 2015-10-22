%$ SCALING-PARAMETER 
% Computes a parameter with which to scale the filters and
% spectrum so as to make integrals smooth and keep all coincidence plots
% have the same vertical scaling.

% inputs:
%       amplitudeMatrix: 2x2 matrix containing the amplitudes of the transition matrix.
%       freq: nx1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: nx1 vector of the first photon's spectrum.
%       phiv2: nx1 vector of the second photon's spectrum.
%       fv1: nx1 vector of the first filter function.
%       fv2: nx1 vector of the second filter function.
% outputs:
%       sp: scaler integral scaling paramater. Multiply all spectra and
%           filters with this parameter. Divide coincidince prob by 
%           sp^(2*number of photons in system) to obtain scaled plot.


function sp = scalingParameter(amplitudeMatrix,freq,phiv1,phiv2,fv1,fv2)

polyn = 20;

% background which is expected
backgroundExpected = amplitudeMatrix(1,1)^2*amplitudeMatrix(2,2)^2 ...
    + amplitudeMatrix(1,2)^2*amplitudeMatrix(2,1)^2;

% fit the spectra and filter functions
[pp1s,~,mup1s] = polyfit(freq,phiv1.^2,polyn);
[pp2s,~,mup2s] = polyfit(freq,phiv2.^2,polyn);

[pf1,~,muf1] = polyfit(freq,abs(fv1).^2,polyn);
[pf2,~,muf2] = polyfit(freq,abs(fv2).^2,polyn);

% functions that are input to the integral function.
backfun11 = @(x) polyval(pp1s,x,[],mup1s).*polyval(pf1,x,[],muf1);
backfun12 = @(x) polyval(pp2s,x,[],mup2s).*polyval(pf2,x,[],muf2);
backfun21 = @(x) polyval(pp1s,x,[],mup1s).*polyval(pf2,x,[],muf2);
backfun22 = @(x) polyval(pp2s,x,[],mup2s).*polyval(pf1,x,[],muf1);

% compute the background with the current scaling.
background = amplitudeMatrix(1,1)^2*amplitudeMatrix(2,2)^2*...
    integral(backfun11,freq(1),freq(end))*integral(backfun12,...
    freq(1),freq(end))+ amplitudeMatrix(1,2)^2*...
    amplitudeMatrix(2,1)^2*integral(backfun21,freq(1),freq(end))...
    *integral(backfun22,freq(1),freq(end));

% the scaling parameter
sp = (backgroundExpected/(background))^0.25;


end



