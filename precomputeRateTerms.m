%% Precompute rate terms
% Compute integrals in the coincidence expression
%% Inputs:
%       tauvals: 1xn vector of delay values.
%       freq: px1 vector of frequencies at which all of filters,and spectra
%               have been measured.
%       phiv1: px1 vector of the first photon's spectrum.
%       phiv2: px1 vector of the second photon's spectrum.
%       fv1: px1 vector of the first filter function.
%       fv2: px1 vector of the second filter function.
%% Outputs:
%       landscapeCosFit: linear fit of the cosine dependent integral.
%       landscapeSinFit: linear fit of the sine dependent integral.
%       background1: the background term corresponding to two reflections.
%       background2: the background term corresponding to two transmissions.
%% Procedure
function [landscapeCosFit,landscapeSinFit,background1,background2] =...
    precomputeRateTerms(tauvals,freq,phiv1,phiv2,fv1,fv2)

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


% compute background term
background1 = integral(backfun11,freq(1),freq(end))*integral(backfun12,...
    freq(1),freq(end));
background2 = integral(backfun21,freq(1),freq(end))*integral(backfun22,...
    freq(1),freq(end));

% compute landscape term
landscapeCos = zeros(1,length(tauvals));
landscapeSin = zeros(1,length(tauvals));
for i=1:length(tauvals)

    % functions for intgrals
    landscapefun1c = @(x) polyval(pf1,x,[],muf1).*polyval(pp1,x,[],mup1)...
        .*polyval(pp2,x,[],mup2).*cos(x*tauvals(i));
    landscapefun2c = @(x) polyval(pf2,x,[],muf2).*polyval(pp1,x,[],mup1)...
        .*polyval(pp2,x,[],mup2).*cos(x*tauvals(i));
    landscapefun1s = @(x) polyval(pf1,x,[],muf1).*polyval(pp1,x,[],mup1)...
        .*polyval(pp2,x,[],mup2).*sin(x*tauvals(i));
    landscapefun2s = @(x) polyval(pf2,x,[],muf2).*polyval(pp1,x,[],mup1)...
        .*polyval(pp2,x,[],mup2).*sin(x*tauvals(i));
    
    % compute all four integrals
    landscape1c = integral(landscapefun1c,freq(1),freq(end));
    landscape2c = integral(landscapefun2c,freq(1),freq(end));
    landscape1s = integral(landscapefun1s,freq(1),freq(end));
    landscape2s = integral(landscapefun2s,freq(1),freq(end));
    
    % the landscape terms
    landscapeCos(i) = landscape2c.*landscape1c + landscape2s.*landscape1s;
    landscapeSin(i) = landscape2s.*landscape1c - landscape2c.*landscape1s;
end

% scale and fit the landscape
landscapeCosFit = fit(tauvals',landscapeCos','linearinterp');
landscapeSinFit = fit(tauvals',landscapeSin','linearinterp');

end

