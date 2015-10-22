function prob = calculateCoincidenceProbability02(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)
% inputs:
%       amplitudeMatrix: 2x2 matrix containting the amplitudes of the transition matrix.
%       phaseMatrix: 2x2 matrix containting the phases of the transition matrix.
%       phiv1: nx1 vector of the first photon's spectrum.
%       phiv2: nx1 vector of the second photon's spectrum.
%       fv1: nx1 vector of the first filter function.
%       fv2: nx1 vector of the second filter function.
%       modeMatch: scaler. 1 means complete spatial mode overlap, 0 means
%                   no spatial mode overlap.
%       tauvals: 1xm vector of delay values
% output:
%     prob: a vector, the same size as tauvals, of the coincidence
%           probability correspoding to the delays values in tauvals

polyn = 20; % degree to which to fit the spectrum and filters.

% fit the spectrum
[pp1,~,mup1] = polyfit(freq,phiv1,polyn);
[pp2,~,mup2] = polyfit(freq,phiv2,polyn);

% fit the filter functions
[pf1,~,muf1] = polyfit(freq,abs(fv1).^2,polyn);
[pf2,~,muf2] = polyfit(freq,abs(fv2).^2,polyn);

ampfun = @(x1,x2) polyval(pf1,x1,[],muf1).*polyval(pf2,x2,[],muf2).*abs(polyval(pp1,x1,[],mup1).*polyval(pp2,x2,[],mup2) + polyval(pp1,x2,[],mup1).*polyval(pp2,x1,[],mup2)).^2;

% compute background term
background = abs(amplitudeMatrix(2,1)*amplitudeMatrix(2,2))^2/2*integral2(ampfun,freq(1),freq(end),freq(1),freq(end));
% compute landscape term
prob = ones(1,length(tauvals))*background;

end
