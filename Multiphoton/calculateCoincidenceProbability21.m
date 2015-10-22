function prob = calculateCoincidenceProbability21(amplitudeMatrix,phaseMatrix,freq,phiv1,phiv2,fv1,fv2,modeMatch,tauvals)
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
A = amplitudeMatrix;
P = phaseMatrix;


% fit the spectrum
[pp1,~,mup1] = polyfit(freq,phiv1,polyn);
[pp2,~,mup2] = polyfit(freq,phiv2,polyn);

% fit the filter functions
[pf1,~,muf1] = polyfit(freq,abs(fv1).^2,polyn);
[pf2,~,muf2] = polyfit(freq,abs(fv2).^2,polyn);


% compute landscape term
prob = zeros(1,length(tauvals));

parfor i=1:length(tauvals)
    
    landscapefun21 = @(x1,x2,x3) polyval(pf1,x1,[],muf1).*polyval(pf1,x2,[],muf1).*polyval(pf2,x3,[],muf2)...
        .*abs(A(1,1)^2*A(2,2)*polyval(pp1,x1,[],mup1).*polyval(pp1,x2,[],mup1).*polyval(pp2,x3,[],mup2).*exp(-1i*(x1+x2)*tauvals(i)+1i*(2*P(1,1)+P(2,2))) ...
        + A(1,1)*A(1,2)*A(2,1)*polyval(pp1,x1,[],mup1).*polyval(pp1,x3,[],mup1).*polyval(pp2,x2,[],mup2).*exp(-1i*(x1+x3)*tauvals(i)+1i*(P(1,1)+P(1,2)+P(2,1))) ...
        + A(1,1)*A(1,2)*A(2,1)*polyval(pp1,x2,[],mup1).*polyval(pp1,x3,[],mup1).*polyval(pp2,x1,[],mup2).*exp(-1i*(x2+x3)*tauvals(i)+1i*(P(1,1)+P(1,2)+P(2,1)))).^2;
    landscapefun12 = @(x1,x2,x3) polyval(pf1,x1,[],muf1).*polyval(pf2,x2,[],muf2).*polyval(pf2,x3,[],muf2)...
        .*abs(A(1,1)*A(1,2)*A(2,2)*polyval(pp1,x1,[],mup1).*polyval(pp1,x2,[],mup1).*polyval(pp2,x3,[],mup2).*exp(-1i*(x1+x2)*tauvals(i)+1i*(P(1,1)+P(1,2)+P(2,2))) ...
        + A(1,1)*A(1,2)*A(2,2)*polyval(pp1,x1,[],mup1).*polyval(pp1,x3,[],mup1).*polyval(pp2,x2,[],mup2).*exp(-1i*(x1+x3)*tauvals(i)+1i*(P(1,1)+P(1,2)+P(2,2))) ...
        + A(1,2)^2*A(2,1)*polyval(pp1,x2,[],mup1).*polyval(pp1,x3,[],mup1).*polyval(pp2,x1,[],mup2).*exp(-1i*(x2+x3)*tauvals(i)+1i*(2*P(1,2)+P(2,1)))).^2;
    
    prob(i) = integral3(landscapefun21,freq(1),freq(end),freq(1),freq(end),freq(1),freq(end)) ...
        + integral3(landscapefun12,freq(1),freq(end),freq(1),freq(end),freq(1),freq(end));
    
end

end
