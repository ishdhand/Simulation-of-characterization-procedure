% Generate simulated single clicks
% Generates simulated single click data.
% Inputs:
%       NumberPulses: the number of photons sent into the setup.
%       eta: positive real number between 0 and 1 inclusive. Probability of
%           multiphoton events.
%       interferometerMatrix: mxm unitary matrix.
%       inputLosses: 1xm vector that contains the losses at the input of
%           interferometer.
%       outputLosses: 1xm vector that contains the losses at the output of
%           interferometer.
% Outputs:
%       SingleClicksData: mxm matrix of single click data with shot noise.

function SingleClicksData = generateSimulatedSingleClicks(NumberPulses,eta,interferometerMatrix,inputLosses,outputLosses)

% size of the interferometer
m = length(interferometerMatrix);


SingleClicksData = zeros(m,m);
for j=1:m
    for k=1:m
        rj = outputLosses(j);
        sk = inputLosses(k);
        tjk = abs(interferometerMatrix(j,k));
        % expression for single click rate correct upto two photons.
        % input channel k, detection channel j
        meanClicks = NumberPulses*(rj*tjk^2*sk + eta^2*rj*tjk^2*(1-sk*tjk^2)*sk + eta^2*rj*(2-rj)*tjk^4*sk^2);
        % add shot noise
        SingleClicksData(j,k) = poissrnd(meanClicks);
        
    end
end


end

