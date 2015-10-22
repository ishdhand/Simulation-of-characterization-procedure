%% SIMULATION-NON-GAUSSIAN-FILTER-FILTER-NON-GAUSSIAN-FIT Simulation of characterization procedure
% This code generates a haar random unitary matrix, simulates one- and
% two-photon data and uses this simulated data to characterize the
% interferometer using our accurate and precise characteriation procedure.
% Please contact Ish Dhand <ishdhand@gmail.com> with questions, comments
% and bug-fixes.

%% Inputs

msizeVec = 5;
% We simulate the characterization experiment for interferometers with
% different number of channels m, choosing values of m from msizeVec.

clickVec = 10.^[8:0.5:8];
% The characterization experiment is performed till a given number of
% coincidences are measured at large time-delay values. We simulate the
% characterization for each of the different values of coincidences given
% in clickVec.

nSimulations = 10;
% We simulate the characterization experiment repeatedly (nSimulations
% times) to estimate inaccuracy in the characterization procedure.

modeMatch = 0.96;
% The mode matching parameter gamma. A calibration experiment is performed
% to estimate its value, which is assumed to be constant.

c = 3e8; % Speed of light. Please don't change.
eta = 0; % Multiphoton emission in the simulated data.
%% Outputs

Inaccuracy = zeros(length(msizeVec),length(clickVec),nSimulations);
% This variable stores the innacuracy in the characterized interferometer
% matrix. The inaccuracy is defined as the trace distance between the
% sampled haar-random unitary matrix and the matrix returned by the
% characterization procedure.

%% Experimental setup
% The following parameters describe the experimental setup, including the
% spectral filers, the calibration beamsplitter, and the translation stage.

filtervector1 = csvread('experimentaldata/filter3.csv');
filtervector2 = csvread('experimentaldata/filter4.csv');
freq = filtervector1(:,1);
fv1 = filtervector1(:,2);
fv2 = filtervector1(:,2);
phiv1 = ones(length(freq),1);
phiv2 = ones(length(freq),1);
% Read off the filter functions. Experimental data courtesy of He Lu at the
% University of Science and Technology of China.

centralwavelengths = 2*pi*c*1e9./[sum(freq.*fv1/sum(fv1)) ...
    sum(freq.*fv2/sum(fv2))];
spectrawidth = sqrt([var(freq,fv1.^2) var(freq,fv2.^2)]).*...
    centralwavelengths.^2*1e-9/(2*pi*c);
% Estimate the central wavelenght and the spectral width of the spectral
% functions. These two estimates are used in the data fitting.

reflectivityCalibBS = 0.5;
theta = asin(sqrt(reflectivityCalibBS));
% Reflectivity of beamslitter used in our calibration procedure.
amplitudeMatrixCalibBS = [cos(theta) sin(theta);sin(theta) cos(theta)];
phaseMatrixCalibBS = [0 pi;0 0];
% Interferometer matrix for the calibrating beamsplitter.

inputLosses = 0.31*ones(1,max(msizeVec));
outputLosses = 0.27*ones(1,max(msizeVec));
% Losses at the input and output ports

spectralWidths = 3;
% The largest value of time delay is equation to spectralWidths number.
taulimit = spectralWidths*max(centralwavelengths)^2.*1e-9...
    /(2*pi*c*min(spectrawidth));
tauN = 50;
% Number of values of time delay at which the experiment is performed.
tauvals = linspace(-taulimit,taulimit,tauN);
% The time delay values at which to simulate experimental data.

% Number of single photon pulses impinged at the interferometer.
numberSinglePulses = 1e8;

%
invExpTrans = zeros(1,3);
invExpTrans(2) = 3;
% Central position of the translation stage, i.e, the position that
% corresponds to zero time delay. Number in cm.
% invExpTrans(1) decides the number of incident photons. Will be set during
% simulation procedure.
invExpTrans(3) = 100*c;
% Speed of light (in cm/second, sorry!)

% Turn the warnings off here
% warning('off','stats:nlinfit:IllConditionedJacobian');
% warning('off','stats:nlinfit:ModelConstantWRTParam');
% warning('off','MATLAB:rankDeficientMatrix');

%% Sample unitary, simulate data and characterize

% filters and spectra
% [freq, phiv1, phiv2, fv1, fv2] = generateGaussianSpectrum...
%    (centralwavelengths,spectrawidth,p);

Counter2 = 0;
for msize = msizeVec
    % Choose number of interferometer channels from msizeVec
    m = msize;
    Counter2 = Counter2+1;
    Counter3 = 0;
    for clicks = clickVec
        % clicks decides the number of coincidences for large time delay
        % values.Basically, more clicks => less relative shot noise.
        Counter3 = Counter3+1;
        % fix click params
        invExpTrans(1) = clicks;
        
        % Perform simulation nSimulations times.
        for Counter1=1:nSimulations
            
            fprintf(['#channels = ',num2str(msize),',\t Log #photons = '...
                ,num2str(log10(clicks)),',\t #Run= ' ,num2str(Counter1)...
                ,'\n']);
            % Print diagnostic information
            
            interferometerMatrix = generateHaarUnitary(m);
            % sample $n \times n$  unitary matrix U from Haar measure.
            interferometerMatrixClean = removeInputOutputPhases...
                (interferometerMatrix);
            % remove input output phases of interferometer matrix
            
            %% Simulate one- and two-photon data
            
            % Two-photon data from calibration beamsplitter
            [simulatedtauvalsCalibBS,simulatedrateCalibBS] = ...
                generateSimulatedHOMData(tauvals,invExpTrans,...
                amplitudeMatrixCalibBS,phaseMatrixCalibBS,modeMatch,...
                freq,phiv1,phiv2,fv1,fv2,eta);
            
            % One-photon data
            singleClicksData = generateSimulatedSingleClicks...
                (numberSinglePulses,eta,interferometerMatrix,...
                inputLosses,outputLosses);
            
            % hom data
            [simulatedtauvals,simulatedRatesInterferometer] = ...
                generateSimulatedHOMDataInterferometer(tauvals,...
                invExpTrans,interferometerMatrix,modeMatch,freq,phiv1,...
                phiv2,fv1,fv2,eta);
            
            %% Reconstruction: Initial conditions and precomputation
            
            % The following vectors are used to set the initial conditions
            % in the curve fitting algorithms.
            scaleVecCalibration = [1/70 clicks/50 1/10 2*c];
            scaleVecPhase = [1/15 clicks/50 1/10 2*c];
            
            % The following code precomputes some integrals that are used
            % in the curve-fitting subroutines.
            [landscapeCosFit,landscapeSinFit,background1,background2] = ...
                precomputeRateTerms(tauvals,freq,phiv1,phiv2,fv1,fv2);
            
            %% Calibration
            
            % integration scaling parameter
            sp = scalingParameter(amplitudeMatrixCalibBS,freq,phiv1,...
                phiv2,fv1,fv2);
            
            % guess the optimization parameters
            beta0 = HOMFittingGuess(simulatedtauvalsCalibBS,...
                simulatedrateCalibBS,amplitudeMatrixCalibBS,'mode match'...
                ,reflectivityCalibBS,freq,phiv1,phiv2,fv1,fv2);
            
            %% Calibrate the mode match
            [modeMatchEstimated,expTransCalib] = Calibration...
                (simulatedrateCalibBS,simulatedtauvalsCalibBS,beta0,...
                reflectivityCalibBS,landscapeCosFit,landscapeSinFit,...
                background1,background2,sp,scaleVecCalibration);
            tauScaling = expTransCalib(3)/scaleVecPhase(4);
            
            %           modeMatchPrecision = BootstrapModeMatch(...
            %                 simulatedrateCalibBS,simulatedtauvalsCalibBS,...
            %                 expTransCalib,modeMatchEstimated,reflectivityCalibBS...
            %                 ,landscapeCosFit,landscapeSinFit,background1,...
            %                 background2,sp,scaleVecCalibration,bootstrapSamples);
            
            %           fprintf('Mode match parameter is %3.2f with precision %g.\n'...
            %                 ,modeMatchEstimated,modeMatchPrecision);
            
            % Amplitude estimation
            [relativeAmplitudesMatrix,relativeAmplitudeErrors] = ...
                AmplitudeEstimation(singleClicksData);
            
            %% Argument estimation
            [phasesMatrix,expTransArray] = PhaseEstimation(...
                simulatedRatesInterferometer,simulatedtauvals,...
                relativeAmplitudesMatrix,modeMatchEstimated,freq,phiv1,...
                phiv2,fv1,fv2,landscapeCosFit,landscapeSinFit,...
                background1,background2,scaleVecPhase,tauScaling);
            
            % precision of phases
            % phasesMatrixPrecision = BootstrapPhases(...
            %     simulatedRatesInterferometer,simulatedtauvals,...
            %     expTransArray,modeMatchEstimated,...
            %     relativeAmplitudesMatrix,phasesMatrix,freq,phiv1,...
            %     phiv2,fv1,fv2,landscapeCosFit,landscapeSinFit,...
            %     background1,background2,sp,scaleVecPhase,...
            %     tauScaling,bootstrapSamples);
            
            %% Unitarity
            TransitionMatrix = UnitaryConstraint(...
                relativeAmplitudesMatrix,phasesMatrix);
               
            delta = TransitionMatrix-interferometerMatrixClean;
            Inaccuracy(Counter2,Counter3,Counter1) = ...
                trace(sqrtm(delta'*delta))/2;

            %% Visualize and write
            % The following code enables visualizing the difference between
            % the matrices
            clf
            % visualizeDelta(sign(angle(interferometerMatrixClean(2:m,...
            %    2:m)))-sign(phasesMatrix(2:m,2:m)));
            visualizeDelta(abs(angle(interferometerMatrixClean(2:m,...
                2:m))-phasesMatrix(2:m,2:m)));
            drawnow
            saveas(gcf,[pwd,'/Datasets/AccurateAndPrecise/Delta_',num2str(msize),'Channel'...
                ,num2str(log10(clicks)),'e10Photons',num2str(Counter1)...
                ,'.eps'],'psc2');
            
            % The following code writes the matrices to csv files in the
            % folder DataSets
            csvwrite([pwd,'/Datasets/AccurateAndPrecise/U_',num2str(msize),'Channel'...
                ,num2str(log10(clicks)),'e10Photons',num2str(Counter1)...
                ,'.csv'],interferometerMatrixClean);
            csvwrite([pwd,'/Datasets/AccurateAndPrecise/V_',num2str(msize),'Channel'...
                ,num2str(log10(clicks)),'e10Photons',num2str(Counter1)...
                ,'.csv'],TransitionMatrix);
            csvwrite([pwd,'/Datasets/AccurateAndPrecise/Delta_',num2str(msize),'Channel'...
                ,num2str(log10(clicks)),'e10Photons',num2str(Counter1)...
                ,'.csv'],delta);
        end 
    end
end
%save('resourcesDataExperimentalFiltersHighRes');
csvwrite([pwd,'/Datasets/AccurateAndPrecise/Inaccuracy.csv'],Inaccuracy);