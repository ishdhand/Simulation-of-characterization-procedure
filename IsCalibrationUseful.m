%% Is calibration useful?
% Simulation of characterization procedure for two cases: Using our
% accurate and precise charcterization procedure with and without
% calibration. In both cases, the coincidence data is obtained by 
% integrating over spectral functions assuming imperfect mode matching.
% In the first case, the mode-matching parameter is esimtated using
% simulated coincidence counts on a beamsplitter and this estimate is used
% to obtain accurate estimates of the interferometer parameters. In the
% second case, the simulated data has mode mismatch, but the
% characterization is performed assuming perfect mode matching. The two 
% cases are run repeatedly for different values of mode-matching
% parameter.
%
% Please contact Ish Dhand <ishdhand@gmail.com> with questions, comments
% and bug-fixes.

tic
%% Inputs
msizeVec = 4;
% We simulate the characterization experiment for interferometers with
% different number of channels m, choosing values of m from msizeVec.

clickVec = 10.^[8:0.5:8];
% The characterization experiment is performed till a given number of
% coincidences are measured at large time-delay values. We simulate the
% characterization for each of the different values of coincidences given
% in clickVec.

nSimulations = 20;
% We simulate the characterization experiment repeatedly (nSimulations
% times) to estimate inaccuracy in the characterization procedure.

modeMatchVec = [1-2.^(-4:-1:-10),1];
% The mode matching parameter gamma. A calibration experiment is performed
% to estimate its value, which is assumed to be constant.

c = 3e8; % Speed of light. Please don't change.
eta = 0; % Multiphoton emission in the simulated data.
%% Outputs

Inaccuracy = zeros(length(modeMatchVec),nSimulations);
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
 warning('off','stats:nlinfit:IllConditionedJacobian');
 warning('off','stats:nlinfit:ModelConstantWRTParam');
 warning('off','MATLAB:rankDeficientMatrix');

%% Sample unitary, simulate data and characterize

% filters and spectra
% [freq, phiv1, phiv2, fv1, fv2] = generateGaussianSpectrum...
%    (centralwavelengths,spectrawidth,p);

Counter2 = 0;
for modeMatch = modeMatchVec
    % Choose number of interferometer channels from msizeVec
    msize = msizeVec;
    Counter2 = Counter2+1;
    Counter3 = 0;
    for clicks = clickVec
        % clicks decides the number of coincidences for large time delay
        % values.Basically, more clicks => less relative shot noise.
        Counter3 = Counter3+1;
        % fix click params
        invExpTrans(1) = clicks;
        
        % Perform simulation nSimulations times.
        parfor Counter1=1:nSimulations
            
            fprintf(['#channels = ',num2str(msize),',\t Log #photons = '...
                ,num2str(log10(clicks)),',\t #Run= ' ,num2str(Counter1)...
                ,'\t ModeMatch = ',num2str(modeMatch),'\n']);
            % Print diagnostic information
            
            interferometerMatrix = generateHaarUnitary(msize);
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
            Inaccuracy(Counter2,Counter1) = ...
                trace(sqrtm(delta'*delta))/2;

            %% Visualize and write
            % The following code enables visualizing the difference between
            % the matrices
            clf
            % visualizeDelta(sign(angle(interferometerMatrixClean(2:m,...
            %    2:m)))-sign(phasesMatrix(2:m,2:m)));
            % visualizeDelta(abs(angle(interferometerMatrixClean(2:msize,...
            %     2:msize))-phasesMatrix(2:msize,2:msize)));
            
            % The following code writes the matrices to csv files in the
            % folder DataSets
            csvwrite(['CalU_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],interferometerMatrixClean);
            csvwrite(['CalV_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],TransitionMatrix);
            csvwrite(['CalDelta_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],delta);

            % plots;drawnow

        end 
    end
end
%save('resourcesDataExperimentalFiltersHighRes');
csvwrite('InaccuracyWCalibration.csv',Inaccuracy);

InaccuracyWCalibration = Inaccuracy;

Counter2 = 0;
for modeMatch = modeMatchVec
    % Choose number of interferometer channels from msizeVec
    msize = msizeVec;
    Counter2 = Counter2+1;
    Counter3 = 0;
    for clicks = clickVec
        % clicks decides the number of coincidences for large time delay
        % values.Basically, more clicks => less relative shot noise.
        Counter3 = Counter3+1;
        % fix click params
        invExpTrans(1) = clicks;
        
        % Perform simulation nSimulations times.
        parfor Counter1=1:nSimulations
            
            fprintf(['#channels = ',num2str(msize),',\t Log #photons = '...
                ,num2str(log10(clicks)),',\t #Run= ' ,num2str(Counter1)...
                ,'\t ModeMatch = ',num2str(modeMatch),'\n']);
            % Print diagnostic information
            
            interferometerMatrix = generateHaarUnitary(msize);
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
            modeMatchEstimated = 1;
            
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
            Inaccuracy(Counter2,Counter1) = ...
                trace(sqrtm(delta'*delta))/2;

            %% Visualize and write
            % The following code enables visualizing the difference between
            % the matrices
            % visualizeDelta(sign(angle(interferometerMatrixClean(2:m,...
            %    2:m)))-sign(phasesMatrix(2:m,2:m)));
            % visualizeDelta(abs(angle(interferometerMatrixClean(2:msize,...
            %      2:msize))-phasesMatrix(2:msize,2:msize)));
            % drawnow
%             figure
%             plots
%             drawnow

            
            % The following code writes the matrices to csv files in the
            % folder DataSets
            csvwrite(['NoCalU_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],interferometerMatrixClean);
            csvwrite(['NoCalV_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],TransitionMatrix);
            csvwrite(['NoCalDelta_',num2str(msize),'Channel',...
                num2str(log10(clicks)),'e10Photons',num2str(Counter1),...
                'ModeMatch',num2str(modeMatch),'.csv'],delta);
        end 
    end
end
%save('resourcesDataExperimentalFiltersHighRes');
csvwrite('InaccuracyWOCalibration.csv',Inaccuracy);
InaccuracyWOCalibration = Inaccuracy;

toc
 
% plot(InaccuracyWCalibration(:))
% hold 
% plot(InaccuracyWOCalibration(:))
