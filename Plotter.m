clc; clf
clear

NoCal = csvread('DataSets/InaccuracyWOCalibration.csv');
Cal = csvread('DataSets/InaccuracyWCalibration.csv');

modeMatchVec = fliplr(100*(1-[1-2.^(-4:-0.5:-10),1]));

NoCalMean = fliplr(mean(abs(NoCal')));
NoCalStd = fliplr(std(abs(NoCal')));
CalMean = fliplr(mean(abs(Cal')));
CalStd = fliplr(std(abs(Cal')));
NoCalMean = NoCalMean(2:14);
NoCalStd= NoCalStd(2:14);
CalMean = CalMean(2:14);
CalStd= CalStd(2:14);

Figure1 = figure(1); clf;
errorbar(log10(CalMean),CalStd,'LineWidth',2);
hold
errorbar(log10(NoCalMean),NoCalStd,'LineWidth',2);
title('Effect of calibration on accuracy of characterization');
set(0,'defaulttextinterpreter','latex');
xlabel('mode mismatch = 1-\gamma')
ylabel('log_{10}(mean error)')
legend('With calibration','Without calibration','Location','northwest')
ax = gca;
ax.FontName = 'Times';
ax.FontSize=16;
ax.XTick = [1:2:13];
ax.XTickLabel = [0.001,0.002,0.004,0.008,0.016,0.032,0.064];
ax.XTickLabelRotation = 45;
ax.YTickLabelRotation = 45;
saveas(gcf,strcat('IsCalibrationUseful','.eps'),'psc2');

% - % - % - %

NoGauss = csvread('DataSets/InaccuracyNonGauss.csv');
Gauss = csvread('DataSets/InaccuracyGaussian.csv');
NoGaussMean = fliplr(mean(abs(NoGauss')));
NoGaussStd = fliplr(std(abs(NoGauss')));
GaussMean = fliplr(mean(abs(Gauss')));
GaussStd = fliplr(std(abs(Gauss')));

clickVec = fliplr(10.^[4:0.5:10]);

Figure2 = figure(2); clf;
errorbar(log10(NoGaussMean),NoGaussStd,'LineWidth',2);
hold
errorbar(log10(GaussMean),GaussStd,'LineWidth',2);
title('Effect of fitting-curve choice on accuracy of characterization');
set(0,'defaulttextinterpreter','latex');
xlabel('log_{10} (shot noise)')
ylabel('log_{10} (mean error)')
legend('Fitting to correct curves','Fitting to Gaussian curves','Location','northwest')
ax = gca;
ax.FontName = 'Times';
ax.FontSize=16;
ax.XTick = 1:2:13;
XaxisLabel = -log10(sqrt(clickVec));
ax.XTickLabel = XaxisLabel(1:2:13);
ax.XTickLabelRotation = 45;
ax.YTickLabelRotation = 45;
saveas(gcf,strcat('GaussianVersusNongauss','.eps'),'psc2');

