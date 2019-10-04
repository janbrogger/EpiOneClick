function varargout = BackgroundPowerByAge()

% varargin:
%           - None
% varargout:
%           - Background power by age. Each output struct contains all epochs for one age group.
%
% Dependencies 
%           - Output from function BackgroundPower() containing 2s-epochs
%           per EEG.

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = 'C:\Midlertidig_Lagring\Epileptiform\eegsnippetsplus';
if ~isdir(rootFolder)
    errorMessage = sprintf('Error, folder %s does not exist', rootFolder);
    uiwait(fprintf(errorMessage));
    return;
end

% Declare save folder for plots.
saveFolder = rootFolder;

% Subfolders/Users containing EEG-snippets
powerperEEG = load(fullfile(rootFolder, 'poweranalysis_study.mat'));

% create matrices for powerbands with gAgeCat as rows and epochs as collumns
gAgeCat = 1:1:10;
totalAgecatpower = cell(length(gAgeCat),1);
alphaAgecatpower = cell(length(gAgeCat),1);
betaAgecatpower = cell(length(gAgeCat),1);
thetaAgecatpower = cell(length(gAgeCat),1);
deltaAgecatpower = cell(length(gAgeCat),1);

% loop through studies and put epochs into matrix agecatpower(m x n) where m =
% gAgeCat and n = epochs
for i = 1 : length(powerperEEG.BackgroundPower.totalpower)
    agecat = powerperEEG.BackgroundPower.gAgeCat(i);
    totalAgecatpower{agecat} = [totalAgecatpower{agecat}, powerperEEG.BackgroundPower.totalpower{i}];
    alphaAgecatpower{agecat} = [alphaAgecatpower{agecat}, powerperEEG.BackgroundPower.powerAlpha{i}];
    betaAgecatpower{agecat} = [betaAgecatpower{agecat}, powerperEEG.BackgroundPower.powerBeta{i}];
    thetaAgecatpower{agecat} = [thetaAgecatpower{agecat}, powerperEEG.BackgroundPower.powerTheta{i}];
    deltaAgecatpower{agecat} = [deltaAgecatpower{agecat}, powerperEEG.BackgroundPower.powerDelta{i}];
end

% loop through agegroups, remove 20% of epochs (coarse artifact removal)
for i = 1 : length(gAgeCat)
    NEpochs = length(totalAgecatpower{i});
    trashepochs = ceil(0.2*NEpochs);
    totalAgecatpower{i} = sort(totalAgecatpower{i});
    totalAgecatpower{i} = totalAgecatpower{i}(1:(NEpochs-trashepochs));
    totalcatpower{i} = median(totalAgecatpower{i});
    alphaAgecatpower{i} = sort(alphaAgecatpower{i});
    alphaAgecatpower{i} = alphaAgecatpower{i}(1:(NEpochs-trashepochs));
    alphacatpower{i} = median(alphaAgecatpower{i});
    betaAgecatpower{i} = sort(betaAgecatpower{i});
    betaAgecatpower{i} = betaAgecatpower{i}(1:(NEpochs-trashepochs));
    betacatpower{i} = median(betaAgecatpower{i});
    thetaAgecatpower{i} = sort(thetaAgecatpower{i});
    thetaAgecatpower{i} = thetaAgecatpower{i}(1:(NEpochs-trashepochs));
    thetacatpower{i} = median(thetaAgecatpower{i});
    deltaAgecatpower{i} = sort(deltaAgecatpower{i});
    deltaAgecatpower{i} = deltaAgecatpower{i}(1:(NEpochs-trashepochs));
    deltacatpower{i} = median(deltaAgecatpower{i}); 
end

%plot median power by age for each power band
powerfig = figure();
set(powerfig, 'Name', 'Median power by age', 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]); %distorted by remote desktop?
subplot(3,2,1)
plot(cell2mat(deltacatpower), 'linewidth', 2);
ylim([0 inf]);
title('Median delta power by age')
subplot(3,2,2)
plot(cell2mat(thetacatpower), 'linewidth', 2);
ylim([0 inf]);
title('Median theta power by age')
subplot(3,2,3)
plot(cell2mat(alphacatpower), 'linewidth', 2);
ylim([0 inf]);
title('Median alpha power by age')
subplot(3,2,4)
plot(cell2mat(betacatpower), 'linewidth', 2);
ylim([0 inf]);
title('Median beta power by age')
subplot(3,2,5)
plot(cell2mat(totalcatpower), 'linewidth', 2);
ylim([0 inf]);
title('Median total power by age')

%save to file
saveplotpattern = fullfile(saveFolder, 'BgPowerByAge.png');
saveas(powerfig, saveplotpattern);
close(powerfig);

result.totalpower = totalAgecatpower;
result.alpha = alphaAgecatpower;
result.beta = betaAgecatpower;
result.theta = thetaAgecatpower;
result.delta = deltaAgecatpower;


varargout = {result};

