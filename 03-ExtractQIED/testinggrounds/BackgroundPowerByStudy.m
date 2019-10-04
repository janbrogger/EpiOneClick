function varargout = BackgroundPowerByStudy()

% varargin:
%           - None
% varargout:
%           - Table with one row per study, mean bandpowers for each study
%           after removal of 20% max power epochs per agecat. 
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

% Declare save folder for plots and files.
saveFolder = rootFolder;

% Subfolders/Users containing EEG-snippets
powerperEEG = load(fullfile(rootFolder, 'poweranalysis_study.mat'));

% create matrices for powerbands with gAgeCat as rows and epochs as collumns
gAgeCat = 1:1:10;
NStudies = length(powerperEEG.BackgroundPower.StudyId);
totalAgecatpower = cell(length(gAgeCat),1);
thresholdpower = zeros(length(gAgeCat),1);

% loop through studies and put epochs into matrix agecatpower(m x n) where m =
% gAgeCat and n = epochs across studies
for i = 1 : length(powerperEEG.BackgroundPower.StudyId)
    agecat = powerperEEG.BackgroundPower.gAgeCat(i);
    totalAgecatpower{agecat} = [totalAgecatpower{agecat}, powerperEEG.BackgroundPower.totalpower{i}];
end

% loop through agegroups, identify 20% power threshold for epochs in each agecat (coarse artifact removal)
for i = 1 : length(gAgeCat)
    NEpochs = length(totalAgecatpower{i});
    epochthresholdidx = ceil(0.8*NEpochs);
    sortedagecat = sort(totalAgecatpower{i});
    thresholdpower(i) = sortedagecat(epochthresholdidx);
end

powerbystudy = zeros(NStudies, 7);
% loop through studies, remove epochs above threshold, compute average, produce table for STATA
for i = 1 : NStudies
    agecat = powerperEEG.BackgroundPower.gAgeCat(i);
    studyid = powerperEEG.BackgroundPower.StudyId(i);
    SRID = powerperEEG.BackgroundPower.SRID(i);
    totalpower = powerperEEG.BackgroundPower.totalpower{i};
    powerDelta = powerperEEG.BackgroundPower.powerDelta{i};
    powerTheta = powerperEEG.BackgroundPower.powerTheta{i};
    powerAlpha = powerperEEG.BackgroundPower.powerAlpha{i};
    powerBeta = powerperEEG.BackgroundPower.powerBeta{i};
    toberemoved = find(totalpower >= thresholdpower(agecat));
    totalpower(toberemoved) = [];
    powerDelta(toberemoved) = [];
    powerTheta(toberemoved) = [];
    powerAlpha(toberemoved) = [];
    powerBeta(toberemoved) = [];
    totalpowermedian = median(totalpower);
    deltamedian = median(powerDelta);
    thetamedian = median(powerTheta);
    alphamedian = median(powerAlpha);
    betamedian = median(powerBeta);   
    
    powerbystudy(i,:) = [studyid SRID totalpowermedian deltamedian thetamedian alphamedian betamedian];
end

T = array2table(powerbystudy, 'VariableNames', {'StudyId', 'SRID', 'medianTotalPower', 'medianDelta', 'medianTheta', 'medianAlpha', 'medianBeta'});
filename_excel = fullfile(saveFolder, 'PowerByStudy.xlsx');
if exist(filename_excel, 'file')
  delete(filename_excel);
end
writetable(T, filename_excel, 'FileType', 'spreadsheet');

varargout = {T};

