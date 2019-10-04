function varargout = BackgroundPower()

% varargin:
%           - None
% varargout:
%           - Background power calculated from total length of annotated IED electrode.
%
% Dependencies 
%           - ScorePipeline episnippets at location
%             C:\Midlertidig_Lagring\Epileptiform\eegsnippetsplus


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
folderPattern = fullfile(rootFolder, 'User*');
allFolders = dir(folderPattern);

%Struct for output
BackgroundPower.powerDelta = [];
BackgroundPower.powerTheta = [];
BackgroundPower.powerAlpha = [];
BackgroundPower.powerBeta = [];
BackgroundPower.totalpower = [];
BackgroundPower.SRID = [];
BackgroundPower.gAgeCat = [];
BackgroundPower.StudyId = [];
BackgroundPower.RMStest = [];

% Loop through subfolders (users)
for l = 1 : length(allFolders)
filePattern = fullfile(rootFolder, allFolders(l).name, 'Epi_SearchResultEventId*.mat');
theFiles = dir(filePattern);

%Loop through EEG-snippets
for k = 1 : length(theFiles)    
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
    fprintf(1, 'Reading file %d %s \n', k, fullFileName);
    
    % Load episnippet signal from one channel.
    % Measure bandpower and total power in consecutive epochs of length 2 seconds (1000
    % samples)
    epiAnno = load(fullFileName);
    srate = epiAnno.epiSnippetPlus.srate;
    % IEDlength = epiAnno.epiSnippetPlus.afterDischargeEnd - epiAnno.epiSnippetPlus.spikeStart;
    SearchResultID = epiAnno.epiSnippetPlus.SearchResultEventId;
    gAgeCat = gAgeCat2num(epiAnno.epiSnippetPlus.gAgeCat);    
    studyid = epiAnno.epiSnippetPlus.StudyId;
    %chop signal in selected IED-channel into epochs.
    epochLength = srate*2;
    NEpochs = floor(length(epiAnno.epiSnippetPlus.EEGIEDchannel(1:end)) / epochLength);
    powerDelta = zeros;
    powerTheta = zeros;
    powerAlpha = zeros;
    powerBeta = zeros;
    totalpower = zeros;
    RMStest = zeros;
    for m = 1 : (NEpochs-1)
        EpochStartpoint = (m-1)*epochLength+1;
        Epoch = (epiAnno.epiSnippetPlus.EEGIEDchannel(EpochStartpoint : (EpochStartpoint+epochLength-1)))';      
        powerDelta(m) = bandpower(Epoch, srate, [1 3.5]);
        powerTheta(m) = bandpower(Epoch, srate, [4 7.5]);
        powerAlpha(m) = bandpower(Epoch, srate, [8 11.5]);
        powerBeta(m)  = bandpower(Epoch, srate, [12 30]);
        totalpower(m) = bandpower(Epoch, srate, [1 30]);
        %totalpower(m) = sum( abs(Epoch).^2 ) / length(Epoch);
    end
    
    BackgroundPower.powerDelta = [BackgroundPower.powerDelta ; {powerDelta}];
    BackgroundPower.powerTheta = [BackgroundPower.powerTheta ; {powerTheta}];
    BackgroundPower.powerAlpha = [BackgroundPower.powerAlpha ; {powerAlpha}];
    BackgroundPower.powerBeta = [BackgroundPower.powerBeta ; {powerBeta}];
    BackgroundPower.totalpower = [BackgroundPower.totalpower ; {totalpower}];
    BackgroundPower.SRID = [BackgroundPower.SRID ; SearchResultID];
    BackgroundPower.gAgeCat = [BackgroundPower.gAgeCat ; gAgeCat];
    BackgroundPower.StudyId = [BackgroundPower.StudyId ; studyid];
end %files/snippets
end %users

savePlotPattern = fullfile(saveFolder, 'poweranalysis_study');        
            save(savePlotPattern, 'BackgroundPower');
varargout = {BackgroundPower};