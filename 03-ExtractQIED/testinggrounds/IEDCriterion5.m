function varargout = IEDCriterion5()

% varargin:
%           - None
% varargout:
%           - Quantitative measures for IED Criterion5 per User.
%
% Dependencies - ScorePipeline episnippets at location C:\Midlertidig_Lagring\Epileptiform\eegsnippets


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
Criterion5.epochOffsets = [];
Criterion5.epochdifferences = [];
Criterion5.powerDiffDelta = [];
Criterion5.powerDiffTheta = [];
Criterion5.powerDiffAlpha = [];
Criterion5.powerDiffBeta = [];
Criterion5.SRID = [];
Criterion5.IEDBGDeltaDiff = [];
Criterion5.IEDBGThetaDiff = [];
Criterion5.IEDBGAlphaDiff = [];
Criterion5.IEDBGBetaDiff = [];

% Loop through subfolders (users)
for l = 1 : length(allFolders)
filePattern = fullfile(rootFolder, allFolders(l).name, 'EpiPlus_SearchResultEventId*.mat');
theFiles = dir(filePattern);

%Loop through EEG-snippets
for k = 1 : length(theFiles)    
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
    fprintf(1, 'Reading file %s \n', fullFileName);
    
    % Load episnippet signal from one channel.
    % Measure total power difference in consecutive E(x) and E(x+1), where  
    % represents one second samples. Ex is spaced by IED-length in samples.
    epiAnno = load(fullFileName);
    srate = epiAnno.epiSnippetPlus.srate;
    IEDlength = epiAnno.epiSnippetPlus.afterDischargeEnd - epiAnno.epiSnippetPlus.spikeStart;
    SearchResultID = epiAnno.epiSnippetPlus.SearchResultEventId;
    snippet = epiAnno.epiSnippetPlus.EEGdata(epiAnno.epiSnippetPlus.clickedChannelIndex, :);
    preBGstart = epiAnno.epiSnippetPlus.spikeStart - (1*srate) - 1;
    postBGend = epiAnno.epiSnippetPlus.afterDischargeEnd+1+srate;
    if(preBGstart >= 1 && postBGend <= length(snippet))
        preBG = snippet(preBGstart : epiAnno.epiSnippetPlus.spikeStart - 1);
        postBG = snippet(epiAnno.epiSnippetPlus.afterDischargeEnd+1 : postBGend );        
        BGDeltadiff = bandpower(postBG, srate, [1 3.5]) - bandpower(preBG, srate, [1 3.5]);
        BGThetadiff = bandpower(postBG, srate, [4 7.5]) - bandpower(preBG, srate, [4 7.5]);
        BGAlphadiff = bandpower(postBG, srate, [8 11.5]) - bandpower(preBG, srate, [8 11.5]);
        BGBetadiff  = bandpower(postBG, srate, [12 20]) - bandpower(preBG, srate, [12 20]);
        %snippetfig = figure();
        %plot([preBG postBG]);
        %close(snippetfig);
    else
        BGDeltadiff = NaN;
    end
    %chop signal in selected IED-channel into epochs "synched" to IED.
    epochLength = srate + IEDlength;
    IEDposition = epiAnno.epiSnippetPlus.positionStart;
    epochOffset = mod(IEDposition, epochLength);
    NEpochs = floor(length(epiAnno.epiSnippetPlus.EEGIEDchannel(1+epochOffset:end)) / epochLength);
    epochs = zeros(NEpochs, srate);
    epochdiff = zeros;
    powerDiffDelta = zeros;
    powerDiffTheta = zeros;
    powerDiffAlpha = zeros;
    powerDiffBeta = zeros;
    for m = 1 : (NEpochs-1)
        E1startpoint = (m-1)*epochLength+1+epochOffset;
        E2startpoint = E1startpoint+epochLength;
        epochs(m,:) = epiAnno.epiSnippetPlus.EEGIEDchannel(E1startpoint : (E1startpoint+srate-1));
        E1 = double( epochs(m,:) );
        E2 = double( epiAnno.epiSnippetPlus.EEGIEDchannel(E2startpoint :(E2startpoint+srate-1)) );
        powerDiffDelta(m) = bandpower(E2, srate, [1 3.5]) - bandpower(E1, srate, [1 3.5]);
        powerDiffTheta(m) = bandpower(E2, srate, [4 7.5]) - bandpower(E1, srate, [4 7.5]);
        powerDiffAlpha(m) = bandpower(E2, srate, [8 11.5]) - bandpower(E1, srate, [8 11.5]);
        powerDiffBeta(m)  = bandpower(E2, srate, [12 20]) - bandpower(E1, srate, [12 20]);
        E1power = sum( abs(E1).^2 ) / length(E1);
        E2power = sum( abs(E2).^2 ) / length(E2);
        epochdiff(m) = E2power-E1power;
        %if abs(epochdiff(m)) > 40000000
        %    epochfig = figure();
        %    plot([E1 E2])
        %    waitfor(epochfig);
        %end
    end
    Criterion5.epochOffsets = [Criterion5.epochOffsets ; epochOffset];
    Criterion5.epochdifferences = [Criterion5.epochdifferences ; {epochdiff}];
    Criterion5.powerDiffDelta = [Criterion5.powerDiffDelta ; {powerDiffDelta}];
    Criterion5.powerDiffTheta = [Criterion5.powerDiffTheta ; {powerDiffTheta}];
    Criterion5.powerDiffAlpha = [Criterion5.powerDiffAlpha ; {powerDiffAlpha}];
    Criterion5.powerDiffBeta = [Criterion5.powerDiffBeta ; {powerDiffBeta}];
    Criterion5.SRID = [Criterion5.SRID ; SearchResultID];
    Criterion5.IEDBGDeltaDiff = [Criterion5.IEDBGDeltaDiff ; BGDeltadiff];
    Criterion5.IEDBGThetaDiff = [Criterion5.IEDBGThetaDiff ; BGThetadiff];
    Criterion5.IEDBGAlphaDiff = [Criterion5.IEDBGAlphaDiff ; BGAlphadiff];
    Criterion5.IEDBGBetaDiff = [Criterion5.IEDBGBetaDiff ; BGBetadiff];
end %files/snippets
end %users
varargout = {Criterion5}