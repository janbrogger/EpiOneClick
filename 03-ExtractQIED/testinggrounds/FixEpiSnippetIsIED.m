function FixEpiSnippetIsIED(userid)
% IMPORTANT: Edit userid and search result id in WHERE-statement.
% Dependencies - ScorePipeline episnippets at location C:\Midlertidig_Lagring\epileptiform\eegsnippetsplus\Paper2
%
% Inserts isIED into
% C:\Midlertidig_Lagring\epileptiform\eegsnippetsplus\Paper2\User*EpiPlot*
% This functionality will be implemented in extractepichannel.m after
% finished annotation and then this file can be deleted. Until then it
% serves as a quick fix.

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = 'C:\Midlertidig_Lagring\epileptiform\eegsnippetsplus\Paper2\';
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

% Loop through subfolders (users)
for l = 1 : length(allFolders)
filePattern = fullfile(rootFolder, allFolders(l).name, 'snippet_U_25*.mat');
theFiles = dir(filePattern);

    %Loop through EEG-snippets
    for k = 1 : length(theFiles)    
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
        %fprintf(1, 'Reading file %s \n', fullFileName);

        % Load episnippet
        epiAnno = load(fullFileName);
        sreID = epiAnno.epiSnippetPlus.SearchResultEventId;        
        
        SQLquery = ([ 
        'SELECT DISTINCT StataPatLevelSharpTransients.gAgeCat ' ...
        'FROM SearchResult_Event ' ...        
        'INNER JOIN SearchResult_Event_Annotation ON SearchResult_Event.SearchResultEventId = SearchResult_Event_Annotation.SearchResultEventId ' ...
        'INNER JOIN Event ON SearchResult_Event.Eventid = Event.EventId  ' ...
        'INNER JOIN Recording ON Event.RecordingId = Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_Recording ON Recording.RecordingId = SearchResult_Recording.RecordingId  ' ...
        'INNER JOIN StataPatLevelSharpTransients ON Recording.StudyId = StataPatLevelSharpTransients.StudyId ' ...
        ' WHERE SearchResult_Event.SearchResultId = 18 AND SearchResult_Recording.SearchResultId = 18 AND SearchResult_Event_Annotation.UserId = 25']);

        SQLquery = strcat(SQLquery, ' AND SearchResult_Event.SearchResultEventId =', num2str(sreID));  

        isIED = ScoreQueryRun(SQLquery); %selects from sharptransients -> no hit means IED.
        if(strcmp(isIED, 'No Data'))
            epiAnno.epiSnippetPlus.isIED = 1;
        else
            epiAnno.epiSnippetPlus.isIED = 0;
        end % if SQL return something it is a sharp transient
        epiSnippetPlus = epiAnno.epiSnippetPlus;
        save(fullFileName, 'epiSnippetPlus');  
    end %EEG-snippets
end %subfolders/users






