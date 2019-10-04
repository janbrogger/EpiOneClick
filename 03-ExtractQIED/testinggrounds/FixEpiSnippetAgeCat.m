function FixEpiSnippetAgeCat(userid, searchresultid)
% IMPORTANT: Edit userid and search result id in WHERE-statement.
% Dependencies - ScorePipeline episnippets at location C:\Midlertidig_Lagring\Epileptiform\eegsnippets
%
% Inserts gAgeCat into Epileptiform\eegsnippets\User*EpiPlot*
%
% STATA gAgeCat needs to be inserted into EpiSnippets to be able to plot  
% the average spike across age groups.

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = 'C:\Midlertidig_Lagring\Epileptiform\eegsnippets';
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
filePattern = fullfile(rootFolder, allFolders(l).name, 'Epi_SearchResultEventId*.mat');
theFiles = dir(filePattern);

    %Loop through EEG-snippets
    for k = 1 : length(theFiles)    
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
        fprintf(1, 'Reading file %s \n', fullFileName);

        % Load episnippet
        epiAnno = load(fullFileName);
        sreID = epiAnno.epiSnippet.SearchResultEventId;        
        
        SQLquery = ([ 
        'SELECT DISTINCT StataPatLevelFocalEpi.gAgeCat ' ...
        'FROM SearchResult_Event ' ...        
        'INNER JOIN SearchResult_Event_Annotation ON SearchResult_Event.SearchResultEventId = SearchResult_Event_Annotation.SearchResultEventId ' ...
        'INNER JOIN Event ON SearchResult_Event.Eventid = Event.EventId  ' ...
        'INNER JOIN Recording ON Event.RecordingId = Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_Recording ON Recording.RecordingId = SearchResult_Recording.RecordingId  ' ...
        'INNER JOIN StataPatLevelFocalEpi ON Recording.StudyId = StataPatLevelFocalEpi.StudyId ' ...
        ' WHERE SearchResult_Event.SearchResultId = 12 AND SearchResult_Recording.SearchResultId = 12 AND SearchResult_Event_Annotation.UserId = 25']);

        SQLquery = strcat(SQLquery, ' AND SearchResult_Event.SearchResultEventId =', num2str(sreID));  

        gAgeCat = ScoreQueryRun(SQLquery);
        epiAnno.epiSnippet.gAgeCat = strtrim(cell2mat(gAgeCat{1,1}));        
        epiSnippet = epiAnno.epiSnippet;
        fileDir = ['C:\Midlertidig_Lagring\epileptiform\eegsnippets\User_' num2str(userid) '\'];
        fileName = ['Epi_SearchResultEventId_' num2str(sreID)];
        filePath = fullfile(fileDir, fileName);
        save(filePath, 'epiSnippet');   

    end %EEG-snippets
end %subfolders/users






