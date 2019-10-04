function SnippetExtract(userid, searchresultid, searchresulteventid, srIDOperator, doBandPass, doNotch)
%
% searchresultid:       ScorePipeline SearchResultId to extract annotated IED or transient.
% userid:               Your userid in ScorePipeline
% searchresulteventid:  0: Ignored. >0:SearchResultEventId or Ids to extract.
% srIDOperator          Operator in SQLquery when extracting specific events.
% doBandPass:           Apply bandpass filter.
% doNotch:              Apply notch filter.
%
% IMPORTANT:            Make sure variable rootfolder is correct
%
%
% out:  
%                       EEG-snippets at location C:\Midlertidig_Lagring\epileptiform\1b-Matlab\data\User_?\snippets\ .
%
% Dependencies - ScorePipeline, searchresult should be completely
% annotated, corresponding EEG-files must exist with correct paths, folder 
% C:\Midlertidig_Lagring\epileptiform\1b-Matlab\data\User_* (where integrated matlab-analyses are saved), 
% and further subfolder ...User_*\snippets, where the snippet m-files are saved. 
% Annotations should have >= 2*sample rate before and after IED-complex.

% Due to errors when plotting in EEGLAB, the code from ScoreMouseDown
% has been internalized and plotting-code removed. This function breaks if 
% EEG-files are missing.

%Get recordingIDs with annotations. 
%Recording is needed for filepath. 
%SearchResultEventID is needed as input to ScoreOpenEegFileInEEglab.m
%AnnotationConfig is needed to know the unit of annotated value (which measure)
%Select order needs to correspond to variables firstcolumn, secondcolumn
%etc.

rootfolder = '\\ihelse.net\Forskning\hbe\2017-01512\2019\EivindArtikkel2\data\03-ExtractQIED';

%Get IDs from already extracted snippets to skip them
fileDir2 = [rootfolder 'User_' num2str(userid) '\snippets\'];    
filePattern2 = fullfile(fileDir2, 'snippet*.mat');
oldFiles = dir(filePattern2);
completed = zeros(length(oldFiles),2);
%Loop through old EEG-snippets
for k = 1 : length(oldFiles)
    baseFileName = oldFiles(k).name;
    fullFileName = fullfile(fileDir2, baseFileName);
    %fprintf(1, 'Reading file %s \n', fullFileName);
    % Load episnippet and save userid/sreid information
    existingsnippet = load(fullFileName);
    oldstudyid = existingsnippet.snippet.StudyId;
    olduserid = existingsnippet.snippet.UserID;
    completed(k,:) = [olduserid oldstudyid];
end

SQLqueryTransients = ([ 
        'SELECT SearchResultRecordingId, SearchResult_Event.SearchResultEventId, SearchResult_Event_Annotation.SearchResultAnnotationConfigId, ' ... 
        'CASE ' ...
        'WHEN SearchResult_AnnotationConfig.HasInteger=1 THEN CAST(SearchResult_Event_Annotation.ValueInt AS varchar(255)) ' ...
        'WHEN SearchResult_AnnotationConfig.HasFloat=1 THEN CAST(SearchResult_Event_Annotation.ValueFloat AS varchar(255)) ' ...
        'WHEN SearchResult_AnnotationConfig.HasString=1 THEN CAST(SearchResult_Event_Annotation.ValueText AS varchar(255))  ' ...
        'WHEN SearchResult_AnnotationConfig.HasBit=1 THEN CAST(SearchResult_Event_Annotation.ValueBit AS varchar(255))  ' ...
        'WHEN SearchResult_AnnotationConfig.HasBlob=1 THEN CAST(SearchResult_Event_Annotation.ValueBlob AS varchar(255))  ' ...
        'ELSE ''unknown'' END ' ...
        'AS Value, ' ...
        'SearchResult_Event_Annotation.UserId, Recording.StudyId, StataPatLevelSharpTransients.gAgeCat ' ...
        'FROM SearchResult_Event ' ...        
        'INNER JOIN SearchResult_Event_Annotation ON SearchResult_Event.SearchResultEventId = SearchResult_Event_Annotation.SearchResultEventId ' ...
        'INNER JOIN Event ON SearchResult_Event.Eventid = Event.EventId  ' ...
        'INNER JOIN Recording ON Event.RecordingId = Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_Recording ON Recording.RecordingId = SearchResult_Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_AnnotationConfig ON SearchResult_Event_Annotation.SearchResultAnnotationConfigId = SearchResult_AnnotationConfig.SearchResultAnnotationConfigId ' ...
        'INNER JOIN StataPatLevelSharpTransients ON Recording.StudyId = StataPatLevelSharpTransients.StudyId ' ...
        ' WHERE SearchResult_Event.SearchResultId =  ' num2str(searchresultid) ' AND SearchResult_Recording.SearchResultId = ' num2str(searchresultid) ' AND SearchResult_Event_Annotation.UserId = ' num2str(userid)]);
%Concatenate searchresulteventid to SQL query if specified as input argument
if(searchresulteventid > 0 && (strcmp(srIDOperator, '>') || strcmp(srIDOperator, '<') || strcmp(srIDOperator, '=')))
    SQLqueryTransients = strcat(SQLqueryTransients, ' AND SearchResult_Event.SearchResultEventId ', srIDOperator, num2str(searchresulteventid));  
end

recordingandeventids = ScoreQueryRun(SQLqueryTransients); 
m = size(recordingandeventids,1);
A = table(zeros(m,1));
recordingandeventids = [recordingandeventids A];
SQLqueryIEDs = ([ 
        'SELECT SearchResultRecordingId, SearchResult_Event.SearchResultEventId, SearchResult_Event_Annotation.SearchResultAnnotationConfigId, ' ... 
        'CASE ' ...
        'WHEN SearchResult_AnnotationConfig.HasInteger=1 THEN CAST(SearchResult_Event_Annotation.ValueInt AS varchar(255)) ' ...
        'WHEN SearchResult_AnnotationConfig.HasFloat=1 THEN CAST(SearchResult_Event_Annotation.ValueFloat AS varchar(255)) ' ...
        'WHEN SearchResult_AnnotationConfig.HasString=1 THEN CAST(SearchResult_Event_Annotation.ValueText AS varchar(255))  ' ...
        'WHEN SearchResult_AnnotationConfig.HasBit=1 THEN CAST(SearchResult_Event_Annotation.ValueBit AS varchar(255))  ' ...
        'WHEN SearchResult_AnnotationConfig.HasBlob=1 THEN CAST(SearchResult_Event_Annotation.ValueBlob AS varchar(255))  ' ...
        'ELSE ''unknown'' END ' ...
        'AS Value, ' ...
        'SearchResult_Event_Annotation.UserId, Recording.StudyId, StataPatLevelFocalEpi.gAgeCat ' ...
        'FROM SearchResult_Event ' ...        
        'INNER JOIN SearchResult_Event_Annotation ON SearchResult_Event.SearchResultEventId = SearchResult_Event_Annotation.SearchResultEventId ' ...
        'INNER JOIN Event ON SearchResult_Event.Eventid = Event.EventId  ' ...
        'INNER JOIN Recording ON Event.RecordingId = Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_Recording ON Recording.RecordingId = SearchResult_Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_AnnotationConfig ON SearchResult_Event_Annotation.SearchResultAnnotationConfigId = SearchResult_AnnotationConfig.SearchResultAnnotationConfigId ' ...
        'INNER JOIN StataPatLevelFocalEpi ON Recording.StudyId = StataPatLevelFocalEpi.StudyId ' ...
        ' WHERE SearchResult_Event.SearchResultId =  ' num2str(searchresultid) ' AND SearchResult_Recording.SearchResultId = ' num2str(searchresultid) ' AND SearchResult_Event_Annotation.UserId = ' num2str(userid)]);
%Concatenate searchresulteventid to SQL query if specified as input argument
if(searchresulteventid > 0 && (strcmp(srIDOperator, '>') || strcmp(srIDOperator, '<') || strcmp(srIDOperator, '=')))
    SQLqueryIEDs = strcat(SQLqueryIEDs, ' AND SearchResult_Event.SearchResultEventId ', srIDOperator, num2str(searchresulteventid));  
end
recordingandeventidsIED = ScoreQueryRun(SQLqueryIEDs); 
m = size(recordingandeventidsIED,1);
A = table(ones(m,1));
recordingandeventidsIED = [recordingandeventidsIED A];

recordingandeventids = [recordingandeventids; recordingandeventidsIED]; 


%Loop through recordingIDs, then userIDs, and lastly eventIDs, to extract annotated snippets
values = table2cell(recordingandeventids(:, [3 4]));
recordingandannotations = table2cell(recordingandeventids);
[recordingstudyids] = unique(cell2mat(recordingandannotations(:,[1 6])), 'rows');
for l = 1 : length(recordingstudyids(:,1))
    recordingID = recordingstudyids(l,1);
    studyID = recordingstudyids(l,2);
    %skip this event/snippet if already extracted
    if(ismember(studyID,completed(:,2)))
        continue;
    end
    [fileExists, filePath] = ScoreCheckOneRecordingFile(recordingstudyids(l,1));
    if fileExists == -1
        fileExistsText = ['EEG FILE NOT FOUND. RecordingID: ' num2str(recordingID)];
        disp(fileExistsText)
        continue;
    elseif fileExists == 1
        fileExistsText = ['EEG file was found. RecordingID: ' num2str(recordingID)];
    else
        fileExistsText = ['Unknown file status. RecordingID: ' num2str(recordingID)];
    end    
    
    disp(fileExistsText)
    
    %all annotations for one recording
    firstColumn = [recordingandannotations{:,1}];
    annotationsOneRecording = recordingandannotations(firstColumn==recordingID,:);
    success = ScoreOpenEegFileInEeglabNOPLOT(filePath, doBandPass, doNotch, annotationsOneRecording{1,2});
    if(~success)
        break;
    end
    EEG = evalin('base','EEG');
    
    %loop through each user
    uniqueusers = unique(cell2mat(annotationsOneRecording(:,5)));
    uniquegAgeCat = strtrim(annotationsOneRecording{1,7}); %getting gAgeCat too (25.06.19)
    isIED = annotationsOneRecording{1,8}; %getting if snippet is IED too (27.06.19)
    for m = 1 : length(uniqueusers)
    userID = uniqueusers(m);    
    fifthcolumn = [annotationsOneRecording{:,5}];
    annotationsOneUser = annotationsOneRecording(fifthcolumn==userID, :); 
        %loop through each annotated event in this recording
        uniqueevents = unique(cell2mat(annotationsOneUser(:,2)));
        for n = 1 : length(uniqueevents)
            eventID = uniqueevents(n);            
            secondColumn = [annotationsOneUser{:,2}];
            annotationsOneEvent = annotationsOneUser(secondColumn==eventID, :);  
            thirdColumn = [annotationsOneEvent{:,3}];
                        
            %Save a separate file with the snippet 
            startSample = max(1, str2double(annotationsOneEvent(thirdColumn==43,4))-EEG.srate*4-1); %21 spikestart, 4 column for valueint
            endSample   = min(size(EEG.data,2), str2double(annotationsOneEvent(thirdColumn==46,4))+EEG.srate*4+1); %24 slowwaveend, 4 column for valueint
            snippet.spikeStart = str2double(annotationsOneEvent(thirdColumn==43,4)) - startSample;
            snippet.spikePeak = str2double(annotationsOneEvent(thirdColumn==44,4)) - startSample; %22 spikePeak
            snippet.positionStart = str2double(annotationsOneEvent(thirdColumn==43,4)); %position of spike in the whole EEG
            snippet.spikeEnd = str2double(annotationsOneEvent(thirdColumn==45,4)) - startSample; %23 spikeEnd
            snippet.afterDischargeEnd = str2double(annotationsOneEvent(thirdColumn==46,4)) - startSample; %24 afterdischargeEnd
            snippet.EEGdata = EEG.data(:,startSample:endSample);            
            snippet.clickedChannelIndex = str2double(annotationsOneEvent(thirdColumn==48,4));
            snippet.clickedChannel = [annotationsOneEvent{thirdColumn==49,4}];
            snippet.localization = [annotationsOneEvent{thirdColumn==50,4}];
            snippet.halford = str2double(annotationsOneEvent(thirdColumn==47,4));
            %snippet.EEGIEDchannel = EEG.data(snippet.clickedChannelIndex, :);
            snippet.chanLocs = EEG.chanlocs;        
            snippet.srate = EEG.srate;
            snippet.bandPass = doBandPass;
            snippet.notch = doNotch;
            snippet.SearchResultEventId = eventID;
            snippet.StudyId = studyID;
            snippet.UserID = userID;
            snippet.gAgeCat = uniquegAgeCat;
            snippet.isIED = isIED;
            txtnotch = 'notchOFF';
            txtbandpass = 'bandpassOFF';
            if(doNotch == 1)
                txtnotch = '';
            end
            if(doBandPass == 1)
                txtbandpass = '';
            end
                        
            fileName = ['snippet_U_' num2str(userID) '_SRId_' num2str(searchresultid) '_SREId_' num2str(eventID) txtnotch txtbandpass];
            filePath = fullfile(fileDir2, fileName);            
            save(filePath, 'snippet');            
        end %events
    end %users
end %recordings    
end

%MATLAB keeps crashing. I suspect it is related to plotting.
%This function is ScoreOpenEegFileInEeglab with plotting-code removed.
function openSuccess = ScoreOpenEegFileInEeglabNOPLOT(newFilePath, doBandPass, doNotch, searchResultEventId)
    disp(['Opening new EEG file ' newFilePath ' SearchResultEventId: ' searchResultEventId] );    
    openSuccess = 0;
    %Start up EEGLAB if pop_fileio is not in the path   
    if not(exist('pop_fileio', 'file'))
        eeglab
    end
    %Close any existing plot
    existingPlot = findobj(0, 'tag', 'EEGPLOT');
    
    if size(existingPlot,1) == 1
        close(existingPlot.Number)
    elseif size(existingPlot,1) > 1
        for i = 1:size(existingPlot,1)    
            close(existingPlot(i).Number)
        end
    end
    
    if not(exist(newFilePath, 'file'))
       ScoreClearEeglabStudy()
    else        
        try 
            EEG = pop_fileio(newFilePath);                
            EEG.setname='test';    
            
            
            disp('Removing unused channels');
            EEG = pop_select( EEG,'nochannel',{'Rate' 'IBI' 'Bursts' 'Suppr'});                                                   
            exclIndex1 = find(strcmp({EEG.chanlocs(:).labels}, '24'));
            exclIndex2 = find(strcmp({EEG.chanlocs(:).labels}, '25'));
            exclIndex3 = find(strcmp({EEG.chanlocs(:).labels}, '26'));
            exclIndex4 = find(strcmp({EEG.chanlocs(:).labels}, '27'));
            exclIndex5 = find(strcmp({EEG.chanlocs(:).labels}, '28'));
            exclIndex6 = find(strcmp({EEG.chanlocs(:).labels}, '29'));
            exclIndex7 = find(strcmp({EEG.chanlocs(:).labels}, '30'));
            exclIndex8 = find(strcmp({EEG.chanlocs(:).labels}, '31'));
            EEG = pop_select( EEG,'nochannel',[exclIndex1  exclIndex2  exclIndex3  exclIndex4  exclIndex5  exclIndex6  exclIndex7  exclIndex8 ]);
            EEG = eeg_checkset( EEG );

            ekgindex = find(strcmp({EEG.chanlocs(:).labels}, 'EKG'));
            photicindex = find(strcmp({EEG.chanlocs(:).labels}, 'Photic'));
            disp('Flipping data');
            EEG.data = -EEG.data;
            disp('Rereferencing data');
            EEG = pop_reref( EEG, [],'exclude',[ekgindex photicindex ]);  
            %notch
            if(doNotch == 1)
                disp('Filtering data');
                EEG = pop_eegfilt(EEG, 48, 52, 3300, 1, [], 0);
            end
            %passband
            if(doBandPass == 1)
                EEG = pop_eegfilt(EEG, 1, 70, 6600, 0, [], 0);
            end
            disp('Fixing some events');
            EEG = SetSomeLongEventsToZero(EEG);
            EEG = InsertSomeEventNames(EEG);
            
            %Downscale EKG
            disp('Downscaling EKG');
            EEG.data(ekgindex,:) = EEG.data(ekgindex,:)/5;
            
            %Order columns
            montageOrder = {'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', ...
                'f11', 'F12', 'C3', 'C4', 'T7', 'T8', 'P3', 'P4', ...
                'P7', 'P8', 'P11', 'P12', 'O1', 'O2', 'A1', 'A2' ...
                'Fz', 'Cz', 'EKG', 'Photic'};
            %EEG = ScoreOrderChannels(EEG, montageOrder);
            
            clear ekgindex
            EEG.filename = ['searchResultEventId = ' num2str(searchResultEventId)];
            assignin('base','EEG',EEG);
            eeg_checkset( EEG );   
            openSuccess = 1;
        catch ME
            ScoreClearEeglabStudy()
        end
    end
end

function ScoreClearEeglabStudy()
        disp('Something went wrong...');
        assignin('base','STUDY', []); 
        assignin('base','CURRENTSTUDY', 0); 
        assignin('base','ALLEEG', []); 
        assignin('base','EEG', []); 
        assignin('base','CURRENTSET', []);  
end

function EEG = SetSomeLongEventsToZero(EEG)
        eventsToChangeDuration = find(strcmp({EEG.event(:).value},'Review progress'));
        for i = eventsToChangeDuration
            EEG.event(i).duration = 0;
        end
end

function EEG = InsertSomeEventNames(EEG)
    for i=1:size(EEG.event,2)
        for j=1:size(ScoreConfig.eventGuids)  
            if length(EEG.event(i).type) >= 36
                if strcmp(EEG.event(i).type(1:36), ScoreConfig.eventGuids(j,1))
                    EEG.event(i).value = ScoreConfig.eventGuids{j,2};
                    EEG.event(i).type = ScoreConfig.eventGuids{j,2};
                end                
            end
        end
    end
end