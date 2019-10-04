function ExtractEpiSnippets(online, doBandPass, doNotch, searchresultid, searchresulteventid, srIDOperator)

% varargin: 
% searchresultid:       ScorePipeline SearchResultId to extract annotated EEG-snippets from.
% doBandPass:           Apply bandpass filter.
% doNotch:              Apply notch filter.
% searchresulteventid:  Optional to specify SearchResultEventId or Ids to extract.
% srIDOperator          Operator in SQLquery when extracting target events.
% runtime:              0: Scorepipeline is not running. 
%                       1: Scorepipeline is running (do not try to open EEG file)
%
% varargout:  
%                       None.
%
% Dependencies - ScorePipeline, searchresult should be completely
% annotated, corresponding EEG-files must exist on disk.
%
% Due to errors when plotting in EEGLAB, the code from ScoreMouseDown
% has been internalized and plotting-code removed. This function breaks if 
% EEG-files are missing.
% Saves episnippets at location
% C:\Midlertidig_Lagring\Epileptiform\eegsnippets. The snippets are
% extended to include 2*samplerate before and after first and last
% annotation point.

%Get recordingIDs with annotations. 
%Recording is needed for filepath. 
%SearchResultEventID is needed as input to ScoreOpenEegFileInEEglab.m
%AnnotationConfig is needed to know the unit of ValueInt (which measure)
%Select order needs to correspond to variables firstcolumn, secondcolumn
%etc.

SQLquery = ([ 
        'SELECT SearchResultRecordingId, SearchResult_Event.SearchResultEventId, SearchResult_Event_Annotation.SearchResultAnnotationConfigId, SearchResult_Event_Annotation.ValueInt, SearchResult_Event_Annotation.UserId, Recording.StudyId ' ...
        'FROM SearchResult_Event ' ...        
        'INNER JOIN SearchResult_Event_Annotation ON SearchResult_Event.SearchResultEventId = SearchResult_Event_Annotation.SearchResultEventId ' ...
        'INNER JOIN Event ON SearchResult_Event.Eventid = Event.EventId  ' ...
        'INNER JOIN Recording ON Event.RecordingId = Recording.RecordingId  ' ...
        'INNER JOIN SearchResult_Recording ON Recording.RecordingId = SearchResult_Recording.RecordingId  ' ...
        ' WHERE SearchResult_Event.SearchResultId =  ' num2str(searchresultid) ' AND SearchResult_Recording.SearchResultId = ' num2str(searchresultid)]);
%Concatenate searchresulteventid to SQL query if specified as input argument
if(searchresulteventid > 0 && (strcmp(srIDOperator, '>') || strcmp(srIDOperator, '<') || strcmp(srIDOperator, '=')))
    SQLquery = strcat(SQLquery, ' AND SearchResult_Event.SearchResultEventId ', srIDOperator, num2str(searchresulteventid));  
end

recordingandeventids = ScoreQueryRun(SQLquery);  

%Loop through recordingIDs, then userIDs, and lastly eventIDs, to extract annotated episnippets
recordingandannotations = table2array(recordingandeventids);
[recordingstudyids] = unique(recordingandannotations(:,[1 6]), 'rows');
for l = 1 : length(recordingstudyids(:,1))
    recordingID = recordingstudyids(l,1);
    studyID = recordingstudyids(l,2);
    [fileExists, filePath] = ScoreCheckOneRecordingFile(recordingstudyids(l,1));
    if fileExists == -1
        fileExistsText = ['EEG FILE NOT FOUND. RecordingID: ' num2str(recordingID)];
        disp(fileExistsText)
        break;
    elseif fileExists == 1
        fileExistsText = ['EEG file was found. RecordingID: ' num2str(recordingID)];
    else
        fileExistsText = ['Unknown file status. RecordingID: ' num2str(recordingID)];
    end    
    
    disp(fileExistsText)
    
    %all annotations for one recording
    firstColumn = recordingandannotations(:,1);
    annotationsOneRecording = recordingandannotations(firstColumn==recordingID,:);
    
    if(online == 0)
        success = ScoreOpenEegFileInEeglabNOPLOT(filePath, doBandPass, doNotch, annotationsOneRecording(1,2));
        if(~success)
            break;
        end
    end
    EEG = evalin('base','EEG');
    
    %loop through each user
    uniqueusers = unique(annotationsOneRecording(:,5));
    for m = 1 : length(uniqueusers)
    userID = uniqueusers(m);
    fifthcolumn = annotationsOneRecording(:,5);
    annotationsOneUser = annotationsOneRecording(fifthcolumn==userID, :); 
        %loop through each annotated event in this recording
        uniqueevents = unique(annotationsOneUser(:,2));
        for n = 1 : length(uniqueevents)
            eventID = uniqueevents(n);
            secondColumn = annotationsOneUser(:,2);
            annotationsOneEvent = annotationsOneUser(secondColumn==eventID, :);  
            thirdColumn = annotationsOneEvent(:,3);
                        
            %Save a separate file with the snippet 
            startSample = max(1, annotationsOneEvent(thirdColumn==29,4)-EEG.srate*2); %21 spikestart, 4 column for valueint
            endSample   = min(size(EEG.data,2), annotationsOneEvent(thirdColumn==30,4)+EEG.srate*2); %24 slowwaveend, 4 column for valueint
            epiSnippet.spikeStart = annotationsOneEvent(thirdColumn==29,4) - startSample;
            epiSnippet.spikePeak = annotationsOneEvent(thirdColumn==30,4) - startSample; %22 spikePeak
            epiSnippet.spikeEnd = annotationsOneEvent(thirdColumn==31,4) - startSample; %23 spikeEnd
            epiSnippet.afterDischargeEnd = annotationsOneEvent(thirdColumn==32,4) - startSample; %24 afterdischargeEnd
            epiSnippet.EEGdata = EEG.data(:,startSample:endSample);
            txtnotch = 'notchOFF';
            txtbandpass = 'bandpassOFF';
            if(doNotch == 1)
                txtnotch = '';
            end
            if(doBandPass == 1)
                txtbandpass = '';
            end
            
            fileDir = ['C:\Midlertidig_Lagring\epileptiform\eegsnippets\User_' num2str(userID) '\'];
            fileName = ['Epi_SID_' num2str(searchresultid) 'SearchResultEventId_' num2str(eventID) txtnotch txtbandpass];
            filePath = fullfile(fileDir, fileName);
            epiSnippet.clickedChannelIndex = annotationsOneEvent(thirdColumn==34,4);
            epiSnippet.chanLocs = EEG.chanlocs;        
            epiSnippet.srate = EEG.srate;
            epiSnippet.bandPass = doBandPass;
            epiSnippet.notch = doNotch;
            epiSnippet.SearchResultEventId = eventID;
            epiSnippet.StudyId = studyID;
            save(filePath, 'epiSnippet');            
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