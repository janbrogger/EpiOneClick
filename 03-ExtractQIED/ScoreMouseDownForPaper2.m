
g = get(gcf,'UserData');
cp = get(gca,'currentpoint');
EEG = evalin('base','EEG');
isOneClick = 1;
[clickedTime, clickedSample, selectedChannelIndex, selectedChannel, clickedEegValue] = ScoreGetClickedTraceFromPoint(cp);    

%disp(['Clicked elapsed time : ' num2str(clickedTime)]);  
%disp(['Clicked sample : ' num2str(clickedSample)]);  
%disp(['Current time position at left is ' num2str(g.time)]);
%disp(['Current time position at right is ' num2str(g.time+g.winlength*g.srate)]);    
%disp([' Selected channel ' EEG.chanlocs(selectedChannelIndex).labels ' ,closest EEG value ' num2str(clickedEegValue) ]);    

if ~isfield(g, 'scoreAnnotationState')  || isempty(g.scoreAnnotationState)
    g.scoreAnnotationState = 'WaitingForFirstClick';        
end

if ~isempty(clickedTime)
    if strcmp(g.scoreAnnotationState, 'WaitingForFirstClick')
        g.scoreClickedChannelIndex= selectedChannelIndex;
        g.scoreClickedEegValue = clickedEegValue;
        g.selectedChannel = selectedChannel.labels;
        if(~isOneClick)
            g.scoreAnnotationState = 'WaitingForClick2';
            g.scoreClickedSample = [];
            g.scoreClickedTime = clickedTime;
            g.scoreClickedSample = [g.scoreClickedSample clickedSample];                
        else
            autoannotation = qIEDScorePipeline(selectedChannelIndex,clickedSample);
            g.scoreClickedSample = autoannotation{1};
            g.scoreClickedSample = [g.scoreClickedSample autoannotation{2}];
            g.scoreClickedSample = [g.scoreClickedSample autoannotation{3}];
            g.scoreClickedSample = [g.scoreClickedSample autoannotation{4}];
            g.scoreAnnotationState = 'Ready';
        end
    elseif strcmp(g.scoreAnnotationState, 'WaitingForClick2')
        g.scoreAnnotationState = 'WaitingForClick3';
        g.scoreClickedTime = clickedTime;
        g.scoreClickedSample = [g.scoreClickedSample clickedSample];

    elseif strcmp(g.scoreAnnotationState, 'WaitingForClick3')
        g.scoreAnnotationState = 'WaitingForClick4';
        g.scoreClickedTime = clickedTime;
        g.scoreClickedSample = [g.scoreClickedSample clickedSample];

    elseif strcmp(g.scoreAnnotationState, 'WaitingForClick4')
        g.scoreAnnotationState = 'Ready';  
        g.scoreClickedSample = [g.scoreClickedSample clickedSample];
    end
    if strcmp(g.scoreAnnotationState, 'Ready')  
        g.scoreAnnotationState = 'WaitingForFirstClick';  
        clickedSamples = g.scoreClickedSample;            
        ClickedChannelIndex = g.scoreClickedChannelIndex;
        ClickedChannel = g.selectedChannel;
        dataSegment = EEG.data(ClickedChannelIndex, clickedSamples(1):clickedSamples(2));
        [ymax ymaxsample] = max(dataSegment);
        allFigures = findall(0,'type','figure');
        oneEventDetails = findobj(allFigures, 'tag', 'oneEventDetails');
        if ~isempty(oneEventDetails)                
            handles = guidata(oneEventDetails);
            if isfield(handles, 'SearchResultEventId')                    
                configData = ScoreGetAnnotationsForOneProject(handles.SearchResultId);
                if strcmp(configData,'No Data')
                    warning(['No annotation configurations have been set up']);
                else
                    OtherChannels = ScoreGetOneEventLocationTextInfo(handles.SearchResultEventId);  
                    spikeStartAID = [];
                    spikeCenterAID = [];
                    spikeEndAID = [];
                    afterDischargeEndAID = [];
                    ClickedChannelAID = [];
                    ClickedChannelIndexAID = [];
                    OtherChannelsAID = [];
                    configData = ScoreGetAnnotationsForOneProject(handles.SearchResultId);
                    if ~strcmp(configData,'No Data')
                        for i=1:size(configData,1)
                            if strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'Start')
                                spikeStartAID = configData.SearchResultAnnotationConfigId(i);
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'Center')
                                spikeCenterAID = configData.SearchResultAnnotationConfigId(i);  
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'End')
                                spikeEndAID = configData.SearchResultAnnotationConfigId(i); 
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'SlowEnd')
                                afterDischargeEndAID = configData.SearchResultAnnotationConfigId(i); 
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'ClickedChannel')
                                ClickedChannelAID = configData.SearchResultAnnotationConfigId(i);   
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'ClickedChannelIndex')
                                ClickedChannelIndexAID = configData.SearchResultAnnotationConfigId(i); 
                            elseif strcmp(configData.AnnotationLevel(i), 'Event') ...
                                && strcmp(configData.FieldName(i), 'OtherChannels')
                                OtherChannelsAID = configData.SearchResultAnnotationConfigId(i); 
                            end
                        end
                    end 
                    
                    if ~isempty(spikeStartAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, spikeStartAID,'ValueInt',clickedSamples(1));
                    else
                        warning(['Annotation configuration SpikeStart not found']);
                    end
                    if ~isempty(spikeCenterAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, spikeCenterAID,'ValueInt',clickedSamples(2));
                    else
                        warning(['Annotation configuration SpikeCenter not found']);
                    end  
                    if ~isempty(spikeEndAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, spikeEndAID,'ValueInt',clickedSamples(3));
                    else
                        warning(['Annotation configuration SpikeEnd not found']);
                    end 
                    if ~isempty(afterDischargeEndAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, afterDischargeEndAID,'ValueInt',clickedSamples(4));
                    else
                        warning(['Annotation configuration AfterDischargeEnd not found']);
                    end
                    if ~isempty(ClickedChannelAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, ClickedChannelAID,'ValueText',ClickedChannel);
                    else
                        warning(['Annotation configuration ClickedChannel not found']);
                    end  
                    if ~isempty(OtherChannelsAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, OtherChannelsAID,'ValueText',OtherChannels);
                    else
                        warning(['Annotation configuration ClickedChannel not found']);
                    end  
                    if ~isempty(ClickedChannelIndexAID)
                        ScoreSetAnnotationForOneEvent(handles.SearchResultEventId, ClickedChannelIndexAID,'ValueInt',ClickedChannelIndex);
                    else
                        warning(['Annotation configuration ClickedChannel not found']);
                    end  
                end  
            end
            handles.UpdateCustomAnnotations(handles);
        end

    end
end    
set(gcf, 'UserData', g);  



