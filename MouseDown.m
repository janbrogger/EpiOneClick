function MouseDown(varargin)
    g = get(gcf,'UserData');
    cp = get(gca,'currentpoint');        
    [clickedTime, clickedSample, selectedChannelIndex, selectedChannel, clickedEegValue] = GetClickedTraceFromPoint(cp);    

    %disp(['Clicked elapsed time : ' num2str(clickedTime)]);  
    %disp(['Clicked sample : ' num2str(clickedSample)]);  
    %disp(['Current time position at left is ' num2str(g.time)]);
    %disp(['Current time position at right is ' num2str(g.time+g.winlength*g.srate)]);    
    %disp([' Selected channel ' EEG.chanlocs(selectedChannelIndex).labels ' ,closest EEG value ' num2str(clickedEegValue) ]);    

    if ~isfield(g, 'scoreAnnotationState')  || isempty(g.annotationState)
        g.annotationState = 'WaitingForFirstClick';        
    end

    if ~isempty(clickedTime)
        if strcmp(g.annotationState, 'WaitingForFirstClick')
            g.scoreClickedChannelIndex= selectedChannelIndex;
            g.scoreClickedEegValue = clickedEegValue;
            g.selectedChannel = selectedChannel.labels;
            
            autoannotation = qIEDAnnotation(selectedChannelIndex,clickedSample);
%             disp('Annotation results:');
%             disp(autoannotation);
%             g.scoreClickedSample = autoannotation{1};
%             g.scoreClickedSample = [g.scoreClickedSample autoannotation{2}];
%             g.scoreClickedSample = [g.scoreClickedSample autoannotation{3}];
%             g.scoreClickedSample = [g.scoreClickedSample autoannotation{4}];
            g.annotationState = 'WaitingForFirstClick';            
        end                     
    end    
    set(gcf, 'UserData', g);  

end