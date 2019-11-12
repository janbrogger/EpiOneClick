function [clickedTime, clickedSample, selectedChannelIndex, selectedChannel, clickedEegValue] = ScoreGetClickedTraceFromPoint(point)    
    clickedTime = [];
    clickedSample = [];
    selectedChannelIndex = [];
    selectedChannel = [];
    clickedEegValue = [];
    
    g = get(gcf,'UserData');
    ax1 = findobj('tag', 'eegaxis', 'parent', gcf);
    XLim = get(ax1, 'XLim');
    YLim = get(ax1, 'YLim');    
    
    if isfield(g, 'time') &&  ...
        point(1,1) > 0 && ...        
        point(1,1) <= g.winlength*g.srate && ...
        point(1,2) >= YLim(1) && ...
        point(1,2) <= YLim(2) 
        
        clickedTime = g.time+point(1)/g.srate;
        clickedSample = max(round(g.time*g.srate, 0), round(g.time*g.srate+point(1),0));
        if clickedSample == 0
            clickedSample = 1;
        end

        clickedYValue = point(1,2);        
        %disp(['Clicked EEG value: ' num2str(clickedYValue)]);

        EEG = evalin('base','EEG');
        thisColumnOfData = EEG.data(:,clickedSample);    

        g = get(gcf,'UserData');        
        for i=1:size(thisColumnOfData, 1)
            thisColumnOfDataWithSpacing(i) = EEG.data(i,clickedSample)+g.spacing*(1+size(thisColumnOfData, 1)-i);
        end    

        diffThisColumnOfDataWithSpacing = thisColumnOfDataWithSpacing-clickedYValue;    
        [diffValue,indexOfClosestMatch] = min(abs(diffThisColumnOfDataWithSpacing));
        selectedChannel = EEG.chanlocs(indexOfClosestMatch);

        selectedChannelIndex = indexOfClosestMatch;
        clickedEegValue = EEG.data(indexOfClosestMatch,clickedSample);        

        %disp(['Clicked elapsed time : ' num2str(clickedTime)]);  
        %disp(['Clicked sample : ' num2str(clickedSample)]);  
        %disp(['Current time position at left is ' num2str(g.time)]);
        %disp(['Current time position at right is ' num2str(g.time+g.winlength*g.srate)]);
        %disp(['Clicked voltage: ' num2str(clickedYValue)]);      
        %disp([' Selected channel ' selectedChannel.labels ' ,closest EEG value ' num2str(EEG.data(indexOfClosestMatch,clickedSample)) ]);    
    end
    
end