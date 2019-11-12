function GotoEvent(seconds)
    disp('Navigating to event');    
    
    if exist('EEG', 'var')
        error('No EEG open in EEGLAB');
    end
    
    existingPlot = GetEeglabPlot();
    if isempty(existingPlot)
        disp('EEG plot not found, cannot goto event');
    else        
        %EEG = evalin('base','EEG');
        
        %Add a little buffer on the left side
        seconds2 = max(1,seconds-1);
        
        EPosition = findobj('tag','EPosition','parent',existingPlot); % ui handle
        if ~isempty(EPosition) 
            set(EPosition, 'String', seconds2);
            eegplot('drawp', 0, '', existingPlot);
        end                
    end
end
