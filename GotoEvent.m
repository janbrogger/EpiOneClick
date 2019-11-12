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
        seconds2 = min(1,seconds-1);
        
        EPosition = findobj('tag','EPosition','parent',existingPlot); % ui handle
        set(EPosition, 'string', num2str(seconds2));

        %evalin('base','eegplot(''drawp'', 0);');    
        eegplot('drawp', 0, '', existingPlot);        
    end
end
