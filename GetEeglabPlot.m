function existingPlot = GetEeglabPlot(printWarning)
        
    if not(exist('printWarning','var'))
        printWarning = true;
    end

    existingPlot = findobj(0, 'tag', 'EEGPLOT');
    if isempty(existingPlot) && printWarning
        warning('No EEG plot open in EEGLAB');
    elseif size(existingPlot,1)>1        
        if printWarning
            warn('More than one EEGPLOT open, using the first one');        
        end
        existingPlot = existingPlot(1);        
    end    
end
