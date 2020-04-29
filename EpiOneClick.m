function EpiOneClick(filePath)

    disp(['Opening new EEG file ' filePath] );    
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
    
    if not(exist(filePath, 'file'))
        assignin('base','STUDY', []); 
        assignin('base','CURRENTSTUDY', 0); 
        assignin('base','ALLEEG', []); 
        assignin('base','EEG', []); 
        assignin('base','CURRENTSET', []);  
    else        
        try 
            EEG = pop_fileio(filePath);                
            EEG.setname=filePath;    
                        
            disp('Locating empty channels');
            exclIndex = [];
            for i=1:size(EEG.chanlocs,1)
                stddev = std(EEG.data(i,:));
                if stddev < 1
                    exclIndex = [exclIndex i];
                end
            end
            if ~isempty(exclIndex)            
                EEG = pop_select( EEG,'nochannel',exclIndex);
                disp(['Empty channels dropped:' exclIndex]);
            end
           
            EEG = eeg_checkset( EEG );

            disp('Looking for EKG and photic...');
            ekgindex = find(strcmp({EEG.chanlocs(:).labels}, 'EKG'));
            photicindex = find(strcmp({EEG.chanlocs(:).labels}, 'Photic'));                        
            
            disp('Flipping data');
            EEG.data = -EEG.data;
            disp('Rereferencing data except EKG and photic');
            EEG = pop_reref( EEG, [],'exclude',[ekgindex photicindex ]);  
            %notch
            disp('Filtering data');
            EEG = pop_eegfilt(EEG, 48, 52, 3300, 1, 1, 0);
            %passband
            EEG = pop_eegfilt(EEG, 1, 70, 6600, 0, 1, 0);

            %Downscale EKG            
            if ~isempty(ekgindex) 
                disp('Downscaling EKG');
                EEG.data(ekgindex,:) = EEG.data(ekgindex,:)/5;
            end            
                        
            %Order columns
            %montageOrder = {'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', ...
            %    'f11', 'F12', 'C3', 'C4', 'T7', 'T8', 'P3', 'P4', ...
            %    'P7', 'P8', 'P11', 'P12', 'O1', 'O2', 'A1', 'A2' ...
            %    'Fz', 'Cz', 'EKG', 'Photic'};
            %EEG = ScoreOrderChannels(EEG, montageOrder);
                                   
            eegplot( EEG.data, ...
                'winlength', 15,  ...
                'selectcommand', {1,1,1}, ...
                'srate', EEG.srate, ...
                'title', 'EEG plot', ...
                'command' , 'disp(''mouse event lead to command'');' , ...
                'events', EEG.event, ...
                'limits', [EEG.xmin EEG.xmax]*1000, ...
                'eloc_file', EEG.chanlocs, ...
                'scale', 'off', ...
                'selectcommand', {'ScoreMouseDown();', 'ScoreMouseMove();', 'ScoreMouseUp();'}, ...
                'ctrlselectcommand', {'disp(''ctrlmousedown'');', '', 'disp(''ctrlmouseup'');'} ...
                );                                
                                    
            clear ekgindex photicindex            
            assignin('base','EEG',EEG);
            eeg_checkset( EEG );                          
            
            evalin('base', 'eeglab redraw');
                              
            % Get EEGLAB plot
            existingPlot = findobj(0, 'tag', 'EEGPLOT');
            if isempty(existingPlot)
                warning('No EEG plot open in EEGLAB');
            elseif size(existingPlot,1)>1                        
                warn('More than one EEGPLOT open, using the first one');                        
                existingPlot = existingPlot(1);        
            end                           
            
            %Attach mouse handler
            set(existingPlot, 'WindowButtonDownFcn', {@MouseDown});            
        catch ME
            error(ME.message);
            assignin('base','STUDY', []); 
            assignin('base','CURRENTSTUDY', 0); 
            assignin('base','ALLEEG', []); 
            assignin('base','EEG', []); 
            assignin('base','CURRENTSET', []);  
        end
    end
end

