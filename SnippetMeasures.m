function varargout = SnippetMeasures(userid, doPlot, doPlotfit, doTitle, saveMeasures)


% IMPORTANT:    Assert correct variable "rootfolder". 
%
% varargin:
% doPlot        - 1: Produce figures/plots. 0: No figures, only measures.
% doPlotfit     - 1: Plot annotations and fitted Gauss. 
%                 0: No annotations or Gauss
% doTitle       - 1: Plot titles with measures 0: No titles 
% saveMeasures  - 1: Save mesaures to excel 0: Don't save
% varargout:
%           - Quantitative measures.
%
% Dependencies  
% -Snippets at location epileptiform\1b-Matlab\data\User_?\snippets\


userstr = num2str(userid);

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = ['C:\Midlertidig_Lagring\EpiOneClick\Kural dataset snippets\User_' userstr '\'];
if ~isfolder(rootFolder)
    errorMessage = sprintf('Error, folder %s does not exist', rootFolder);
    uiwait(fprintf(errorMessage));
    return;
end

% Debugging missing gAgeCat
missingagecats = 0;

% Declare save folder for plots.
saveFolder = rootFolder;
FFTsaveFolder = [saveFolder '\FFTplots\'];
SnippetsignalsaveFolder = [saveFolder '\SnippetSignals\'];

%Gauss fit options
gfoptions = fitoptions('gauss1');
gfoptions.Upper = [Inf Inf 1000];
gfoptions.Lower = [-Inf -Inf 0];

% Subplot properties for basic measures
subplotrows = 10;
subplotcollumns = 8;
subplotn = subplotrows * subplotcollumns;
subplotindex = 1;
plotpagessindex = 1;

%struct for varargout
UsersSnippets.userid = [];
UsersSnippets.SearchResultEventIds = [];
UsersSnippets.StudyId = [];
UsersSnippets.gAgeCat = [];
UsersSnippets.onsetamplitudes = []; 
UsersSnippets.descamplitudes = [];
UsersSnippets.durations = [];
UsersSnippets.criterion2FFT = [];
UsersSnippets.precedingRMS = [];
UsersSnippets.onsetslopes = [];
UsersSnippets.onsethalfslopes1 = [];
UsersSnippets.onsethalfslopes2 = [];
UsersSnippets.descslopes = [];
UsersSnippets.onset_desc_ratio = [];
UsersSnippets.afterDGaussAreas = [];
UsersSnippets.afterDTrapzAreas = [];
UsersSnippets.afterDGaussTrapzAreas = [];
UsersSnippets.afterDDuration = [];
UsersSnippets.Rsquared1F = [];
UsersSnippets.frostsharpness = [];
UsersSnippets.annoStartx = [];UsersSnippets.annoStarty = [];
UsersSnippets.annoPeakx = [];UsersSnippets.annoPeaky = [];
UsersSnippets.annoEndx = [];UsersSnippets.annoEndy = [];
UsersSnippets.annoSlowEndx = [];UsersSnippets.annoSlowEndy = [];
UsersSnippets.halford = [];
UsersSnippets.isIED = [];
UsersSnippets.numberOfSelected = [];
userstring = ['User_' userstr];
% Abort if folder is not named "User_*"

filePattern = fullfile(rootFolder, '\snippets\', 'S*.mat');
theFiles = dir(filePattern);
% Exit if no snippet-files.
if(isempty(theFiles))
    fprintf('No snippets to process.');
    return;
end

%Loop through EEG-snippets
for k = 1 : length(theFiles)    
    try 
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(rootFolder, '\snippets\', baseFileName);
        %fprintf(1, 'Reading file %s \n', fullFileName);

        % Load snippet signal from one channel and measure
        % Spike amplitude, slope from spike onset to peak, slope from spike
        % onset to half maximum to peak.
        epiAnno = load(fullFileName);
        srate = epiAnno.snippet.srate;
        samplescale = (1/500)*srate; %This variable replaces originally project specific hardcoded sample-related calculations.
        t_spikestart = epiAnno.snippet.spikeStartX;
        t_spikepeak = epiAnno.snippet.spikePeak;
        t_spikeend = epiAnno.snippet.spikeEnd;
        t_afterslowend = epiAnno.snippet.afterDischargeEnd;
        afterdduration = t_afterslowend - t_spikeend;
        signal = double(epiAnno.snippet.EEGdata(epiAnno.snippet.clickedChannelIndex, :));
        signal_maxvalue = max(signal);
        xaxis = 1 : length(signal);
        x_signal = [xaxis; signal];
        signalHalfWave1 = x_signal(:, t_spikestart:t_spikepeak);
        signalHalfWave2 = x_signal(:, t_spikepeak:t_spikeend);
        signalSpike = [signalHalfWave1 signalHalfWave2(:,2:length(signalHalfWave2))];

        %Get annotated IED-points for scatter plot. Peak set as time 0.
        annoStartx = t_spikestart - t_spikepeak;
        annoStarty = signal(t_spikestart);
        annoPeakx = t_spikepeak - t_spikepeak; %0
        annoPeaky = signal(t_spikepeak);
        annoEndx = t_spikeend - t_spikepeak;
        annoEndy = signal(t_spikeend);
        annoSlowEndx = t_afterslowend - t_spikepeak;    
        annoSlowEndy = signal(t_afterslowend);   

        if (annoEndx - annoStartx < 4)
            error('Very short spike, shorter than 4 samples');
        end
        if (annoStartx >= annoPeakx)
            error('Wrong ordering of annotations, start is on or after peak');
        end
        if (annoPeakx >= annoEndx)
            error('Wrong ordering of annotations, peak is on or after end');        
        end

        duration = t_spikeend - t_spikestart;
        %First halfwave (FH_)
        onsetpeakexactAmplitude = signal(t_spikepeak) - signal(t_spikestart);
        onsetpeakamplitude = round(onsetpeakexactAmplitude);
        descpeakexactAmplitude = signal(t_spikepeak) - signal(t_spikeend);
        descpeakamplitude = round(descpeakexactAmplitude);
        halfamplitude = onsetpeakexactAmplitude/2 + signal(t_spikestart);
        temp = find(signalSpike(2,:) > halfamplitude);
        FH_FWHM = x_signal(:, (temp(1)+t_spikestart));   
        SH_FWHM = x_signal(:, (temp(length(temp))+t_spikestart));
        FH_duration = t_spikepeak - t_spikestart;    
        FH_exactslope = onsetpeakamplitude/FH_duration;
        FH_fullslope = round(FH_exactslope,1);     
        FH_exacthalfslope1 = (FH_FWHM(2,1) - signal(t_spikestart)) / (FH_FWHM(1,1) - t_spikestart);
        FH_exacthalfslope2 = (signal(t_spikepeak) - FH_FWHM(2,1)) / (t_spikepeak - FH_FWHM(1,1));
        FH_halfslope1 = round(FH_exacthalfslope1, 1);
        FH_halfslope2 = round(FH_exacthalfslope2, 1);
        %Second halfwave (SH_) (halfslope 1 and 2 are left to right)
        descpeakexactAmplitude = signal(t_spikepeak) - signal(t_spikeend);
        descpeakamplitude = round(descpeakexactAmplitude);
        SH_duration = t_spikeend - t_spikepeak;
        SH_exactfullslope = descpeakamplitude/SH_duration;
        SH_fullslope = round(SH_exactfullslope,1);     

        % Asymmetry, ascending slope to descending slope ratio
        onset_desc_ratio = FH_fullslope/SH_fullslope;

        %J.D.Frost 2nd derivate 8 ms around peak as suggested by PGL.    
        %Frost used 1 sample per 4 msec and calculated 2nd derivate aroun 5
        %samplepoints. We use 500 samples/sec, 1sample/2msec. 
        frostms = 4; %divide by 2 to get samples according to 500 sample/sec
        frostx = (t_spikepeak-round(4*samplescale)):round(2*samplescale):(t_spikepeak+round(4*samplescale);
        frosty = signal(frostx);
        frostsharpness = abs( (frosty(5)-(2*frosty(3))+frosty(1))/4 );

        % Gaussian1 fit to afterdischarge. Coefficient bounds to increase
        % performance. Gauss area is computed on a gauss curve shiftet
        % vertically by subtracting 'minimum' leading to overestimates in cases
        % of high voltage difference between start and stop.
        % Exaxtareagfitsubtrapz subtracts area under line between start and
        % stop voltage.
        signalAfterDischarge = x_signal(:, t_spikeend:epiAnno.snippet.afterDischargeEnd);
        %signalAfterDischarge(2,:) = smoothdata(signalAfterDischarge(2,:),2,'movmean',50);
        afterDischShift = -min(signalAfterDischarge(2,:));
        sAfterDischargeLength = length(signalAfterDischarge(2, :));
        sAfterDischargeLastidx = signalAfterDischarge(1, sAfterDischargeLength);
        sAfterDischargeFirstidx = signalAfterDischarge(1, 1);
        sAfterDischargeCenteridx = signalAfterDischarge(1, round(sAfterDischargeLength/2));
        gfoptions.Lower = [-2000 sAfterDischargeFirstidx 0];
        gfoptions.Upper = [2000 sAfterDischargeLastidx (sAfterDischargeLength)];
        afterDischargeChecksOut = (sAfterDischargeLength > 3 && sAfterDischargeFirstidx < sAfterDischargeLastidx);
        sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;


        % Calculate areas only if annotated aftercharge durations is > 3 
        % samples and ordered start/stop      
        exactareagfit = 0;
        integraltrapz = 0;
        exactareagfitsubtrapz = 0;
        if afterDischargeChecksOut
            gfit = fit(signalAfterDischarge(1, :).', (signalAfterDischarge(2, :)+afterDischShift).','gauss1', gfoptions);
            YHat = gfit(signalAfterDischarge(1, :).')-afterDischShift;
            integralgfit = integrate(gfit, signalAfterDischarge(1, :).', 1);
            exactareagfit = integralgfit(length(integralgfit)) - integralgfit(1);
            exactareagfitsubtrapz = exactareagfit - trapz([1 sAfterDischargeLength], [(signalAfterDischarge(2, 1)+afterDischShift).', (signalAfterDischarge(2, sAfterDischargeLength)+afterDischShift)]);
            sADareaforplot = round(exactareagfitsubtrapz);    

            % Integration of afterdischarge using trapezoids. Area could be
            % overestemitated due to high-frequency activity. Smoothing? Filtering?
            integraltrapz = trapz((signalAfterDischarge(1, :)+afterDischShift).', (signalAfterDischarge(2, :)+afterDischShift)) ... 
                            - trapz([1 sAfterDischargeLength], [(signalAfterDischarge(2, 1)+afterDischShift).', (signalAfterDischarge(2, sAfterDischargeLength)+afterDischShift)]);

        end %afterdischarge duration check

        %Try FFT-analysis of preceding epoch if signal has enough samples
        %(2*srate).     
        fftpowerratio = NaN;
        ongoingBG_RMS = NaN;    
        if(t_spikestart >= (2*srate+1))
            %Save 2*srate samples preceding IED
            ongoingBG = signal(t_spikestart-1 - (2*srate) : (t_spikestart-1));
            spikesignal = signal(t_spikestart : t_spikeend);
            Nspikesignal = length(spikesignal);   
            NongoingBG = length(ongoingBG);

            % X-axis for frequency domain plots and max hz at our lowpass
            %500 bins, should be spaced 0,5 Hz.
            hz = linspace(0, srate/2, floor(NongoingBG/2)+1); 
            %Find index in hz-vector for where element is 50Hz
            [~, hzmaxidx] = find(hz==50); 

            % Convert IED-duration to band in Hz. We use +- 10% of duration,
            % then divide sample rate by upper/lower duration bound.
            spikedurshort = duration * 0.9;
            spikedurlong = duration * 1.1;
            spikehzhigh = srate/spikedurshort;
            spikehzlow = srate/spikedurlong;
            %Locate the index in hz-vector. For integration and plotting.
            [~, lowidx] = min(abs(hz-spikehzlow));
            [~, highidx] = min(abs(hz-spikehzhigh));        

            % FFT 2 second window of background activity preceding IED.
            ongoingBGX = fft(ongoingBG);
            %Total power starting arbitrarily at index 5, or 2,5 Hz.
            totalpower    = (trapz(hz(5:hzmaxidx), (2*abs(ongoingBGX(5:hzmaxidx))/length(ongoingBG)).^2)); %Total AUC/power for preceding background  
            fftcriterion2 = (trapz(hz(lowidx:highidx),(2*abs(ongoingBGX(lowidx:highidx))/length(ongoingBG)).^2)); %Power within band corresponding to IED
            fftpowerratio = fftcriterion2/totalpower;

            % Root-mean-square of background activity preceding IED
            ongoingBG_RMS = rms(ongoingBG);
        end
        
        %BEMS
        %BEMS Age
    input_age = epiAnno.snippet.patientage;
    
    BEMS_age = 0;
    if(input_age < 10)
        BEMS_age = 16;
    elseif(input_age < 20)
        BEMS_age = 0;
    elseif(input_age < 60)
        BEMS_age = 12;
    elseif(input_age >= 60)
        BEMS_age = 25;
    end
    
    %BEMS Descending amplitude
    BEMS_descamp = 0;
    if(descpeakamplitude < 70)
        BEMS_descamp = 1;
    elseif(descpeakamplitude >= 70 && deschampl(chan) < 90)
        BEMS_descamp = 0;
    elseif(descpeakamplitude >= 90 && deschampl(chan) < 120)
        BEMS_descamp = 7;
    elseif(descpeakamplitude >= 120)
        BEMS_descamp = 17;
    end
    %BEMS onset slope
    BEMS_onsslope = 0;
    
    if(onsetslope(chan) < 1) %samples, not ms.
        BEMS_onsslope = 0;
    elseif(onsetslope(chan) >= 1 && onsetslope(chan) < 1.5)
        BEMS_onsslope = 4;
    elseif(onsetslope(chan) >= 1.5 && onsetslope(chan) < 2)
        BEMS_onsslope = 5;
    elseif(onsetslope(chan) >= 2)
        BEMS_onsslope = 11;
    end
    %BEMS spike to background power
    BEMS_spiketobg = 0;
    if(fftpowerratio*100 >= 8.6) %convert ratio to percent
        BEMS_spiketobg = 0;
    elseif(fftpowerratio*100 >= 4.7)
        BEMS_spiketobg = 9;
    elseif(fftpowerratio*100 >= 2.6)
        BEMS_spiketobg = 6;
    elseif(fftpowerratio*100 < 2.6)
        BEMS_spiketobg = 14;
    end
    %BEMS slow after-wave area
    BEMS_slow = 0;
    if(exactareagfitsubtrapz/srate < 5)
        BEMS_slow = 0;
    elseif(exactareagfitsubtrapz/srate < 10)
        BEMS_slow = 6;
    elseif(exactareagfitsubtrapz/srate < 20)
        BEMS_slow = 11;
    elseif(exactareagfitsubtrapz/srate >= 20)
        BEMS_slow = 19;
    end   
    %BEMS total
    BEMS_score = BEMS_age + BEMS_descamp + BEMS_onsslope + BEMS_spiketobg + BEMS_slow;

        % Plot original signal, annotations, gauss model.
        if(doPlot == 1)
            if(subplotindex == 1)
                %Basic measures: basicfig. 
                figuretitle = strcat('User_', userstr, ' Basic measures', ' plot number ', num2str(plotpagessindex));   
                basicfig = figure('Name', figuretitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]); 
                %Waveletfigure: waveletfig.
                waveletfigtitle = strcat('Wavelet measures for SRID ', num2str(SearchResultID));
                FFTfig = figure('Name', waveletfigtitle, 'NumberTitle', 'off','Position', [0 0 800 600]);
                set(FFTfig,'Visible', 'off');            
                figure(basicfig);
            else
                figure(basicfig);
            end

            currentsubplotn = mod(k, subplotn);
            cutsignal = signal(t_spikestart - srate : t_spikeend + srate);
            Lagcutsignal = t_spikestart - srate;
            cutx = t_spikestart-srate:t_spikeend + srate;
            cutsstart = t_spikestart - Lagcutsignal;
            cutpeak = t_spikepeak - Lagcutsignal;
            cutend = t_spikeend - Lagcutsignal;
            cutslowend = epiAnno.snippet.afterDischargeEnd - Lagcutsignal;
            if currentsubplotn == 0
                currentsubplotn = subplotn;
            end        
            subplot(subplotrows, subplotcollumns, currentsubplotn)

            if(doPlotfit == 0)
                plot(cutx,cutsignal); 
                if(doTitle == 1)
                    str = sprintf('sreID %s \t', num2str(SearchResultID));
                    title(str);
                end
            else
                plot(cutx,cutsignal, '-p', 'MarkerIndices', [cutsstart cutpeak cutend cutslowend], 'MarkerFaceColor', 'red', 'MarkerSize', 5);
                hold;
                line([FH_FWHM(1,1) SH_FWHM(1,1)], [halfamplitude halfamplitude], 'Color', 'r', 'LineWidth', 2);    
                if afterDischargeChecksOut
                    plot((t_spikeend:epiAnno.snippet.afterDischargeEnd), YHat); 
                end
                if(doTitle == 1)
                    str = sprintf('sreID %s \t', num2str(SearchResultID));
                    title(str)
                end

                figure(FFTfig);
                waveletfigtitle = strcat('FFT measures for SRID ', num2str(SearchResultID));
                set(FFTfig, 'Name', waveletfigtitle, 'NumberTitle', 'off', 'Position', [50 50 1000 300]);
                set(FFTfig,'Visible', 'off'); 
                if(length(signal) >= 8*srate)
                    tempsecrange = 4;
                elseif(length(signal) >= 4*srate)
                    tempsecrange = 2;
                else
                    tempsecrange = 0;
                end
                %Some snippets might be too short
                if(tempsecrange >= 2)
                    subplot(2,2,[1,2])            
                    plot(xaxis(1:t_spikestart-(tempsecrange*0.5*srate))/srate-tempsecrange, signal(1:t_spikestart-(tempsecrange*0.5*srate)),'g' ... 
                        , xaxis(t_spikestart-(tempsecrange*0.5*srate):t_spikestart)/srate-tempsecrange, signal(t_spikestart-(tempsecrange*0.5*srate):t_spikestart), 'b' ...
                        , xaxis(t_spikestart:t_afterslowend)/srate-tempsecrange, signal(t_spikestart:t_afterslowend), 'r' ... 
                        , xaxis(t_afterslowend:end)/srate-tempsecrange, signal(t_afterslowend:end), 'k', 'linew', 1);
                    title('Green: Baseline. Blue: Preceding background activity. Red: IED');
                    ymax = min(ceil(signal_maxvalue / 100)*100, 300);
                    set(gca, 'ylim', [-ymax ymax]);

                    subplot(2,2,3)
                    BGXplot = ((2*abs(ongoingBGX(5:hzmaxidx))).^2)/length(ongoingBG);
                    plotfactor = max(BGXplot);
                    plot(hz(5:hzmaxidx), BGXplot, 'b');
                    line([hz(lowidx) hz(lowidx)],[0 plotfactor], 'color', 'r');
                    line([hz(highidx) hz(highidx)],[0 plotfactor], 'color', 'r');
                    set(gca, 'xlim', [0 50]);
                    title('Preceding background activity FFT');

                    subplot(2,2,4);
                    title('Results');
                    %str_u = sprintf('%0.2f', Rsquared);
                    str_v = sprintf('%0.1f', fftpowerratio*100);
                    str_v2 = sprintf('%0.1f', hz(lowidx));
                    str_v3 = sprintf('%0.1f', hz(highidx));
                    resultline0 = strcat('How much power in preceding activity corresponds to spike duration?');            
                    resultline2 = strcat(['FFT,high pass ', str_v2, ', low pass ', str_v3, ', percent of total power:', str_v]);
                    text(0, 1.5, resultline0);
                    text(0,0.5, resultline2);
                    set(gca, 'ylim', [0 2]);
                    set(gca, 'Visible', 'off'); 
                else
                    text(0, 0.5, 'Snippet too short for background analysis');
                end

                %Save plot to disk for post-analysis visual verification
                str = strcat(userstring, 'Criterion2FFT_SREID_', num2str(SearchResultID), '.png');
                saveFFTPlotPattern = fullfile(FFTsaveFolder, str);        
                saveas(figure(FFTfig), saveFFTPlotPattern);
                clf(FFTfig);
            end
            hold off;        

            % Save figure and initialize a new figure for every N subplots.     
            if (currentsubplotn == subplotn || k == length(theFiles))
                str = strcat(userstring, '_plot ', num2str(plotpagessindex), '.png');
                savePlotPattern = fullfile(SnippetsignalsaveFolder, str);        
                saveas(figure(basicfig), savePlotPattern);
                close(basicfig);
                plotpagessindex = plotpagessindex + 1;
               if k ~= length(theFiles)
                figuretitle = strcat(userstring, ' Basic measures', ' plot number ', num2str(plotpagessindex));   
                basicfig = figure('Name', figuretitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
               end
            end
            subplotindex = subplotindex + 1;
        end %doPlot true/false

        % Concatenate snippet measures to UsersSnippets
        %debugging = epiAnno.snippet.StudyId
        if ~isempty(theFiles)
            UsersSnippets.userid = [UsersSnippets.userid; userid];
            UsersSnippets.onsetamplitudes = [UsersSnippets.onsetamplitudes; onsetpeakexactAmplitude];
            UsersSnippets.descamplitudes = [UsersSnippets.descamplitudes; descpeakexactAmplitude];
            UsersSnippets.durations = [UsersSnippets.durations; duration];
            %UsersSnippets.criterion2CMW = [UsersSnippets.criterion2CMW; cmwpowerratio];
            UsersSnippets.criterion2FFT = [UsersSnippets.criterion2FFT; fftpowerratio];
            %UsersSnippets.criterion2XCOR = [UsersSnippets.criterion2XCOR; xcorratio];
            UsersSnippets.onsetslopes = [UsersSnippets.onsetslopes; FH_fullslope];
            UsersSnippets.onsethalfslopes1 = [UsersSnippets.onsethalfslopes1; FH_halfslope1];
            UsersSnippets.onsethalfslopes2 = [UsersSnippets.onsethalfslopes2; FH_halfslope2];
            UsersSnippets.descslopes = [UsersSnippets.descslopes; SH_fullslope];
            UsersSnippets.onset_desc_ratio = [UsersSnippets.onset_desc_ratio; onset_desc_ratio];
            %UsersSnippets.afterDGaussIntegrals = [UsersSnippets.afterDGaussIntegrals single(integralgfit)];
            UsersSnippets.afterDGaussAreas = [UsersSnippets.afterDGaussAreas; exactareagfit];
            UsersSnippets.afterDTrapzAreas = [UsersSnippets.afterDTrapzAreas; integraltrapz];
            UsersSnippets.afterDGaussTrapzAreas = [UsersSnippets.afterDGaussTrapzAreas; exactareagfitsubtrapz];
            UsersSnippets.afterDDuration = [UsersSnippets.afterDDuration; afterdduration];
            UsersSnippets.SearchResultEventIds = [UsersSnippets.SearchResultEventIds; SearchResultID];
            UsersSnippets.StudyId = [UsersSnippets.StudyId; epiAnno.snippet.StudyId];        
            %UsersSnippets.Rsquared1F = [UsersSnippets.Rsquared1F; Rsquared];
            UsersSnippets.frostsharpness = [UsersSnippets.frostsharpness; frostsharpness];
            UsersSnippets.annoStartx = [UsersSnippets.annoStartx; annoStartx];
            UsersSnippets.annoStarty = [UsersSnippets.annoStarty; annoStarty];
            UsersSnippets.annoPeakx = [UsersSnippets.annoPeakx; annoPeakx];
            UsersSnippets.annoPeaky = [UsersSnippets.annoPeaky; annoPeaky];
            UsersSnippets.annoEndx = [UsersSnippets.annoEndx; annoEndx];
            UsersSnippets.annoEndy = [UsersSnippets.annoEndy; annoEndy];
            UsersSnippets.annoSlowEndx = [UsersSnippets.annoSlowEndx; annoSlowEndx];
            UsersSnippets.annoSlowEndy = [UsersSnippets.annoSlowEndy; annoSlowEndy];
            UsersSnippets.precedingRMS = [UsersSnippets.precedingRMS; ongoingBG_RMS];
            UsersSnippets.isIED = [UsersSnippets.isIED; epiAnno.snippet.isIED];
            UsersSnippets.numberOfSelected = [UsersSnippets.numberOfSelected; numberOfSelected];
            if(~isnan(halford))
                UsersSnippets.halford = [UsersSnippets.halford; halford];
            else
                UsersSnippets.halford = [UsersSnippets.halford; NaN];
            end
            if isfield(epiAnno.snippet, 'gAgeCat')
                agecategory = gAgeCat2num(epiAnno.snippet.gAgeCat);
                UsersSnippets.gAgeCat = [UsersSnippets.gAgeCat; agecategory];            
            else
                UsersSnippets.gAgeCat = [UsersSnippets.gAgeCat; NaN];
                missingagecats = missingagecats + 1;
            end    
        end
    catch ME
        warning(ME.message);
    end
end %EEG-snippets


% Save UsersSnippets struct to save rootfolder
filenamepart = strcat(userstring, '_Paper2_EpiMeasures');
str = strcat(filenamepart, '.mat');
saveStructPattern = fullfile(saveFolder, str);
save(saveStructPattern, 'UsersSnippets');

%Convert struct UsersSnippets to table and save to rootfolder
outputtable(:,1) = UsersSnippets.userid;
outputtable(:,2) = UsersSnippets.SearchResultEventIds;
outputtable(:,3) = UsersSnippets.StudyId;
outputtable(:,4) = UsersSnippets.gAgeCat;
outputtable(:,5) = UsersSnippets.onsetamplitudes;
outputtable(:,6) = UsersSnippets.descamplitudes;
outputtable(:,7) = UsersSnippets.durations;
%outputtable(:,8) = UsersSnippets.criterion2CMW;
outputtable(:,8) = UsersSnippets.criterion2FFT;
%outputtable(:,10) = UsersSnippets.criterion2XCOR;
outputtable(:,9) = UsersSnippets.onsetslopes;
outputtable(:,10) = UsersSnippets.onsethalfslopes1;
outputtable(:,11) = UsersSnippets.onsethalfslopes2;
outputtable(:,12) = UsersSnippets.descslopes;
outputtable(:,13) = UsersSnippets.onset_desc_ratio;
outputtable(:,14) = UsersSnippets.afterDGaussAreas;
outputtable(:,15) = UsersSnippets.afterDTrapzAreas;
outputtable(:,16) = UsersSnippets.afterDGaussTrapzAreas;
outputtable(:,17) = UsersSnippets.afterDDuration;
%outputtable(:,18) = UsersSnippets.Rsquared1F;
outputtable(:,18) = UsersSnippets.frostsharpness;
outputtable(:,19) = UsersSnippets.annoStartx;
outputtable(:,20) = UsersSnippets.annoStarty;
outputtable(:,21) = UsersSnippets.annoPeakx;
outputtable(:,22) = UsersSnippets.annoPeaky;
outputtable(:,23) = UsersSnippets.annoEndx;
outputtable(:,24) = UsersSnippets.annoEndy;
outputtable(:,25) = UsersSnippets.annoSlowEndx;
outputtable(:,26) = UsersSnippets.annoSlowEndy;
outputtable(:,27) = UsersSnippets.halford;
outputtable(:,28) = UsersSnippets.precedingRMS;
outputtable(:,29) = UsersSnippets.numberOfSelected;
%outputtable(:,29) = UsersSnippets.isIED;
T = array2table(outputtable, 'VariableNames', {'UserId', 'SearchResultEventId', 'StudyId', 'gAgeCat', 'Amplitude_onset', 'Amplitude_desc', 'Duration', 'Crit2FFT', 'Onsetslope', 'OnsetHalfSlope1', 'OnsetHalfSlope2', 'Descslope', 'OnsetDescendingRatio', 'AfterSlowGArea', 'AfterSlowTArea', 'AfterSlowGTArea', 'AfterSlowDuration', 'FrostSharpness', 'annoStartx','annoStarty','annoPeakx', 'annoPeaky', 'annoEndx', 'annoEndy','annoSlowEndx','annoSlowEndy', 'halford', 'precedingRMS', 'numberOfSelected'});
if(saveMeasures)
    filenamepart = strcat(userstring, '_Paper2_EpiMeasures', '.xlsx');
    filename_excel = fullfile(saveFolder, filenamepart);
    if exist(filename_excel, 'file')
      delete(filename_excel);
    end
    writetable(T, filename_excel, 'FileType', 'spreadsheet');
end
varargout = {T};



