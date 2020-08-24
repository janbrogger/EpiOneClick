% qIEDScorePipeline()  
%
% Purpose: To annotate IEDs automatically only from % manually marked peaks.
% Usage:
%   >>  out = qIED( electrode, time );
%
% Inputs:
%   electrode       - clicked channel/electrode
%   time            - clicked timesample
%    
% Outputs:
%   Start, peak, end, slowend.
%
% See also: 
%   POP_qIED, EEGLAB 

% Copyright (C) 2019  Eivind Aanestad and Jan Brogger
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function out = qIEDScorePipeline(electrode,sample)
    settings = EpiOneClickSettings();
    EEG = evalin('base','EEG');
    srate = EEG.srate;
    samplescale = 500/srate; %This variable replaces originally project specific hardcoded sample-related calculations.
    
    %Output folder
    outputfolder = pwd;
    if(isfield(settings, 'SaveFolder'))
        outputfolder = settings.SaveFolder;
    end
    
    % Declare qIED parameters
    sharpness = [];
    onsetslope = [];
    onsetampl = [];
    descslope = [];
    descampl = [];
    sduration = [];
    swarea = [];

    %tmprej = evalin('base','TMPREJ');
    eegdata = evalin('base', 'EEG.data');
    chan = electrode;
    %clicktime = str2double(time.String); %time in seconds
    clicksample = sample;
    %startIED = round(tmprej(1))
    %endIED = round(tmprej(2))


    %% Find all peaks. Select most prominent and narrow as spike peak.
    % x- and y-values for signal around peak.

    %Handle edge case of click very close to the beginning or end
   
    siglengthBefore=srate*2+1+round(100*samplescale);
    siglengthAfter=siglengthBefore;
    edgeCase=0; 
    if (clicksample-(siglengthBefore)) < 1 %Require 200 ms + 2 seconds before peak (for background calculations)
        edgeCase=1;
        siglengthBefore=clicksample-1;
    end
    if (clicksample+(siglengthBefore)) > size(eegdata,2)
        edgeCase=1;
        siglengthAfter=size(eegdata,2)-clicksample;
    end
    if edgeCase==1
        disp(['clicksample: ' num2str(clicksample)]);
        disp(['siglengthBefore: ' num2str(siglengthBefore) ' siglengthAfter: ' num2str(siglengthAfter)]);
        disp(['clicksample-siglengthBefore: ' num2str(clicksample-siglengthBefore) ' clicksample+siglengthAfter: ' num2str(clicksample+siglengthAfter)]);
        error('This is an edge case at the beginning or end of EEG. It does not work. Skip it for now.');
    end
    signal = eegdata(chan,(clicksample-siglengthBefore):(clicksample+siglengthAfter));
    EEGdata = eegdata(:,(clicksample-siglengthBefore):(clicksample+siglengthAfter));
    x = 1:1:length(signal);
    t_spikepeak = ceil(length(signal)/2);
    searchdistance = round(25/samplescale);
    xaroundpeak = t_spikepeak-searchdistance:1:t_spikepeak+searchdistance;
    yaroundpeak = signal(xaroundpeak);
    %locate peak indexes, then convert x-values back to global index
    [~, peaksidx, ~, ~] = findpeaks(yaroundpeak);
    if(all(peaksidx<1))
        return; %did not find a peak
    end
    peakpos = peaksidx + t_spikepeak - searchdistance;
    [~, minindx] = min(abs(peakpos-t_spikepeak));
    %Now we have the peak closest to manual mark.
    IEDspikepeak = peakpos(minindx);
    IEDspikepeaky = signal(IEDspikepeak);

    %Find local minima before peak. Assume spike start within 200ms or 100 samples from peak. 
    %First, go with the local minimum closest to peak. If the next local minimum
    %has a lower voltage and delta x is lower than 0.625 delta y, then
    %replace spike start with this minimum. 0.625 relates to aspect ratio
    %in visual analysis (1cm per ?V and 3cm per second).
    searchdistance = round(100/samplescale);
    xaroundstart = t_spikepeak-searchdistance:1:t_spikepeak;
    yaroundstart = signal(xaroundstart);
    [minimalogic, minimaprominence] = islocalmin(yaroundstart);
    if(all(minimalogic<1))
        return; %did not find spike start
    end
    xstartminima = xaroundstart(minimalogic); 
    pminima = minimaprominence(minimalogic);
    aiterationlength = length(xstartminima);
    tempminimum = xstartminima(aiterationlength);
    tempminimumy = signal(tempminimum);
    for i = aiterationlength-1:-1:1
        %Replace currentminimum if next local minimum has
        %double the slope and greater matlab-prominence.
        atempslope = (IEDspikepeaky - tempminimumy) / (IEDspikepeak - tempminimum);
        acurrentslope = (IEDspikepeaky - signal(xstartminima(i))) / (IEDspikepeak - xstartminima(i));
        ainterslope = (tempminimumy - signal(xstartminima(i))) / (tempminimum - xstartminima(i));
        acurrentprominence = pminima(i);
        acurrentVolt = signal(xstartminima(i));
        tempVolt = signal(tempminimum);
        if(acurrentVolt > tempVolt)
            break;        
        elseif(ainterslope > 0.6*srate/500)
            tempminimum = xstartminima(i);
            tempminimumy = signal(tempminimum);
        end        
    end

    %Now we have the spike start.
    IEDspikestart = tempminimum;
    IEDspikestarty = signal(IEDspikestart);

    %Find local minima. Assume spike end within 200ms or 100 samples from peak. 
    %First, go with the local minimum furthest away from peak. If there is 
    %another local minimum closer to peak for which the slope doubles and higher
    %matlabprominence, or if voltage is lower, choose this minimum instead.
    lesserend = min(IEDspikepeak+(100/samplescale), length(signal)); 
    xaroundend = IEDspikepeak:1:lesserend;
    yaroundend = signal(xaroundend);
    [minimalogic, minimaprominence] = islocalmin(yaroundend);
    if(all(minimalogic<1))
        return; %did not find spike end
    end
    xendminima = xaroundend(minimalogic); 
    pminima = minimaprominence(minimalogic);
    iterationlength = length(xendminima);
    tempminimum = xendminima(1);
    tempminimumy = signal(tempminimum);
    for i = 2:1:iterationlength
        %Replace currentminimum if next local minimum has
        %double the slope and greater matlab-prominence.
        btempslope = (tempminimumy-IEDspikepeaky) / (tempminimum-t_spikepeak);
        bcurrentslope = ((signal(xendminima(i))-IEDspikepeaky)) / (xendminima(i)-t_spikepeak);
        binterslope = (signal(xendminima(i))- tempminimumy) / (xendminima(i)-tempminimum);
        bcurrentprominence = pminima(i);
        bcurrentVolt = signal(xendminima(i));
        tempVolt = signal(tempminimum);
        if(bcurrentVolt > tempVolt)
            break;
        elseif( binterslope < -0.6*srate/500)
            tempminimum = xendminima(i);
            tempminimumy = signal(tempminimum);
        end    
    end

    %Now we have the spike end.
    IEDspikeend = tempminimum;
    IEDspikeendy = signal(IEDspikeend);

    % Slow after-wave. Assuming duration less than 800ms or 400 samples, 
    % and duration > 166 ms (<6Hz or >83 samples).
    % When sub-optimal duration, Gauss will still approximate
    % area well, e.g. in the case of JME fast SW.
    % Gaussian1 fit to afterdischarge. Coefficient bounds to increase
    % performance. Gauss area is computed on a gauss curve shiftet
    % vertically by subtracting 'minimum' leading to overestimates in cases
    % of high voltage difference between start and stop.
    % Exaxtareagfitsubtrapz subtracts area under line between start and
    % stop voltage.
    minstartpos = IEDspikeend+83/samplescale;
    lesserend = min(IEDspikeend+400/samplescale, length(signal));
    hasSlow = minstartpos < length(signal);
    %calculate slow-wave only if the signal is long enough for a slow wave

    slowX = (IEDspikeend:1:lesserend);
    slowsignal = [slowX; signal(slowX)];
    slowXformin = IEDspikeend:1:lesserend;
    slowYformin = smoothdata(signal(slowXformin),2,'movmean',(83/samplescale));
    [slowminima, slowminimaprominences] = islocalmin(slowYformin); 
    %inverseslowY = max(slowYformin)-slowYformin; %try findpeaks method
    %[Minima, MinIdx] = findpeaks(inverseslowY); %try findpeaks method
    slowminimalprominencesY = slowminimaprominences(slowminima);
    slowminimalprominencesX = slowXformin(slowminima);
    %slowminimalprominencesX = slowXformin(MinIdx); %try findpeaks method
    hasSlow = ~isempty(slowminimalprominencesX);
    if(~hasSlow)
        [~,idx] = min(slowYformin);
        slowminimalprominencesX = slowXformin(idx);
        slowminimalprominencesY = slowYformin(idx);
        hasSlow = 1;
    end    
    %try to optimize identification of slow-wave end. 
    %trying harder smoothing and choosing first minimum occuring after
    %spikeend+offset. Replace with next minimum if lower voltage within reason.
    %Experimenting with prominence. Replace if difference to next considered
    %prominence is greater than a fourth of the total distribution width of
    %prominences.
    %slowpromlaggedX = slowminimalprominencesX(slowminimalprominencesX > minstartpos);
    %slowpromlaggedY = slowminimalprominencesY(slowminimalprominencesX > minstartpos);
    slowXminima = slowXformin(slowminima);
    slowYminima = slowYformin(slowminima);
    slowXminimalagged = slowXminima(slowXminima > minstartpos);
    slowYminimalagged = slowYminima(slowXminima > minstartpos);
    if(~isempty(slowXminimalagged))
        %[~, idx] = max(slowpromlaggedY);
        currentmin = slowYminimalagged(1);
        IEDslowend = slowXminimalagged(1);
        minspread = max(slowYminimalagged) - min(slowYminimalagged);
        for i=2:length(slowYminimalagged)
            mindiff = currentmin - slowYminimalagged(i) ;
            if(slowYminimalagged(i)<currentmin && mindiff>0.25*minspread)
                IEDslowend = slowXminimalagged(i);
                currentmin = slowYminimalagged(i);
            end
        end
    else
        IEDslowend = IEDspikeend; %slowminimalprominencesX(1); 
    end
    IEDslowendy = signal(IEDslowend); 
    
    %% Calculate qIED
    hasSlow = (minstartpos < length(signal) && (IEDslowend-IEDspikeend>3));
    if(hasSlow)           
        %adjust slow-wave signal length so it ends at prominent minimum.
        slowX = (IEDspikeend:1:IEDslowend);
        slowsignal = [slowX; signal(slowX)];
        gfoptions = fitoptions('gauss1');     
        afterDischShift = -min(slowsignal(2,:));
        sAfterDischargeLength = length(slowsignal(2, :));
        sAfterDischargeLastidx = slowsignal(1,sAfterDischargeLength);
        sAfterDischargeFirstidx = slowsignal(1,1);
        gfoptions.Lower = [(-2000/samplescale) sAfterDischargeFirstidx 0];
        gfoptions.Upper = [(2000/samplescale) sAfterDischargeLastidx (sAfterDischargeLength)];
        sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;

        gfit = fit(slowsignal(1, :).', (slowsignal(2, :)+afterDischShift).','gauss1', gfoptions);
        YHat = gfit(slowsignal(1, :).')-afterDischShift;
        integralgfit = integrate(gfit, slowsignal(1, :).', 1);
        exactareagfit = integralgfit(length(integralgfit)) - integralgfit(1);
        exactareagfitsubtrapz = exactareagfit - trapz([1 sAfterDischargeLength], [(slowsignal(2, 1)+afterDischShift).', (slowsignal(2, sAfterDischargeLength)+afterDischShift)]);
        slowwaveweber = exactareagfitsubtrapz/srate; %Convert area from µV*samples to µV*seconds
        
        % Integration of afterdischarge using trapezoids. Area could be
        % overestemitated due to high-frequency activity. Smoothing? Filtering?
        integraltrapz = trapz((slowsignal(1, :)+afterDischShift).', (slowsignal(2, :)+afterDischShift)) ... 
                        - trapz([1 sAfterDischargeLength], [(slowsignal(2, 1)+afterDischShift).', (slowsignal(2, sAfterDischargeLength)+afterDischShift)]);
    else
        sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;
        IEDslowend = IEDspikeend;
        IEDslowendy = IEDspikeendy;
        slowwaveweber = 0;
    end %hasslow
    sduration(chan) = (IEDspikeend - IEDspikestart);
    %First halfwave (FH_)
    onsetpeakexactAmplitude = signal(IEDspikepeak) - signal(IEDspikestart);
    onsetpeakamplitude = round(onsetpeakexactAmplitude); onsetampl(chan) = onsetpeakamplitude;
    descpeakexactAmplitude = signal(IEDspikepeak) - signal(IEDspikeend);
    descpeakamplitude = round(descpeakexactAmplitude); deschampl(chan) = descpeakamplitude;
    FH_duration = IEDspikepeak - IEDspikestart;    
    FH_exactslope = onsetpeakamplitude/FH_duration; onsetslope(chan) = FH_exactslope*(srate/1000);
    FH_fullslope = round(FH_exactslope,1);   
    SH_duration = IEDspikeend - IEDspikepeak;
    SH_exactfullslope = descpeakamplitude/SH_duration; descslope(chan) = SH_exactfullslope*(srate/1000);
    SH_fullslope = round(SH_exactfullslope,1);
    % Asymmetry, ascending slope to descending slope ratio
    onset_desc_ratio = FH_fullslope/SH_fullslope;

    %J.D.Frost 2nd derivate 8 ms around peak as suggested by PGL.

    %Frost used 1 sample per 4 msec and calculated 2nd derivate aroun 5
    %samplepoints. We use 500 samples/sec, 1sample/2msec. 
    frostx = (IEDspikepeak-round(4/samplescale)):round(2/samplescale):(IEDspikepeak+round(4/samplescale));
    frosty = signal(frostx);
    frostsharpness(chan) = (abs( (frosty(5)-(2*frosty(3))+frosty(1))/2 ))/4;   
    
    %Try FFT-analysis of preceding epoch if signal has enough samples
    %(2*srate).     
    fftpowerratio = NaN;
    ongoingBG_RMS = NaN;    
    if(edgeCase==0)
        %Save 2*srate samples preceding IED
        ongoingBG = signal(IEDspikestart-1 - (2*srate) : (IEDspikestart-1));
        spikesignal = signal(IEDspikestart : IEDspikeend);
        Nspikesignal = length(spikesignal);   
        NongoingBG = length(ongoingBG);

        % X-axis for frequency domain plots and max hz at our lowpass
        %500 bins, should be spaced 0,5 Hz.
        hz = linspace(0, srate/2, floor(NongoingBG/2)+1); 
        %Find index in hz-vector for where element is 50Hz
        [~, hzmaxidx] = find(hz==50); 

        % Convert IED-duration to band in Hz. We use +- 10% of duration,
        % then divide sample rate by upper/lower duration bound.
        spikedurshort = sduration(chan) * 0.9;
        spikedurlong = sduration(chan) * 1.1;
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
        spiketobackground = fftpowerratio*100;
        % Root-mean-square of background activity preceding IED
        ongoingBG_RMS = rms(ongoingBG);
    end   
    sdurationMS = sduration(chan) *(1000/srate); %milliseconds.
    %% Result strings
    str_x = sprintf('%0.2f',frostsharpness(chan));
    str_y = sprintf('%0.2f', onsetslope(chan));
    str_z = sprintf('%0.2f', descslope(chan)); 
    str_xx = sprintf('%0.2f', sdurationMS); 
    str_yy = sprintf('%0.2f', onsetampl(chan));
    str_zz = sprintf('%0.2f', deschampl(chan));    
    resultsharp = strcat(['Sharpness: ',str_x,]); 
    resultonslope = strcat(['Onsetslope: ',str_y]);
    resultdeslope = strcat(['Descslope: ',str_z]);
    resultdur = strcat(['Duration: ',str_xx]); 
    resultonsamp = strcat(['Onsetampl: ', str_yy]);
    resultdescamp = strcat(['Descampl: ', str_zz]);
    resultslow = '';
    if(hasSlow)
        str_ss = sprintf('%0.2f', exactareagfitsubtrapz/srate);
        resultslow = strcat(['Slowwave: ', str_ss]);
    end
    resultagestr = '';
    resultnumchanstr = '';
    resultfocusidstr = '';
    resultregionsstr = '';
    input_age = 0;
    input_focusid = 0;
    input_regionstr = "";
    input_numchan = 0;
    BEMS_age = 0;
    BEMS_descamp = 0;
    BEMS_onsslope = 0;
    BEMS_spiketobg = 0;
    BEMS_slow = 0;
    BEMS_score = 0;
    
    %% Plot the IED template and results
    Nxtime = length(signal);
    ymin = min(signal);
    ymax = max(signal);
    if(hasSlow)
        IEDsignal = (IEDspikestart:IEDslowend);
    else
        IEDsignal = (IEDspikestart:IEDspikeend);
    end
    xaspectr = (Nxtime*2/1000)*3;
    yaspectr = (ymax - ymin)/100;
    fheight=700;
    f = figure('Name', 'Is this the IED you are looking for?','Position',[360,300,450,fheight]);    
    resultinfo1 = uicontrol('Style','text', 'String',resultsharp,'Position',[50,fheight-200,400,20]);
    resultinfo2 = uicontrol('Style','text', 'String',resultonslope,'Position',[50,fheight-225,400,20]);
    resultinfo3 = uicontrol('Style','text', 'String',resultdeslope,'Position',[50,fheight-250,400,20]);
    resultinfo4 = uicontrol('Style','text', 'String',resultdur,'Position',[50,fheight-275,400,20]);
    resultinfo5 = uicontrol('Style','text', 'String',resultonsamp,'Position',[50,fheight-300,400,20]);
    resultinfo6 = uicontrol('Style','text', 'String',resultdescamp,'Position',[50,fheight-325,400,20]);
    resultinfo7 = uicontrol('Style','text', 'String',resultslow,'Position',[50,fheight-350,400,20]);    
    resultinfo8a = uicontrol('Style','text', 'String','Enter age (integer):','Position',[50,fheight-425,200,20]);
    resultinputage = uicontrol('Style','edit', 'String','Skip age','Position',[50,fheight-450,200,20],'Callback',@agebutton_Callback);
    resultinfo9a = uicontrol('Style','text', 'String','Enter number of channels (integer):','Position',[50,fheight-475,240,20]);
    resultinputnumchan = uicontrol('Style','edit', 'String',resultnumchanstr,'Position',[50,fheight-500,240,20],'Callback',@numchanbutton_Callback);
    resultinfo10a = uicontrol('Style','text', 'String','Enter focus ID (integer):','Position',[50,fheight-525,240,20]);
    resultinputfocusid = uicontrol('Style','edit', 'String',resultfocusidstr,'Position',[50,fheight-550,240,20],'Callback',@focusidbutton_Callback);
    resultinfoBEMS = uicontrol('Style','text', 'String','BEMS: ','Position',[50,fheight-575,400,20]);
    savebutton = uicontrol('Style','pushbutton', 'String','Save and close','Position',[50,fheight-600,200,20],'Callback',@savebutton_Callback);
    ha = axes('Units','pixels','Position',[5,fheight-220,440,200]);
    align([savebutton,resultinfo1,resultinfo2,resultinfo3,resultinfo4,resultinfo5,resultinfo6,resultinfo7,resultinfo8a,resultinputage,resultinfo9a,resultinputnumchan, resultinfo10a, resultinputfocusid, resultinfoBEMS],'Center','None');
    f.Units = 'normalized';
    ha.Units = 'normalized';
    savebutton.Units = 'normalized';
    resultinfo1.Units = 'normalized';
    resultinfo2.Units = 'normalized';
    resultinfo3.Units = 'normalized';
    resultinfo4.Units = 'normalized';
    resultinfo5.Units = 'normalized';
    resultinfo6.Units = 'normalized';
    resultinfo7.Units = 'normalized';
    resultinfo8a.Units = 'normalized';
    resultinputage.Units = 'normalized';
    resultinfo9a.Units = 'normalized';
    resultinputnumchan.Units = 'normalized';
    resultinfo10a.Units = 'normalized';
    resultinputfocusid.Units = 'normalized';
    resultinfoBEMS.Units = 'normalized';
    fig_signal = plot(x,signal);
    title(num2str(chan));
    hold on;
    fig_IED = plot(IEDsignal, signal(IEDsignal),'LineWidth', 1, 'color', 'k');
    fig_sp = scatter(IEDspikepeak, IEDspikepeaky,7, 'red', 'DisplayName', 'Spike peak');
    %fig_minima_o = scatter(xstartminima, signal(xstartminima),5,'green', 'DisplayName', 'Spike onset')
    %fig_minima_e = scatter(xendminima, signal(xendminima),5,'green', 'DisplayName', 'Spike end')
    if(hasSlow)
        %plot(slowXformin, smoothdata(slowYformin,2,'movmean',50));
        %scatter(slowminimalprominencesX, signal(slowminimalprominencesX),5,'green')
        fig_slow_e = scatter(IEDslowend, IEDslowendy,9, 'blue');
        fig_slow_aprox = plot(slowsignal(1, :), YHat, 'red');
    end
    fig_so = scatter(IEDspikestart, IEDspikestarty,9, 'blue');
    fig_se = scatter(IEDspikeend, IEDspikeendy,9,'blue');
    legend([fig_signal fig_IED fig_slow_aprox fig_sp fig_so],'EEG signal', 'IED signal', 'Slow-wave approximation', 'Spike peak', 'Spike onset, spike end and slow-wave end', 'Location', 'southoutside')
    pbaspect([xaspectr yaspectr 1]);        
    set(gca, 'xtick', [], 'ytick', []);
    out = {IEDspikestart+clicksample-siglengthBefore, IEDspikepeak+clicksample-siglengthBefore, IEDspikeend+clicksample-siglengthBefore, IEDslowend+clicksample-siglengthBefore};


%% Supporting functions. (Save,..)
    
    function savebutton_Callback(source,eventdata) 
    % Save snippet and data to mat-file and close window.          
        snippet.srate = EEG.srate;
        snippet.chanLocs = EEG.chanlocs;
        %snippet.clickedChannel = EEG.chanloc(
        snippet.clickedSample = sample;
        snippet.clickedChannelIndex = chan;
        snippet.signal = signal;
        snippet.EEGdata = EEGdata;
        snippet.spikeStartX = IEDspikestart;
        snippet.spikeStartY = IEDspikestarty;
        snippet.spikePeakX = IEDspikepeak;
        snippet.spikePeakY = IEDspikepeaky;
        snippet.spikeEndX = IEDspikeend;
        snippet.spikeEndY = IEDspikeendy;    
        snippet.Amplitude_onset = onsetpeakamplitude;
        snippet.Amplitude_desc = descpeakamplitude;
        snippet.FrostSharpness = frostsharpness(chan);
        snippet.Slope_Onset =  onsetslope(chan);
        snippet.Slope_Desc = descslope(chan);
        snippet.SpikeToBackground = spiketobackground;
        snippet.PrecedingRMS = ongoingBG_RMS;
        snippet.SlowWaveWeber = slowwaveweber;
        snippet.hasSlow = hasSlow;
        snippet.slowEndX = IEDslowend;
        snippet.sduration = sduration(chan);
                      
        resultagestr = resultinputage.String;        
        input_age = str2double(resultagestr);
        if isnan(input_age) || fix(input_age) ~= input_age
          resultinputage.String = "Invalid";
        else
            snippet.patientage = input_age;
            snippet.BEMS_age = BEMS_age;
            snippet.BEMS_descamp = BEMS_descamp;
            snippet.BEMS_onsslope = BEMS_onsslope;
            snippet.BEMS_spiketobg = BEMS_spiketobg;
            snippet.BEMS_slow = BEMS_slow;  
        end        
        
        %Only proceed to save when data have been entered correctly
        numchanvisual  = str2double(resultinputnumchan.String);
        if(isnan(numchanvisual))
            resultnumchanstr.String = "Invalid";
            return;
        end        
        snippet.numchanvisual = numchanvisual;
        
        focusIDvisual = str2double(resultinputfocusid.String);
        if(isnan(focusIDvisual))
            resultfocusidstr.String = "Invalid";
            return;
        end
        snippet.focusIDvisual = focusIDvisual;     

        snippet.EEGfilepath = EEG.setname;
        
        %Use clickedsample/srate as postfix (store one unique IED-candidate per
        %second.
        [filepath,name,ext] = fileparts(EEG.setname);  
        snippetpostfix = num2str(round(clicksample/srate));
        fileName = ['snippet_' name '_' snippetpostfix ];
        snippet.BEMS_filename = fileName;
        filePath = fullfile(outputfolder, fileName);            
        save(filePath, 'snippet');
        close(f)
    end
    
    function agebutton_Callback() 
        % Once age is entered BEMS can be calculated.
        %BEMS Age
        resultagestr = resultinputage.String;        
        input_age = str2double(resultagestr);
        if isnan(input_age) || fix(input_age) ~= input_age
          resultinputage.String = "Invalid";
          return;
        end
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
        if(deschampl(chan) < 70)
            BEMS_descamp = 1;
        elseif(deschampl(chan) >= 70 && deschampl(chan) < 90)
            BEMS_descamp = 0;
        elseif(deschampl(chan) >= 90 && deschampl(chan) < 120)
            BEMS_descamp = 7;
        elseif(deschampl(chan) >= 120)
            BEMS_descamp = 17;
        end
        %BEMS onset slope
        BEMS_onsslope = 0;

        if(onsetslope(chan) < 1) 
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
        if(spiketobackground >= 8.6)
            BEMS_spiketobg = 0;
        elseif(spiketobackground >= 4.7)
            BEMS_spiketobg = 9;
        elseif(spiketobackground >= 2.6)
            BEMS_spiketobg = 6;
        elseif(spiketobackground < 2.6)
            BEMS_spiketobg = 14;
        end
        %BEMS slow after-wave area
        BEMS_slow = 0;
        if(slowwaveweber < 5)
            BEMS_slow = 0;
        elseif(slowwaveweber < 10)
            BEMS_slow = 6;
        elseif(slowwaveweber < 20)
            BEMS_slow = 11;
        elseif(slowwaveweber >= 20)
            BEMS_slow = 19;
        end   
        %BEMS total
        BEMS_score = BEMS_age + BEMS_descamp + BEMS_onsslope + BEMS_spiketobg + BEMS_slow;
        str_bems = sprintf('%0.0f', BEMS_score);
        resultinfoBEMS.String = strcat(['BEMS: ', str_bems]);
    end

    function numchanbutton_Callback(src,~)     
        input_numchan = str2double(resultinputnumchan.String);
        if isnan(input_numchan) || fix(input_numchan) ~= input_numchan
          resultinputnumchan.String = "Invalid";
          resultinfoBEMS.String = "BEMS: Invalid";
          return;
        end
    end

    function focusidbutton_Callback(src,~)    
        input_focusid = str2double(resultinputfocusid.String);
        if isnan(input_focusid) || fix(input_focusid) ~= input_focusid
          resultinputfocusid.String = "Invalid";
          return;
        end
    end

    function regionsbutton_Callback(src,~)    
%       Skipping this check for now.
%         input_regionstr = resultinputregions.String;
%         if ~(contains(input_regionstr, 'F') || contains(input_regionstr, 'C') || contains(input_regionstr, 'T') || contains(input_regionstr, 'P') || contains(input_regionstr, 'O')) ...
%                 || ~(contains(input_regionstr, 'left') || contains(input_regionstr, 'right') || contains(input_regionstr, 'mid'))
%           resultinputregions.String = "Invalid";
%           return;
%         end
    end
end 