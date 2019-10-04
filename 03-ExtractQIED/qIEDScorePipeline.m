% qIEDScorePipeline() - Purpose: To annotate IEDs automatically only from manually marked peaks, 
%   and then compare the resulting measures with results from manual
%   annotations.
% Usage:
%   >>  out = qIED( electrode, time );
%
% Inputs:
%   electrode       - first input of the function
%   time            - timesample clicked
%    
% Outputs:
%   Start, peak, end, slowend.
%
% See also: 
%   POP_qIED, EEGLAB 

% Copyright (C) 2019  Eivind Aanestad
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


% Find all peaks. Select most prominent and narrow as spike peak.
% x- and y-values for signal around peak.

%handle edge case of click very close to the beginning or end
siglengthBefore=500;
siglengthAfter=siglengthBefore;
edgeCase=0;
if (clicksample-siglengthBefore) < 1
    edgeCase=1;
    siglengthBefore=clicksample-1;
end
if (clicksample+siglengthAfter) > size(eegdata,2)
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

x = 1:1:length(signal);
t_spikepeak = ceil(length(signal)/2);
searchdistance = 25;
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
searchdistance = 100;
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
    elseif(ainterslope > 0.6)
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
lesserend = min(IEDspikepeak+100, length(signal)); 
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
    elseif( binterslope < -0.6)
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
minstartpos = IEDspikeend+83;
lesserend = min(IEDspikeend+400, length(signal));
hasSlow = minstartpos < length(signal);
%calculate slow-wave only if the signal is long enough for a slow wave

slowX = (IEDspikeend:1:lesserend);
slowsignal = [slowX; signal(slowX)];
slowXformin = IEDspikeend:1:lesserend;
slowYformin = smoothdata(signal(slowXformin),2,'movmean',83);
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
    gfoptions.Lower = [-2000 sAfterDischargeFirstidx 0];
    gfoptions.Upper = [2000 sAfterDischargeLastidx (sAfterDischargeLength)];
    sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;

    gfit = fit(slowsignal(1, :).', (slowsignal(2, :)+afterDischShift).','gauss1', gfoptions);
    YHat = gfit(slowsignal(1, :).')-afterDischShift;
    integralgfit = integrate(gfit, slowsignal(1, :).', 1);
    exactareagfit = integralgfit(length(integralgfit)) - integralgfit(1);
    exactareagfitsubtrapz = exactareagfit - trapz([1 sAfterDischargeLength], [(slowsignal(2, 1)+afterDischShift).', (slowsignal(2, sAfterDischargeLength)+afterDischShift)]);

    % Integration of afterdischarge using trapezoids. Area could be
    % overestemitated due to high-frequency activity. Smoothing? Filtering?
    integraltrapz = trapz((slowsignal(1, :)+afterDischShift).', (slowsignal(2, :)+afterDischShift)) ... 
                    - trapz([1 sAfterDischargeLength], [(slowsignal(2, 1)+afterDischShift).', (slowsignal(2, sAfterDischargeLength)+afterDischShift)]);
else
    sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;
    IEDslowend = IEDspikeend;
    IEDslowendy = IEDspikeendy;
end %hasslow
sduration(chan) = (IEDspikeend - IEDspikestart)*2;
%First halfwave (FH_)
onsetpeakexactAmplitude = signal(IEDspikepeak) - signal(IEDspikestart);
onsetpeakamplitude = round(onsetpeakexactAmplitude); onsetampl(chan) = onsetpeakamplitude;
descpeakexactAmplitude = signal(IEDspikepeak) - signal(IEDspikeend);
descpeakamplitude = round(descpeakexactAmplitude); deschampl(chan) = descpeakamplitude;
FH_duration = IEDspikepeak - IEDspikestart;    
FH_exactslope = onsetpeakamplitude/FH_duration; onsetslope(chan) = FH_exactslope/2;
FH_fullslope = round(FH_exactslope,1);   
SH_duration = IEDspikeend - IEDspikepeak;
SH_exactfullslope = descpeakamplitude/SH_duration; descslope(chan) = SH_exactfullslope/2;
SH_fullslope = round(SH_exactfullslope,1);
% Asymmetry, ascending slope to descending slope ratio
onset_desc_ratio = FH_fullslope/SH_fullslope;

%J.D.Frost 2nd derivate 8 ms around peak as suggested by PGL.

%Frost used 1 sample per 4 msec and calculated 2nd derivate aroun 5
%samplepoints. We use 500 samples/sec, 1sample/2msec. 
frostx = (IEDspikepeak-4):2:(IEDspikepeak+4);
frosty = signal(frostx);
frostsharpness(chan) = (abs( (frosty(5)-(2*frosty(3))+frosty(1))/2 ))/4;             

%Plot the IED template
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
f = figure('Name', 'Is this the IED you are looking for?');
subplot(2,1,1)
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
str_x = sprintf('%0.2f',frostsharpness(chan));
str_y = sprintf('%0.2f', onsetslope(chan));   
str_z = sprintf('%0.2f', descslope(chan));  
str_xx = sprintf('%0.2f', sduration(chan));
str_yy = sprintf('%0.2f', onsetampl(chan));
str_zz = sprintf('%0.2f', deschampl(chan));
legend([fig_signal fig_IED fig_slow_aprox fig_sp fig_so],'EEG signal', 'IED signal', 'Slow-wave approximation', 'Spike peak', 'Spike onset, spike end and slow-wave end', 'Location', 'southoutside')
pbaspect([xaspectr yaspectr 1]);

subplot(2,1,2)
resultline1 = strcat(['Sharpness: ', str_x,', Onsetslope: ', str_y, ', Descslope: ', str_z]);
resultline2 = strcat(['Duration: ', str_xx,', Onsetampl: ', str_yy, ', Descampl: ', str_zz]);
text(0.1, 0.8, resultline1);
hold on;
text(0.1, 0.6, resultline2);
if(hasSlow)
    str_ss = sprintf('%0.2f', exactareagfitsubtrapz/siglengthAfter);
    resultline3 = strcat(['Slowwave: ', str_ss]);
    text(0.1, 0.4, resultline3);
end
set(gca, 'xtick', [], 'ytick', []);
out = {IEDspikestart+clicksample-siglengthBefore, IEDspikepeak+clicksample-siglengthBefore, IEDspikeend+clicksample-siglengthBefore, IEDslowend+clicksample-siglengthBefore};


