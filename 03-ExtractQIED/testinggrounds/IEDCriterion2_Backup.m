    % Try using Morlet wavelets (#1) to estimate IED criterion 2 (different
    % wave-duration than the ongoing background activity, either shorter or
    % longer). Code developed from Mike X Cohen's online lecturelets. In
    % addition, xcorr (#2) and standard FFT "bandpass"(#3)
    
    % Wavelet cycles determines width of Gaussian. Few cycles: specific in
    % the time domain. Many cycles: specific in the freq domain and an 
    % effective band pass. 
    wavelet_cycles = 8;
    wavelet_freq = srate/length(signalSpike);
    wfreqround = round(wavelet_freq*2)/2; %Trying rounded frequency to neares 0.5 Hz.
    wavelet_time = -0.5:(1/srate):0.5;
    wavelet_half = (length(wavelet_time)-1) / 2;    
    normalizefactor = abs(max(signal)) + abs(min(signal));
    baselineBG = signal(1:t_spikestart - (1*srate));
    baselineBG_t = 1:1:length(baselineBG);
    ongoingBG = signal(t_spikestart - (1*srate) : (t_spikestart-1));
    ongoingBG_t = 1:1:length(ongoingBG);
    spikesignal = signal(t_spikestart : t_spikeend);
    spikesignal_t = 1:1:length(spikesignal);
    signal_t = 1:1:length(signal); 
    
    % FFT parameters
    nKern = length(wavelet_time);
    NongoingBG = length(ongoingBG);
    NbaselineBG = length(baselineBG);
    NconvongoingBG = nKern+NongoingBG-1;
    NconvbaselingeBG = nKern+NbaselineBG -1;
    Nsignal = length(signal);
    Nconvsignal = nKern+Nsignal-1;
    
     % X-axis for frequency domain plots and "arbitrary" max hz to get
    % reasonable power ratios.
    hz = linspace(0, srate/2, floor(NconvongoingBG/2)+1);
    [hzmax hzmaxidx] = find(hz==50);
    
    %Parameters to determine bandpass or Gauss SD limits.
    %Ghz = linspace(0, srate, NongoingBG);
    bpasshalf = round( ( 2^(log2(wfreqround)-3) ) *2)/2; 
    hpass = wfreqround-bpasshalf;
    lpass = wfreqround+bpasshalf;
    [~, pfidx] = min(abs(hz-wfreqround));
    [~, idxlp] = min(abs(hz-lpass));
    [~, idxhp] = min(abs(hz-hpass));
        
    % FFT 1 and 2 second window of background activity preceding IED
    ongoingBGX = (fft(ongoingBG, NconvongoingBG));
    %ongoingBGforGaussX = (fft(ongoingBG));
    baselineBGX = (fft(baselineBG, NconvbaselingeBG));
    
    %FFT of the whole epiSnippetPlus-signal. To help visual interpretation.
    signalX = (fft(signal, Nconvsignal));
    
    %Power law (1/frequency)    
    %Find N in 1/(F^N) from baseline. Ommit <1 Hz.
    y = log10((2*abs(baselineBGX(3:hzmaxidx))).^2);
    x = (hz(3:hzmaxidx));    
    y = y';
    x = x';
    X = [ones(length(x),1) x];
    regressor = X\y;
    yCalc1 = X*regressor;
    %linearfit = figure();
    %scatter(x,y)
    %hold on
    %plot(x,yCalc1)    
    %grid on
    a = regressor(2,1);
    b = regressor(1,1);       
    F = [ones(1,2) double(10.^(a*hz(3:hzmaxidx) + b))];    
    Rsquared = 1 - sum((y-yCalc1).^2) / sum((y - mean(y)).^2); 
    str_rsq1 = sprintf('%0.2f', Rsquared);
    str_rsq2 = strcat(['R-squared: ', str_rsq1]);
    %text(x(floor(length(x)/2)), double(yCalc1(floor(length(x)/2)))*1.3, str_rsq2);   
    %xlabel('frequency, 1-50 Hz.')
    %ylabel('log power')
    %title('Linear Regression Relation Between frequency & log(power)')
    %str = strcat(allFolders(l).name, 'Criterion2PowerLaw ', num2str(SearchResultID), '.png');
    %savePlotPattern = fullfile(saveFolder, allFolders(l).name, str);        
    %saveas(figure(linearfit), savePlotPattern);    
    %close(linearfit);
    
   
    % Morlet wavelet defined in the time domain.
    s           = wavelet_cycles / (2*pi*wfreqround); %((2*pi) * 8 * (log2(wfreqround))); 
    sine_wave   = exp(2*1i*pi*wfreqround.*wavelet_time);
    gauss_win   = exp( (-wavelet_time.^2) ./ (2*s^2) );
    cmw         = sine_wave .* gauss_win;
    cmwX        = fft(cmw,NconvongoingBG);
    cmwX        = cmwX./max(cmwX);  
    cmwXongoingBGXproduct   = cmwX.*ongoingBGX;
    %cmwXongoingBGX1Fproduct = cmwX(3:hzmaxidx).*ongoingBGX1F;
    cmwsignalX                   = fft(cmw,Nconvsignal);
    cmwsignalX                   = cmwsignalX./max(cmwsignalX);
    
    % Morlet wavelet defined by Gaussian in the frequency domain.
    % This to verify log2-relationship to bandwidth.
    s = hz(pfidx)-hz(idxhp);
    xGauss = hz-hz(pfidx);
    freqGaussX = exp(-.5*(xGauss/s).^2);
    freqGaussX = freqGaussX./max(freqGaussX);
    freqGauss = fftshift( ifft(freqGaussX, NconvongoingBG ) );
    freqGaussplot = ( freqGauss./max(freqGauss));
    freqGaussXhalf = floor( length(freqGaussX) / 2 );
    freqGaussplot = freqGaussplot(freqGaussXhalf:end-freqGaussXhalf);
    freqGaussXongoingBGforGaussX = freqGaussX.*ongoingBGX(1:length(freqGaussX));    
    
    similardurationpower    = (trapz(hz(3:hzmaxidx), (2*abs(cmwXongoingBGXproduct(3:hzmaxidx))/length(ongoingBG)).^2)); %AUC/power. *2 for negative FFT coeff, ^2 for power.
    freqGausspower          = (trapz(hz(3:hzmaxidx), (2*abs(freqGaussXongoingBGforGaussX(3:hzmaxidx))/length(ongoingBG)).^2));
    totalpower              = (trapz(hz(3:hzmaxidx), (2*abs(ongoingBGX(3:hzmaxidx))/length(ongoingBG)).^2)); %Total AUC/power for preceding background
    cmwpowerratio           = (similardurationpower/totalpower);
    freqGausspowerratio     = (freqGausspower/totalpower);
         
    % 2: FFT    
    %Gaussian sd from Morlet parameters. Asked MXC about it on his youtube.    
    fftcriterion2 = (trapz(hz(idxhp:idxlp),(2*abs(ongoingBGX(idxhp:idxlp))/length(ongoingBG)).^2));
    fftpowerratio = fftcriterion2/totalpower;
    
    % 3: XCORR
    r = xcorr(spikesignal, ongoingBG);
    ar = xcorr(spikesignal);
    meanr = mean(r);
    meanar = mean(ar);
    xcorratio = meanr/meanar;   