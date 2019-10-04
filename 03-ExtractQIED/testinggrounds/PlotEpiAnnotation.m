function varargout = PlotEpiAnnotation(doPlot, SearchResultEventId)

% varargin -    doPlot : 1: Produce figures/plots. 2: No figures, only measures.
%               SearchResultEventId: Specific event to be plotted/measured.
% varargout -   Quantitative measures per User.
%
% Dependencies - ScorePipeline episnippets at location C:\Midlertidig_Lagring\Epileptiform\eegsnippets
%
% Saves figures into Epileptiform\eegsnippets\User*\
% Saves varargout into Epileptiform\eegsnippets\User*EpiPlot*
%
% PlotEpiAnnotations de

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = 'C:\Midlertidig_Lagring\Epileptiform\eegsnippets';
if ~isdir(rootFolder)
    errorMessage = sprintf('Error, folder %s does not exist', rootFolder);
    uiwait(fprintf(errorMessage));
    return;
end

% Declare save folder for plots.
saveFolder = rootFolder;

% Subfolders/Users containing EEG-snippets
folderPattern = fullfile(rootFolder, 'User*');
allFolders = dir(folderPattern);

%Gauss fit options
gfoptions = fitoptions('gauss1');
gfoptions.Upper = [Inf Inf 1000];
gfoptions.Lower = [-Inf -Inf 0];

% Subplot properties
subplotrows = 5;
subplotcollumns = 4;
subplotn = subplotrows * subplotcollumns;
plotsindex = 1;

%struct for varargout
UsersSnippets.userid = single([]);
UsersSnippets.SearchResultEventIds = [];
UsersSnippets.amplitudes = []; 
UsersSnippets.onsetslopes = [];
UsersSnippets.onsethalfslopes1 = [];
UsersSnippets.onsethalfslopes2 = [];
UsersSnippets.descslopes = [];
UsersSnippets.onset_desc_ratio = [];
%UsersSnippets.afterDGaussCoefficients = [];
%UsersSnippets.afterDGaussIntegrals = single([]);
UsersSnippets.afterDGaussAreas = [];
UsersSnippets.afterDTrapzAreas = [];
UsersSnippets.afterDGaussTrapzAreas = [];

% Loop through subfolders (users)
for l = 1 : length(allFolders)
filePattern = fullfile(rootFolder, allFolders(l).name, ['Epi_SearchResultEventId_' num2str(SearchResultEventId) '.mat' ]);
theFiles = dir(filePattern);

% Initialize UsersSnippets and figure if *mat files exist.
if ~isempty(theFiles)
    figuretitle = strcat(allFolders(l).name, '_single_plot_', num2str(SearchResultEventId), num2str(plotsindex));
    if(doPlot == 1)
        figure('Name', figuretitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
        fig = gcf;
    end  
end

%Loop through EEG-snippets
for k = 1 : length(theFiles)    
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
    fprintf(1, 'Reading file %s \n', fullFileName);
    
    % Load episnippet signal from one channel and measure
    % Spike amplitude, slope from spike onset to peak, slope from spike
    % onset to half maximum to peak.
    epiAnno = load(fullFileName);
    signal = epiAnno.epiSnippet.EEGdata(epiAnno.epiSnippet.clickedChannelIndex, :);
    x = 1 : length(signal);
    x_signal = [x; signal];
    signalHalfWave1 = x_signal(:, epiAnno.epiSnippet.spikeStart:epiAnno.epiSnippet.spikePeak);
    signalHalfWave2 = x_signal(:, epiAnno.epiSnippet.spikePeak:epiAnno.epiSnippet.spikeEnd);
    signalSpike = [signalHalfWave1 signalHalfWave2(:,2:length(signalHalfWave2))];
   
    
    duration = epiAnno.epiSnippet.spikeEnd - epiAnno.epiSnippet.spikeStart;
    %First halfwave (FH_)
    onsetpeakexactAmplitude = signal(epiAnno.epiSnippet.spikePeak) - signal(epiAnno.epiSnippet.spikeStart);
    onsetpeakamplitude = round(onsetpeakexactAmplitude);
    halfamplitude = onsetpeakexactAmplitude/2 + signal(epiAnno.epiSnippet.spikeStart);
    temp = find(signalSpike(2,:) > halfamplitude);
    FH_FWHM = x_signal(:, (temp(1)+epiAnno.epiSnippet.spikeStart));   
    SH_FWHM = x_signal(:, (temp(length(temp))+epiAnno.epiSnippet.spikeStart));
    FH_duration = epiAnno.epiSnippet.spikePeak - epiAnno.epiSnippet.spikeStart;    
    FH_exactslope = onsetpeakamplitude/FH_duration;
    FH_fullslope = round(FH_exactslope,1);     
    FH_exacthalfslope1 = (FH_FWHM(2,1) - signal(epiAnno.epiSnippet.spikeStart)) / (FH_FWHM(1,1) - epiAnno.epiSnippet.spikeStart);
    FH_exacthalfslope2 = (signal(epiAnno.epiSnippet.spikePeak) - FH_FWHM(2,1)) / (epiAnno.epiSnippet.spikePeak - FH_FWHM(1,1));
    FH_halfslope1 = round(FH_exacthalfslope1, 1);
    FH_halfslope2 = round(FH_exacthalfslope2, 1);
    %Second halfwave (SH_) (halfslope 1 and 2 are left to right)
    descpeakexactAmplitude = signal(epiAnno.epiSnippet.spikePeak) - signal(epiAnno.epiSnippet.spikeEnd);
    descpeakamplitude = round(descpeakexactAmplitude);
    SH_duration = epiAnno.epiSnippet.spikeEnd - epiAnno.epiSnippet.spikePeak;
    SH_exactfullslope = descpeakamplitude/SH_duration;
    SH_fullslope = round(SH_exactfullslope,1);     
    
    % Asymmetry, ascending slope to descending slope ratio
    onset_desc_ratio = FH_fullslope/SH_fullslope;
    
    % Gaussian1 fit to afterdischarge. Coefficient bounds to increase
    % performance. Gauss area is computed on a gauss curve shiftet
    % vertically by subtracting 'minimum' leading to overestimates in cases
    % of high voltage difference between start and stop.
    % Exaxtareagfitsubtrapz subtracts area under line between start and
    % stop voltage.
    signalAfterDischarge = x_signal(:, epiAnno.epiSnippet.spikeEnd:epiAnno.epiSnippet.afterDischargeEnd);
    afterDischShift = -min(signalAfterDischarge(2,:));
    sAfterDischargeLength = length(signalAfterDischarge(2, :));
    sAfterDischargeLastidx = signalAfterDischarge(1, sAfterDischargeLength);
    sAfterDischargeFirstidx = signalAfterDischarge(1, 1);
    sAfterDischargeCenteridx = signalAfterDischarge(1, round(sAfterDischargeLength/2));
    gfoptions.Lower = [-2000 sAfterDischargeFirstidx 0];
    gfoptions.Upper = [2000 sAfterDischargeLastidx (sAfterDischargeLength)];
    afterDischargeChecksOut = (sAfterDischargeLength > 3 && sAfterDischargeFirstidx < sAfterDischargeLastidx);
    sADareaforplot = 0; integralgfit = 0; exactareagfit = 0; exactareagfitsubtrapz = 0; integraltrapz = 0;
    
    % Calculate areas only if aftercharge durations is > 3 samples and
    % ordered start/stop      
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
        
    % Plot original signal, annotations, gauss model.
    if(doPlot == 1)
        currentsubplotn = mod(k, subplotn);
        if currentsubplotn == 0
            currentsubplotn = subplotn;
        end
        subplot(subplotrows, subplotcollumns, currentsubplotn)
        plot(x,signal, '-p', 'MarkerIndices', [epiAnno.epiSnippet.spikeStart epiAnno.epiSnippet.spikePeak epiAnno.epiSnippet.spikeEnd, epiAnno.epiSnippet.afterDischargeEnd], 'MarkerFaceColor', 'red', 'MarkerSize', 5)    
        hold;
        line([FH_FWHM(1,1) SH_FWHM(1,1)], [halfamplitude halfamplitude], 'Color', 'r', 'LineWidth', 2);    
        if afterDischargeChecksOut
        plot(signalAfterDischarge(1, :), YHat); 
        end
        hold off;
        %plot(gfit, epiAnno.epiSnippet.spikeStart:epiAnno.epiSnippet.spikeEnd, signalSpike); 
        str = sprintf('srID %s \t A: %s \t S: %s \t SlowArea: %s', num2str(epiAnno.epiSnippet.SearchResultEventId), num2str(onsetpeakamplitude), num2str(FH_fullslope), num2str(sADareaforplot));
        title(str)

        % Save figure and initialize a new figure for every N subplots.     
        if (currentsubplotn == subplotn || k == length(theFiles))
            str = strcat(allFolders(l).name, 'plot ', num2str(plotsindex), '.png');
            savePlotPattern = fullfile(saveFolder, allFolders(l).name, str);        
            saveas(fig, savePlotPattern);
          % uiwait(gcf);
           plotsindex = plotsindex + 1;
           if k ~= length(theFiles)
            figuretitle = strcat(allFolders(l).name, ' plot number ', num2str(plotsindex));
            figure('Name', figuretitle, 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
            fig = gcf;
           end
        end
    end %doPlot true/false
    
    % Concatenate snippet measures to UsersSnippets
    useridasint = str2num(erase(allFolders(l).name, "User_"))
    if ~isempty(theFiles)
        UsersSnippets.userid = [UsersSnippets.userid; single(useridasint)];
        UsersSnippets.amplitudes = [UsersSnippets.amplitudes; onsetpeakamplitude];
        UsersSnippets.onsetslopes = [UsersSnippets.onsetslopes; FH_fullslope];
        UsersSnippets.onsethalfslopes1 = [UsersSnippets.onsethalfslopes1; FH_halfslope1];
        UsersSnippets.onsethalfslopes2 = [UsersSnippets.onsethalfslopes2; FH_halfslope2];
        UsersSnippets.descslopes = [UsersSnippets.descslopes SH_fullslope];
        UsersSnippets.onset_desc_ratio = [UsersSnippets.onset_desc_ratio; onset_desc_ratio];
        %UsersSnippets.afterDGaussIntegrals = [UsersSnippets.afterDGaussIntegrals single(integralgfit)];
        UsersSnippets.afterDGaussAreas = [UsersSnippets.afterDGaussAreas; exactareagfit];
        UsersSnippets.afterDTrapzAreas = [UsersSnippets.afterDTrapzAreas; integraltrapz];
        UsersSnippets.afterDGaussTrapzAreas = [UsersSnippets.afterDGaussTrapzAreas; exactareagfitsubtrapz];
        UsersSnippets.SearchResultEventIds = [UsersSnippets.SearchResultEventIds; single(epiAnno.epiSnippet.SearchResultEventId)];
    end
end %EEG-snippets

end %subfolders/users
% Save UsersSnippets struct to save rootfolder
str = 'PlotEpiAnnotationUsersSnippet.mat';
saveStructPattern = fullfile(saveFolder, str);
save(saveStructPattern, 'UsersSnippet');

%Convert struct UsersSnippets to table and save to rootfolder
outputtable(:,1) = UsersSnippets.userid;
outputtable(:,2) = UsersSnippets.SearchResultEventIds;
outputtable(:,3) = UsersSnippets.amplitudes;
outputtable(:,4) = UsersSnippets.onsetslopes;
outputtable(:,5) = UsersSnippets.onsethalfslopes1;
outputtable(:,6) = UsersSnippets.onsethalfslopes2;
outputtable(:,7) = UsersSnippets.descslopes;
outputtable(:,8) = UsersSnippets.onset_desc_ratio;
outputtable(:,9) = UsersSnippets.afterDGaussAreas;
outputtable(:,10) = UsersSnippets.afterDTrapzAreas;
outputtable(:,11) = UsersSnippets.afterDGaussTrapzAreas;
T = array2table(outputtable, 'VariableNames', {'UserId', 'SearchResultEventId', 'Amplitude', 'Onsetslope', 'OnsetHalfSlope1', 'OnsetHalfSlope2', 'Descslope', 'OnsetDescendingRatio', 'AfterSlowGArea', 'AfterSlowTArea', 'AfterSlowGTArea'});
filename_excel = fullfile(saveFolder, 'UserSnippet.xlsx');
if exist(filename_excel, 'file');
  delete(filename_excel);
end
writetable(T, filename_excel, 'FileType', 'spreadsheet');

varargout = {T};



