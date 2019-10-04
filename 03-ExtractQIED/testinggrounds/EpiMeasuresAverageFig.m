function varargout = EpiMeasuresAverageFig(isIED, groupby)

% isIED     -   1: IEDs 2: sharp transients
% groupby   -   'age':Age categories on the form '<1' '1-9' etc.
%           -   'halford': Halford scale 1,2,3,4,5.
% varargout -   Matrix with average measures by age.
%
% Dependencies - 'PlotEpiAnnotationsUsersSnippets.mat' at location C:\Midlertidig_Lagring\Epileptiform\eegsnippets
%               This is an output file from PlotEpiAnnotations.m,
%               containing measures for each patient/EEG.

% Declare root folder containing EEG-snippets and check that it exists.
rootFolder = 'C:\Midlertidig_Lagring\Epileptiform\eegsnippets';
if ~isfolder(rootFolder)
    errorMessage = sprintf('Error, folder %s does not exist', rootFolder);
    uiwait(fprintf(errorMessage));
    return;
end

% Sample rate and conversion rate to milliseconds.
srate = 500;
sample2ms = 1000/srate;

% Declare and load input file
filePattern = fullfile(rootFolder, 'PlotEpiAnnotationsUsersSnippets.mat');
loadedstruct = load(filePattern);
ma = loadedstruct.UsersSnippets;
% Converting to units with milliseconds instead of samples.
measures = [ma.gAgeCat ma.durations*2 ma.onsetslopes/2 ma.onsetamplitudes ma.descslopes/2 ma.descamplitudes ma.afterDGaussAreas*2 ma.afterDDuration*2]
Nmeasures = length(measures)

% Declare save folder for plots.
saveFolder = rootFolder;

%varargin set to default if no input age categories are given.
if(strcmp(groupby, 'age'))
    gAgeCat = 1:1:10;
else
    halford = 1:1:5;
end

Nagecats = length(gAgeCat);
counters = ones(1,Nagecats);

% Calculate mean for each age category
for i = 1: Nagecats
    ind{i} = measures(:,1) == i; %indices of rows with correct agecat
    acmeasures{i} = measures(ind{i},:);
    acmeasures{i} = mean(acmeasures{i});
end

% plot (same plot window as in EpiSnippetsAverage)
xmin = -100*1000/srate;
xmax = 200*1000/srate;
ymin = -75;
ymax = 150;
xtime = xmin : xmax;
Nxtime = length(xtime);
Nsubplots = length(gAgeCat);
Nsubplotcollumns = 5;
Nsubplotrows = ceil(Nsubplots/Nsubplotcollumns);
xaspectr = (Nxtime/1000)*3;
yaspectr = (ymax - ymin)/100;
measurefig = figure();
set(measurefig, 'Name', 'AverageSpikeByAge', 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
for i = 1: Nsubplots
    subplot(Nsubplotrows,Nsubplotcollumns,i)
    spikestartx = 0;
    sstarty = 0;
    speakx = acmeasures{i}(4) / acmeasures{i}(3);
    speaky = acmeasures{i}(4);
    sendx = speakx + ( acmeasures{i}(6) / acmeasures{i}(5) );
    sendy = speaky - acmeasures{i}(6);    
    slowendx = sendx + acmeasures{i}(8);
    slowendy = 0;
    slowtime = sendx:1:slowendx;
    slowtimehalf = slowtime(ceil(end/2));
    slowplot = slowfig(slowtime, length(slowtime)*0.34, slowtimehalf, acmeasures{i}(7));
    slowtranslate = slowplot(1)-sendy; 
    plot(slowtime, slowplot-slowtranslate, 'LineWidth',2);
    hold on;
    line([xmin spikestartx speakx sendx], [0 sstarty speaky sendy], 'LineWidth',2);
    line([slowtime(end) xmax],[slowplot(end)-slowtranslate 0], 'LineWidth',2);   
    hold off;    
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    pbaspect([xaspectr yaspectr 1]);
    xticks([-200 0 200 400 600]);
    yticks([-50 0 50 100 150]);
    ax = gca;
    ax.FontSize = 6;
    xlabel('ms');
    ylabel('µV');      
    str = strcat({'Age '}, {gAgeCat2char(gAgeCat(i))});
    title(str);
end
figurestring = 'AverageMeasureFigByAge.png';
savePlotPattern = fullfile(saveFolder, 'User_25', '\epifig_byage\', figurestring);        
saveas(measurefig, savePlotPattern);
varargout = {acmeasures};

% Helping function to find gaussian with area like measured average
% Calling function with sigma "s" like 34% of slow
% after-wave length ensures a good fit.
function f = slowfig(x,s,b, area)

for i = 1:(area*5)
    f = i* ( gaussmf(x, [s b]) );
    y = trapz(x, f); 
    if(y >= area)
        break;
    end
end





