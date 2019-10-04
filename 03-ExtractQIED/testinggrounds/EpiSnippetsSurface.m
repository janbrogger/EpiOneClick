function varargout = EpiSnippetsSurface(varargin)

% varargin -    Age categories on the form '<1' '1-9' etc.
% varargout -   Matrix with all signals.
%
% Dependencies - ScorePipeline episnippets at location C:\Midlertidig_Lagring\Epileptiform\eegsnippets
%              - Requires EpiSnippets to be of same length

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

% Loop through subfolders (users)
for l = 1 : length(allFolders)
filePattern = fullfile(rootFolder, allFolders(l).name, 'Epi_SearchResultEventId*.mat');
theFiles = dir(filePattern);

%struct for varargout
if(nargin < 1)
    gAgeCat = {'<1', '1-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-101'};
else
    for i = 1 : nargin
        gAgeCat{i} = varargin{i};
    end
end

%Properties for episnippet-plots. Time relative to peak
srate = 500;
xmin = -100;
xmax = 300;
xtime = xmin : xmax;
Nxtime = length(xtime);

%Loop through EEG-snippets
counters = [1 1 1 1 1 1 1 1 1 1];
for k = 1 : length(theFiles)    
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(rootFolder, allFolders(l).name, baseFileName);
    
    % Load episnippet signal from one channel and measure
    % Spike amplitude, slope from spike onset to peak, slope from spike
    % onset to half maximum to peak.
    epiAnno = load(fullFileName);
    if isfield(epiAnno.epiSnippet, 'gAgeCat')
        trimmedsnippetagecat = strtrim(epiAnno.epiSnippet.gAgeCat); 
        for m = 1 : length(gAgeCat)
            if(strcmp(trimmedsnippetagecat, gAgeCat{m})) %|| strcmp('all', gAgeCat))
                peaktime = epiAnno.epiSnippet.spikePeak;
                peakchan = epiAnno.epiSnippet.clickedChannelIndex;
                subsig = epiAnno.epiSnippet.EEGdata(peakchan, (peaktime+xmin):(peaktime+xmax));
                outmatrix{m}(counters(m),:) = subsig;
                counters(m) = counters(m) + 1;
            end %If wanted AgeCat  
        end %loop through agegroups
    end %if agecat exists
end %EEG-snippets
end %subfolders/users


outmatrixMEAN = zeros(xmax-xmin+1, length(gAgeCat));
for i = 1 : length(gAgeCat)
    tempmatrix = squeeze(outmatrix{i}(:,:));
    outmatrixMEAN(:,i) = mean(tempmatrix);    
end

%3D-plot, average spike by age categories            
avgplot3d = figure();
set(avgplot3d, 'Name', 'AverageSpikeByAge3D', 'NumberTitle', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]); %distorted by remote desktop?
bigplotstr = 'AverageSpikeByAge3d.png';
z3d = length(gAgeCat):-1:1;
[X,Y] = meshgrid(xtime*1000/srate, z3d);
surfmatrix = outmatrixMEAN.';
s = surf(X, Y, surfmatrix);
az = -6;
el = 18;
view(az, el);
grid off;
axis off;
ylabel('Age group');
xlabel('Time in millseconds from peak.');
zlabel('µV');
for i = 1:90
    az = az + 4;
    view(az,el);
    pause(0.2);
    str_fnum = sprintf('%02d', i);
    str_frame = strcat('avgspikeanim_', str_fnum, '.png');
    savePlotPattern = fullfile(saveFolder, allFolders(l).name, '\epifig_byage\', str_frame); 
    saveas(avgplot3d, savePlotPattern);
end
az = -6;
el = 18;
grid on;
axis on;
colorbar;
view(az, el);       
savePlotPattern = fullfile(saveFolder, allFolders(l).name, '\epifig_byage\', bigplotstr); 
saveas(avgplot3d, savePlotPattern);

varargout = {outmatrix};




