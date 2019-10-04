function varargout = EpiSnippetDTW(templateIED, signal, spikeend)

% IN
% episnippet1: IED template.
% episnippet2: signal to compare to IED template.
%
% OUT
% varargout: Vector "IEDdtw" of euclidian distances between warped templateIED at
% each samplepoint. IEDdtw > threshold = matches.
%
% TODO
% Consider splitting IED into spike and slow-wave and sum result. This to
% have an appropriately narrow sakoeChibaBand for the spike component.

%Lenght of template
lengthIED = length(templateIED);
spiketemplate = templateIED(1:spikeend);
slowtemplate = templateIED(spikeend+1:end);
signallength = length(signal);
IEDdtw = zeros(1,lengthIED);

%Specify width of SakoeChibaBand
DmatrixWidth = round(lengthIED/10);
DspikeWidth = round(length(spiketemplate)/10);
DslowWidth = round(length(slowtemplate)/10);

%FOR TESTING PHASE. Maximum 2 minutes.
signallength = min(signallength, 60000);
fprintf('0%');
prc100 = round(signallength/100);
prc10 = round(signallength/10);
%iterate through signal concatenating euclidian distance to output IEDdtw
for i=1:(signallength-lengthIED)
    %signalsnippet = signal(i:(i+lengthIED)); 
    spikesignal = signal(i:(i+spikeend));
    slowsignal = signal((i+spikeend):i+lengthIED);
    %Try to allign on 1st sample.
    %difference = signalsnippet(1)-templateIED(1);
    %signalsnippet = signalsnippet-difference;    
    if(mod(i,prc10)==0)       
        fprintf('%u', i/prc10*10);  
    elseif(mod(i,prc100)==0)         
        fprintf('.');
    end
    %Execute dynamic time warping to get distance between template and signal
    IEDdtw(i) = dtw(spiketemplate, spikesignal, DspikeWidth) + dtw(slowtemplate, slowsignal, DslowWidth);  
end


varargout = {IEDdtw}