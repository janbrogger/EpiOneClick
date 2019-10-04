function varargout = EpiSnippetxcorr(templateIED, signal)

% IN
% episnippet1: IED template.
% episnippet2: signal to compare to IED template.
%
% OUT
% varargout: Vector "IEDdtw" of euclidian distances between warped templateIED at
% each samplepoint. IEDdtw > threshold = matches.

%Lenght of template
lengthIED = length(templateIED);
signallength = length(signal);
IEDxcorr = zeros(1,lengthIED);

%FOR TESTING PHASE. Maximum 2 minutes.
signallength = min(signallength, 60000);

%This is not the way to do it. Testing.
for i=1:(signallength-lengthIED)
    signalsnippet = signal(i:(i+lengthIED)); 
    if(mod(i,5000)==0) %debug
        i
    end
    [xc lag] = xcorr(templateIED,signalsnippet);
end

%[xc lag] = xcorr(templateIED,signal(1:signallength));
%IEDxcorr = [xc; lag];
IEDxcorr = max(xc);
varargout = {IEDxcorr};