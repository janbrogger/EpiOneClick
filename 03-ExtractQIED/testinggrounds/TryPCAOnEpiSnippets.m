function TryPCAOnEpiSnippets(A)

% Testing primary component analysis on matrix A, where matrix A is varargout
% from function AverageEpiSnippets, i.e. measures from epileptiform events.
% Rows = patients. Columns = measures.


%Store relevant columns/measures in matrix Araw
Araw = table2array(A(:, 4:12));

%Hack to remove NaN column 6 (?)
%Araw(:,6) = [];

%Mean of rows of Araw
meanAraw = mean(Araw,1);

%Standard deviation of Araw
stdAraw = std(Araw,0,1);

%Subtract mean from Araw and divide by standard deviation
Astandardized = bsxfun(@minus, Araw, meanAraw);
Astandardized = bsxfun(@rdivide, Astandardized, stdAraw);

[coeff,score,latent] = pca(Astandardized);
str = sprintf("\t\t%5.2g\n",latent);
fprintf("Variances of principal components:\n");
fprintf(str);

biplot(coeff(:,1:3),'scores',score(:,1:3),'varlabels',{'Amplitude_ons','Amplitude_desc','Duration','SlopeOns', 'SlopeOns1','SlopeOns2','SlopeDesc','OnsDescRatio', 'AfterSlowGauss'});





