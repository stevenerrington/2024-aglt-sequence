% Initialize empty matrices for spike density function (SDF) data
% Each matrix is nNeurons x nTime and contains trial-averaged responses
% condA_SDF might represent a "congruent" condition, and condB_SDF an "incongruent" one
condA_SDF = []; 
condB_SDF = [];

% Clear any previously existing inputSDF variable and define new input SDF data
% The data is wrapped into a cell array for compatibility with the clustering function
clear inputSDF
inputSDF = {condA_SDF, condB_SDF};

% Define the time axis for each SDF input (in ms)
% This represents the full time window of the SDF, e.g., from -200 ms to +800 ms relative to some event
sdfTimes = {[-200:800], [-200:800]};

% Define the epoch of interest within the SDF time window
% This is typically the analysis window, such as 0 ms to 500 ms post-event
sdfEpoch = {[0:500], [0:500]};

% Define how each condition should be treated in clustering
% In this case, [1, 1] means both conditions are used with equal weight
colorMapping = [1, 1];

% Perform consensus clustering across neurons based on their activity during the defined epoch
% Arguments:
% - inputSDF: cell array of condition-wise SDFs
% - sdfTimes: time vector for each condition
% - '-e': epoch used for clustering
% - '-ei': input color mapping
% - '-er': epoch used for response extraction
% - '-c': consensus threshold (e.g., 0.5)
[sortIDs, idxDist, raw, respSumStruct, rawLink, myK] =...
    consensusCluster(inputSDF, sdfTimes, '-e', sdfEpoch, '-ei', colorMapping, '-er', sdfEpoch, '-c', 0.5);

% Close any open figures from the clustering process
close all

% Store the number of clusters discovered by the consensus clustering
nClusters_manual = myK;

% Initialize a cell array to store neuron indices for each cluster
clusterNeurons = [];

% Loop through each cluster and find which neurons belong to it
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:, nClusters_manual) == cluster_i);
end

% Plot the dendrogram to visualize the hierarchical clustering result
figuren('Renderer', 'painters', 'Position', [100 100 500 400]);
subplot(1,5,5)
[h, ~, outPerm] = dendrogram(rawLink, 0, 'Orientation', 'right');  % Create right-oriented dendrogram
set(gca, 'YDir', 'Reverse');  % Reverse Y-axis to match clustering order

% Color dendrogram branches based on cluster assignments
klDendroClustChange(h, rawLink, sortIDs(:, nClusters_manual))

% Tidy up dendrogram axes
set(gca, 'YTick', [], 'YLim', [1 length(rawLink)]); 
xlabel('Similarity')

% Plot the similarity matrix between neurons
subplot(1,5,[1:4]);

% Fill in the lower triangle of the symmetric matrix for visualization
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic, ir) = raw(ir, ic);
    end
end

% Display the similarity matrix with reordered neurons
imagesc(raw(outPerm, outPerm));  % Reorder by dendrogram permutation
colormap(flipud(cbrewer2('PuBu')));  % Use a perceptually uniform blue colormap
xlabel('Unit Number'); 
set(gca, 'YAxisLocation', 'Left');  % Align Y-axis to th
