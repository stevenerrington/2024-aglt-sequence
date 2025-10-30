function mds_neuron_pipeline(sdf, varargin)
% MDS_NEURON_PIPELINE   Explore common motifs across neurons using MDS.
%   mds_neuron_pipeline(sdf) expects sdf: [nNeurons x nSamples].
%
% Optional name-value arguments:
%   'Zscore'      (true)  - z-score each neuron across time (recommended)
%   'SmoothWin'   (1)     - moving average window (samples); 1 = no smoothing
%   'Distance'    ('corr')- 'corr' (1-corr), 'euclidean', 'cosine', 'cityblock'
%   'MDSmethod'   ('classic') - 'classic' (cmdscale) or 'nonmetric' (mdscale)
%   'Dims'        (2)     - embedding dims to compute
%   'MaxK'        (8)     - max clusters to try for silhouette

% Parse inputs
p = inputParser;
addRequired(p,'sdf',@(x) ismatrix(x) && size(x,1)>1);
addParameter(p,'Zscore',true,@islogical);
addParameter(p,'SmoothWin',10,@(x) isnumeric(x) && x>=1);
addParameter(p,'Distance','corr',@(s) ismember(s,{'corr','euclidean','cosine','cityblock'}));
addParameter(p,'MDSmethod','classic',@(s) ismember(s,{'classic','nonmetric'}));
addParameter(p,'Dims',5,@(x) isnumeric(x) && x>=1);
addParameter(p,'MaxK',10,@(x) isnumeric(x) && x>=2);
parse(p,sdf,varargin{:});
opts = p.Results;

[nNeurons, nSamples] = size(sdf);

%% 1) Preprocess
data = double(sdf);

% Optional smoothing
if opts.SmoothWin > 1
    kernel = ones(1,opts.SmoothWin)/opts.SmoothWin;
    for i = 1:nNeurons
        data(i,:) = conv(data(i,:), kernel, 'same');
    end
end

% Optional z-score each neuron across time (emphasizes shape)
if opts.Zscore
    data = (data - mean(data,2))./std(data,[],2);
    % handle zero-variance rows
    zeroSD = std(data,[],2)==0;
    if any(zeroSD)
        data(zeroSD,:) = 0; % or leave as raw
    end
end

%% 2) Distance matrix (neurons x neurons)
switch opts.Distance
    case 'corr'    % correlation distance (1 - corr)
        R = corr(data');              % neurons x neurons
        D = 1 - R;
    case {'euclidean','cosine','cityblock'}
        D = squareform(pdist(data, opts.Distance));
    otherwise
        error('Unknown distance.');
end

% Ensure D is valid (non-negative, zero diagonal)
D = max(D,0);
D(1:size(D,1)+1:end) = 0;

%% 3) Run MDS
dims = opts.Dims;
switch opts.MDSmethod
    case 'classic'
        [Y, eigvals] = cmdscale(D, dims);
        stress = NaN;
        fprintf('Classic MDS (cmdscale). Eigenvalues (first %d): %s\n', dims, mat2str(eigvals(1:dims),3));
    case 'nonmetric'
        % mdscale wants condensed pdist OR full distance matrix; pass full
        % Use default options, try several starts to avoid local minima
        options = statset('MaxIter',300,'Display','final'); 
        [Y,stress] = mdscale(D, dims, 'Start','random','Criterion','stress', 'Options',options);
        fprintf('Nonmetric MDS (mdscale). Final stress: %.4f\n', stress);
end

%% 4) Shepard plot (Distance vs embedded distance) & stress (if available)
% Compute pairwise distances in embedded space
D_embedded = squareform(pdist(Y));
% Flatten upper triangle
iu = triu(true(nNeurons),1);
X = D(iu);   % original distances
Yd = D_embedded(iu); % embedded distances

figure('Name','Shepard plot'); 
scatter(X, Yd, 8, 'filled'); hold on;
plot([min(X) max(X)], [min(X) max(X)], '--k');
xlabel('Original distance'); ylabel('Embedded distance');
title('Shepard plot: original vs embedded distances');
grid on;

%% 5) Choose number of clusters (k-means on MDS coords) with silhouette
maxK = min(opts.MaxK, nNeurons-1);
sil_means = nan(maxK-1,1);
for k = 2:maxK
    idxk = kmeans(Y, k, 'Replicates', 20, 'Distance', 'sqeuclidean','Display','off');
    s = silhouette(Y, idxk);
    sil_means(k-1) = mean(s);
end

% Plot silhouette vs k
figure('Name','Silhouette for choosing k');
plot(2:maxK, sil_means, '-o','LineWidth',1.5);
xlabel('k (clusters)'); ylabel('Mean silhouette');
title('Silhouette vs k (use the k with max silhouette)');
grid on;

% Best k
[~, best_rel] = max(sil_means);
best_k = best_rel + 1;
fprintf('Recommended k (highest mean silhouette): %d\n', best_k);

% Final clustering
cluster_idx = kmeans(Y, best_k, 'Replicates', 50, 'Distance','sqeuclidean');

%% 6) Visualize MDS embedding colored by cluster
figure('Name','MDS embedding');
scatter3(Y(:,1), Y(:,2), Y(:,3), 40, cluster_idx, 'filled');
xlabel('MDS1'); ylabel('MDS2'); title('MDS embedding colored by cluster');
colorbar; colormap(jet);
text(Y(:,1)+0.01, Y(:,2), cellstr(num2str((1:nNeurons)')), 'FontSize',7);

%% 7) Compute cluster-average motifs & neuron contributions
cluster_mean = nan(best_k, nSamples);
for k = 1:best_k
    cluster_mean(k,:) = mean(sdf(cluster_idx==k, :), 1);   % use raw sdf for motif shape
end

% How well does the cluster mean explain each neuron? (R^2 per neuron)
rsq_per_neuron = nan(nNeurons,1);
for i = 1:nNeurons
    c = cluster_idx(i);
    y = sdf(i,:)';
    yhat = cluster_mean(c,:)';
    % simple R^2
    SSres = sum((y - yhat).^2);
    SStot = sum((y - mean(y)).^2);
    rsq_per_neuron(i) = 1 - SSres/SStot;
end

% Plot cluster means
figure('Name','Cluster-average motifs');
hold on;
for k = 1:best_k
    plot(-200:800, cluster_mean(k,:), 'LineWidth', 2);
end
legend(arrayfun(@(x) sprintf('Cluster %d',x), 1:best_k, 'UniformOutput', false));
xlabel('Time (samples)'); ylabel('SDF (Hz or raw)');
title('Cluster-average firing motifs');
grid on; xlim([-200 800])
hold off;

%% 8) Summary table
T = table((1:nNeurons)', cluster_idx, rsq_per_neuron, 'VariableNames', {'Neuron','Cluster','R2_to_cluster_mean'});
disp('First 20 mapping (Neuron -> Cluster, R^2 to cluster mean):');
disp(T(1:min(20,height(T)),:));

%% Optional: Save results
results.Y = Y;
results.D = D;
results.cluster_idx = cluster_idx;
results.cluster_mean = cluster_mean;
results.rsq_per_neuron = rsq_per_neuron;
results.sil_means = sil_means;
results.best_k = best_k;
assignin('base','mds_results',results); % pushes to workspace

fprintf('Done. Results saved to workspace variable ''mds_results''.\n');

end
