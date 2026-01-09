function tsne_neuron_pipeline(sdf, varargin)
% TSNE_NEURON_PIPELINE  t-SNE + clustering pipeline for neuron SDFs
%   tsne_neuron_pipeline(sdf) where sdf is [nNeurons x nSamples]
%
% Optional name-value args:
%   'SmoothWin'    (5)    - moving average window in samples (1 = no smooth)
%   'Zscore'       (true) - z-score each neuron across time (emphasize shape)
%   'PCAreds'      (30)   - # of PCA components to keep before t-SNE (0 = no PCA)
%   'Perplexity'   (30)   - t-SNE perplexity (try 5-50)
%   'Dims'         (2)    - output embedding dims (2 or 3)
%   'MaxK'         (8)    - max k to try when choosing k with silhouette
%   'RandSeed'     (42)   - seed for reproducibility
%   'RunDBSCAN'    (true) - run DBSCAN (eps/minpts need tuning)
%
% Example:
%   tsne_neuron_pipeline(sdf, 'SmoothWin',5, 'PCAreds',20, 'Perplexity',30);

p = inputParser;
addRequired(p,'sdf',@(x) ismatrix(x) && size(x,1)>1);
addParameter(p,'SmoothWin',5,@(x) isnumeric(x) && x>=1);
addParameter(p,'Zscore',true,@islogical);
addParameter(p,'PCAreds',30,@(x) isnumeric(x) && x>=0);
addParameter(p,'Perplexity',30,@(x) isnumeric(x) && x>0);
addParameter(p,'Dims',5,@(x) isnumeric(x) && ismember(x,[2,3]));
addParameter(p,'MaxK',8,@(x) isnumeric(x) && x>=2);
addParameter(p,'RandSeed',42,@isnumeric);
addParameter(p,'RunDBSCAN',true,@islogical);
parse(p,sdf,varargin{:});
opts = p.Results;

% Setup
data_raw = double(sdf);
[nNeurons, nSamples] = size(data_raw);
fprintf('Running t-SNE pipeline on %d neurons × %d samples\n', nNeurons, nSamples);

%% 1) Preprocess
data = data_raw;

% Smooth
if opts.SmoothWin > 1
    kernel = ones(1,opts.SmoothWin)/opts.SmoothWin;
    for i = 1:nNeurons
        data(i,:) = conv(data(i,:), kernel, 'same');
    end
end

% Z-score each neuron across time (emphasize shape)
if opts.Zscore
    mu = mean(data,2);
    sd = std(data,[],2);
    sd(sd==0) = 1; % avoid divide-by-zero
    data = (data - mu) ./ sd;
end

% Optional PCA pre-reduction (recommended for speed/noise)
if opts.PCAreds > 0 && opts.PCAreds < min([nNeurons, nSamples])
    [~, score, ~, ~, explained] = pca(data);
    pcaDims = min(opts.PCAreds, size(score,2));
    X = score(:,1:pcaDims); % neurons x pcaDims
    fprintf('PCA pre-reduction: keeping %d PCs (%.1f%% variance)\n', pcaDims, sum(explained(1:pcaDims)));
else
    X = data; % use time samples directly (may be large)
end

%% 2) t-SNE
rng(opts.RandSeed);
disp('Running t-SNE (this may take a moment)...');

% Choose distance metric appropriately: 'cosine' or 'correlation' often good for shapes
tsneOpts = {'Algorithm','barneshut','NumDimensions',opts.Dims,'Perplexity',opts.Perplexity,'Distance','cosine','Standardize',false};
% Note: for older MATLAB versions remove 'Algorithm','barneshut' if unsupported.
Y = tsne(X, tsneOpts{:});  % Y: nNeurons x opts.Dims

%% 3) Plot embedding
figure('Name','t-SNE embedding');
if opts.Dims == 2
    scatter(Y(:,1),Y(:,2),36,'filled');
    xlabel('tSNE1'); ylabel('tSNE2');
else
    scatter3(Y(:,1),Y(:,2),Y(:,3),36,'filled');
    xlabel('tSNE1'); ylabel('tSNE2'); zlabel('tSNE3');
end
title(sprintf('t-SNE embedding (perplexity=%g)', opts.Perplexity));
grid on;

%% 4) Choose k via silhouette (kmeans in t-SNE space)
maxK = min(opts.MaxK, nNeurons-1);
if maxK >= 2
    sil_means = nan(maxK-1,1);
    k_range = 2:maxK;
    for ki = 1:numel(k_range)
        k = k_range(ki);
        idxk = kmeans(Y, k, 'Replicates', 20, 'Distance', 'sqeuclidean', 'Display','off');
        s = silhouette(Y, idxk);
        sil_means(ki) = mean(s);
    end

    % Plot silhouette vs k
    figure('Name','Silhouette for k selection');
    plot(k_range, sil_means, '-o','LineWidth',1.5);
    xlabel('k (clusters)'); ylabel('Mean silhouette');
    title('Mean silhouette across k (higher = better)');
    grid on;

    [~, best_rel] = max(sil_means);
    best_k = k_range(best_rel);
    fprintf('Recommended k by silhouette: %d (mean silhouette = %.3f)\n', best_k, sil_means(best_rel));
else
    best_k = 2;
    fprintf('Too few neurons to explore many k; defaulting best_k = 2\n');
end

% Final kmeans with best_k
cluster_idx = kmeans(Y, best_k, 'Replicates', 50, 'Distance','sqeuclidean');

% Overlay embedding with clusters
figure('Name','t-SNE embedding colored by kmeans cluster');
cmap = jet(best_k);  % best_k × 3 RGB values
gscatter(Y(:,1), Y(:,2), cluster_idx, cmap, 'o', 8);
xlabel('tSNE1'); ylabel('tSNE2'); title(sprintf('t-SNE embedding (k=%d)', best_k));
grid on;

%% 5) DBSCAN (optional) - density based; eps needs tuning
if opts.RunDBSCAN
    % Heuristic: set MinPts = max(3, round(log(nNeurons)))
    MinPts = max(3, round(log(nNeurons)));
    % A starting eps: median of 5-nearest neighbor distances in Y
    kNN = min(nNeurons-1, 5);
    Dp = pdist2(Y,Y);
    sortedD = sort(Dp,2);
    d5 = sortedD(:, kNN+1); % kNN+1 because first column is zero (self)
    eps0 = median(d5);
    dbidx = dbscan(Y, eps0, MinPts);
    fprintf('DBSCAN run with eps=%.4g, MinPts=%d; found %d clusters (0 = noise)\n', eps0, MinPts, numel(unique(dbidx(dbidx>0))));
    figure('Name','t-SNE embedding colored by DBSCAN');
    cmap = parula(numel(unique(dbidx)));
    gscatter(Y(:,1), Y(:,2), dbidx, cmap, 'o', 8);
    title(sprintf('DBSCAN (eps=%.3g, MinPts=%d)', eps0, MinPts));
    xlabel('tSNE1'); ylabel('tSNE2'); grid on;
else
    dbidx = [];
end

%% 6) Compute cluster-average motifs & R^2 mapping (use raw sdf for motif shapes)
cluster_mean = nan(best_k, nSamples);
for k = 1:best_k
    cluster_members = find(cluster_idx==k);
    if isempty(cluster_members)
        cluster_mean(k,:) = nan(1,nSamples);
    else
        cluster_mean(k,:) = mean(data_raw(cluster_members, :), 1); % raw SDFs
    end
end

rsq_per_neuron = nan(nNeurons,1);
for i = 1:nNeurons
    c = cluster_idx(i);
    if any(isnan(cluster_mean(c,:)))
        rsq_per_neuron(i) = NaN;
    else
        y = data_raw(i,:)';
        yhat = cluster_mean(c,:)';
        SSres = sum((y - yhat).^2);
        SStot = sum((y - mean(y)).^2);
        rsq_per_neuron(i) = 1 - SSres/SStot;
    end
end

% Plot cluster-mean motifs
figure('Name','Cluster-average motifs');
hold on;
for k = 1:best_k
    plot(cluster_mean(k,:), 'LineWidth', 2);
end
legend(arrayfun(@(x) sprintf('Cluster %d',x), 1:best_k, 'UniformOutput', false), 'Location','best');
xlabel('Time (samples)'); ylabel('SDF (raw)');
title('Cluster-average firing motifs');
grid on; hold off;

%% 7) Summary table + save results
T = table((1:nNeurons)', cluster_idx, rsq_per_neuron, 'VariableNames', {'Neuron','Cluster','R2_to_cluster_mean'});
disp('First 20 neuron -> cluster mapping (Neuron, Cluster, R2):');
disp(T(1:min(20,height(T)),:));

results.Y = Y;
results.cluster_idx = cluster_idx;
results.cluster_mean = cluster_mean;
results.rsq_per_neuron = rsq_per_neuron;
results.dbscan_idx = dbidx;
results.params = opts;
assignin('base','tsne_results',results);

fprintf('Done. Results saved to workspace variable ''tsne_results''. Inspect tsne_results for details.\n');
end
