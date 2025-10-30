% Identify significant neurons from auditory and frontal brain regions
sig_neurons = []; 
sig_neurons = [neuron_class.auditory.all; neuron_class.frontal.all]; % Combine auditory and frontal neurons

% Clear and define input spike density function (SDF) data
clear inputSDF inputSDF_plot
inputSDF = {normFR_in.norm_fr_soundAll(sig_neurons,:)}; % Get normalized firing rates for all sound trials

% Smooth SDF data for each neuron
for neuron_i = 1:size(inputSDF{1},1)
    inputSDF{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),10); % Light smoothing for clustering
    inputSDF_plot{1}(neuron_i,:) = smooth(inputSDF{1}(neuron_i,:),50); % Heavier smoothing for plotting
end

sdf = inputSDF{1};
mds_neuron_pipeline(sdf)

[Y,stress] = mdscale(mds_results.D,5);  % D = distance matrix from SDFs
fprintf('Stress = %.3f\n', stress);

plot_lineclust_color = cool(mds_results.best_k);


% Plot cluster means
figure('Renderer', 'painters', 'Position', [100 100 1250 800]);
hold on;

% Plot average SDF for each cluster
for k = 1:mds_results.best_k
    a = subplot(2,ceil(mds_results.best_k/2),k); hold on
    plot(-200:800, mds_results.cluster_mean(k,:),'Color',plot_lineclust_color(k,:));
    vline([0 563], 'k-');  % Stimulus window
    vline([413], 'k--');   % Possible event of interest
    xlim([-150 600])
    title(['Cluster ' int2str(k)])
end


glm_encoding_flag(:,[1 2 3 4 5])
