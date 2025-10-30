function pc_out = perform_pca_and_plot(neurons_in, pca_sdf_out)
    % Function to perform PCA analysis and plot the results

    % Clear variables from workspace to avoid conflicts
    clear sdf_in_regular sdf_in_shuffled

    % Define PCA window
    pca_window = -100:10:2750;  % Can be overridden by input if needed
    sdf_in_regular = pca_sdf_out(neurons_in, 1000 + pca_window);

    sdf_in_regular = rmmissing(sdf_in_regular, 1);


    % Shuffle data for comparison
    for neuron_i = 1:size(sdf_in_regular, 1)
        sdf_in_shuffled(neuron_i, :) = sdf_in_regular(neuron_i, randperm(size(sdf_in_regular, 2)));
    end

    % Clear previous PCA variables
    clear pcs* var*

    % Perform PCA on regular and shuffled data
    [coeff, pcs, latent, ~, var_exp, ~] = pca(sdf_in_regular');
    [coeff_shuffled, pcs_shuffled, latent_shuffled, ~, var_exp_shuffled, ~] = pca(sdf_in_shuffled');

    pc_out.obs.coeff = coeff;
    pc_out.obs.pcs = pcs;
    pc_out.obs.latent = latent;
    pc_out.obs.var_exp = var_exp;

    pc_out.shuffled.coeff = coeff_shuffled;
    pc_out.shuffled.pcs = pcs_shuffled;
    pc_out.shuffled.latent = latent_shuffled;
    pc_out.shuffled.var_exp = var_exp_shuffled;

    pc_out.window = pca_window;

    % Define color schemes for plotting
    pca_color = [162 13 30; 172 40 55; 183 67 80; 193 94 105; 203 121 130; 214 147 155; 224 174 180; 234 201 205] ./ 255;
    pca_color_shuff = [140 140 140; 153 153 153; 166 166 166; 178 178 178; 191 191 191; 204 204 204; 217 217 217] ./ 255;

    % Display variance explained
    % disp([var_exp(1:20)'; var_exp_shuffled(1:20)']);

    % Clear PCA components to avoid conflict
    clear pc1 pc2 pc3
    pc1 = pcs(:, 1); pc2 = pcs(:, 2); pc3 = pcs(:, 3);
    pc4 = pcs(:, 4); pc5 = pcs(:, 5); pc6 = pcs(:, 6);

    % Define onset time index
    onset_time_idx = find(pca_window == 0);

    % % Create figure for plotting
    % figuren('Renderer', 'painters', 'Position', [100 100 1000 250]); hold on;
    % 
    % % Variance explained plot
    % n_vars = 6;
    % colorscale = abs(flipud(cbrewer('seq', 'PuRd', n_vars)));
    % colorscale_shuf = abs(flipud(cbrewer('seq', 'Greys', n_vars)));
    % 
    % nsubplot(3, 10, [1 2 3], [1 2]); hold on
    % b_obs = bar([1:n_vars], var_exp([1:n_vars]), 'LineStyle', 'None', 'FaceAlpha', 0.5);
    % b_shuf = bar([1:n_vars], var_exp_shuffled([1:n_vars]), 'LineStyle', 'None', 'FaceAlpha', 0.5);
    % xlabel('Principal Component')
    % ylabel('Cumulative Variance Explained (%)')
    % xticks([1:1:n_vars]);
    % ylim([0 20]);
    % 
    % for var_i = 1:n_vars
    %     b_obs.FaceColor = 'flat';
    %     b_shuf.FaceColor = 'flat';
    %     b_obs.CData(var_i, :) = colorscale(var_i, :);
    %     b_shuf.CData(var_i, :) = colorscale_shuf(1, :);
    % end
    % 
    % % Define sound times
    % sound_times = [0, 563, 1126, 1689, 2252];
    % 
    % closest_vals = zeros(size(sound_times));
    % closest_idx  = zeros(size(sound_times));
    % 
    % for i = 1:length(sound_times)
    %     [~, idx] = min(abs(pca_window - sound_times(i)));
    %     closest_vals(i) = pca_window(idx);
    %     closest_idx(i)  = idx;
    % end
    % 
    % % PCA x time plots for PC1, PC2, PC3
    % nsubplot(3, 10, [1], [4 5 6]); hold on
    % plot(pca_window, pc1, 'color', colorscale(1, :), 'LineWidth', 1.5);
    % xlim([min(pca_window) max(pca_window)]);
    % vline(sound_times, 'k');
    % set(gca, 'Xcolor', [1 1 1]);
    % 
    % nsubplot(3, 10, [2], [4 5 6]); hold on
    % plot(pca_window, pc2, 'color', colorscale(2, :), 'LineWidth', 1.5);
    % vline(sound_times, 'k'); ylabel('PC');
    % xlim([min(pca_window) max(pca_window)]);
    % set(gca, 'Xcolor', [1 1 1]);
    % 
    % nsubplot(3, 10, [3], [4 5 6]); hold on
    % plot(pca_window, pc3, 'color', colorscale(3, :), 'LineWidth', 1.5);
    % xlabel('Time from stimulus onset (ms)');
    % xlim([min(pca_window) max(pca_window)]);
    % vline(sound_times, 'k');
    % 
    % % 3D PCA Plot
    % nsubplot(3, 10, [1 2 3], [8 9 10]); hold on
    % color_line3(pc1, pc2, pc3, pca_window, 'LineWidth', 2);
    % view(34.2409, 7.6800);
    % xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    % scatter3(pc1(closest_idx), pc2(closest_idx), pc3(closest_idx), 100, [0 0 0], '^', 'filled');
    % grid on;
end