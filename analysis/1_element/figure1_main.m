


clu_exp_neurons = [22 24 36 70];
plot_ylim_example = {[0 15], [14 24], [0 10], [0 10]};
plot_ylim_pop = {[-1 7], [-2 6], [-6 6], [-3 3]};
plot_startxpos = [0.15 0.35 0.55 0.75];

color_pal.clu1 = [255 115 38]./255;
color_pal.clu2 = [255 25 77]./255;
color_pal.clu3 = [191 38 105]./255;
color_pal.clu4 = [112 42 140]./255;


% Clear variables starting with 'example_'
clear example_* single_unit_fig
for cluster_i = 1:4
    example_neuron_idx = neuron_class.cluster_idx.(['clu' int2str(cluster_i)])(clu_exp_neurons(cluster_i));
    xlim_vals = [-180 600];

    % Filter example neuron data for non-baseline and non-viol data points
    example_valid_idx = strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,5)) ,'nonviol') & ~strcmp(string(sdf_soundAlign_data{example_neuron_idx}(:,3)) ,'Baseline');
    example_sdf_in = cell2mat(sdf_soundAlign_data{example_neuron_idx}(example_valid_idx, 1));  % Extract SDFs for valid trials
    example_raster_in = sdf_soundAlign_data{example_neuron_idx}(example_valid_idx, 2);

    for trial_i = 1:size(example_sdf_in,1)
        example_sdf_in(trial_i,:) = smooth(example_sdf_in(trial_i,:),100)';
    end

    % Plot 1: Example neuron raster plot
    single_unit_fig(1, cluster_i) = gramm('x', example_raster_in);
    single_unit_fig(1, cluster_i).geom_raster('geom', {'point'});
    single_unit_fig(1, cluster_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));

    % Plot 2: Example neuron SDF plot
    single_unit_fig(2, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in);
    single_unit_fig(2, cluster_i).stat_summary();
    single_unit_fig(2, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_example{cluster_i});
    single_unit_fig(2, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));

    % Plot 3: Population SDF plot
    single_unit_fig(3, cluster_i) = gramm('x', ops.sound_sdf_window, 'y', inputSDF{1}(neuron_class.(['cluster' int2str(cluster_i)]),:));
    single_unit_fig(3, cluster_i).stat_summary();
    single_unit_fig(3, cluster_i).axe_property('XLim', xlim_vals, 'YLim', plot_ylim_pop{cluster_i});
    single_unit_fig(3, cluster_i).set_color_options('map',color_pal.(['clu' int2str(cluster_i)]));


    % Figure layout ////////////////////////////////////////////////////////////
    single_unit_fig(1, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.8 0.15 0.1], 'legend', false, 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(2, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.5 0.15 0.25], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(3, cluster_i).set_layout_options('Position', [plot_startxpos(cluster_i) 0.2 0.15 0.25], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);


    single_unit_fig(1, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(2, cluster_i).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 0, 'style', 'k-');
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 413, 'style', 'k--');
    single_unit_fig(3, cluster_i).geom_vline('xintercept', 563, 'style', 'k-');

end

% Create the figure and draw the plots
figure('Renderer', 'painters', 'Position', [100 100 1200 600]);
single_unit_fig.draw();

