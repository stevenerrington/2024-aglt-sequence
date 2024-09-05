function plot_element_resp_fig(sdf_soundAlign_data, neuron_example, pop_neurons_in, normFR_in, plot_beta_weight, ops)
    % Function to plot example neuron raster, example neuron SDF, and population SDF
    % Parameters:
    %   - sdf_soundAlign_data: Cell array containing sound-aligned spike density function data
    %   - aud_neuron_example: Index of the example neuron to plot
    %   - aud_mod_neurons_pos: Indices of positively modulated neurons
    %   - aud_mod_neurons_mixed: Indices of mixed modulated neurons
    %   - normFR_in: Struct containing normalized firing rates
    %   - plot_beta_weight: Cell array containing beta weights for different sounds
    %   - ops: Struct containing additional plotting options such as sound SDF window
    
    % Define list of sound stimuli
    sound_list = {'A', 'C', 'D', 'F', 'G'};

    % Set axis limits for the plots
    xlim_vals = [-100 500];
    ylim_vals = [0 60];

    % Clear any existing figure data
    clear single_unit_fig;

    % Clear variables starting with 'example_'
    clear example_*

    % Filter example neuron data for non-baseline and non-viol data points
    example_valid_idx = strcmp(string(sdf_soundAlign_data{neuron_example}(:,5)) ,'nonviol') & ~strcmp(string(sdf_soundAlign_data{neuron_example}(:,3)) ,'Baseline');
    example_sdf_in = cell2mat(sdf_soundAlign_data{neuron_example}(example_valid_idx, 1));  % Extract SDFs for valid trials

    % Smooth each row (trial) of SDF data with a smoothing window of 50
    for i = 1:size(example_sdf_in, 1)
        example_sdf_in(i, :) = smooth(example_sdf_in(i, :), 50);
    end

    % Extract raster, sound condition, and order condition data for example neuron
    example_raster_in = sdf_soundAlign_data{neuron_example}(example_valid_idx, 2);
    example_sound_cond = sdf_soundAlign_data{neuron_example}(example_valid_idx, 3);
    example_order_cond = sdf_soundAlign_data{neuron_example}(example_valid_idx, 4);

    % Clear any existing population data variables
    clear pop_

    % Concatenate normalized SDF data for population neurons across different sounds
    pop_sdf_in = [normFR_in.norm_fr_soundA(pop_neurons_in, :);
                  normFR_in.norm_fr_soundC(pop_neurons_in, :);
                  normFR_in.norm_fr_soundD(pop_neurons_in, :);
                  normFR_in.norm_fr_soundF(pop_neurons_in, :);
                  normFR_in.norm_fr_soundG(pop_neurons_in, :)];

    % Smooth each row (neuron) of the population SDF data with a smoothing window of 50
    for i = 1:size(pop_sdf_in, 1)
        pop_sdf_in(i, :) = smooth(pop_sdf_in(i, :), 1);
    end

    % Concatenate beta weights for population neurons across different sounds
    pop_beta_in = [plot_beta_weight{1}(pop_neurons_in, :);
                   plot_beta_weight{2}(pop_neurons_in, :);
                   plot_beta_weight{3}(pop_neurons_in, :);
                   plot_beta_weight{4}(pop_neurons_in, :);
                   plot_beta_weight{5}(pop_neurons_in, :)];

    % Create a condition label vector for population data
    pop_sound_cond = [repmat({'A'}, length(pop_neurons_in), 1);
                      repmat({'C'}, length(pop_neurons_in), 1);
                      repmat({'D'}, length(pop_neurons_in), 1);
                      repmat({'F'}, length(pop_neurons_in), 1);
                      repmat({'G'}, length(pop_neurons_in), 1)];

    % Plot 1: Example neuron raster plot
    single_unit_fig(1, 1) = gramm('x', example_raster_in, 'color', example_sound_cond);
    single_unit_fig(1, 1).geom_raster('geom', {'point'});
    single_unit_fig(1, 1).axe_property('XLim', xlim_vals);

    % Plot 2: Example neuron SDF plot
    single_unit_fig(2, 1) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in, 'color', example_sound_cond);
    single_unit_fig(2, 1).stat_summary();
    single_unit_fig(2, 1).axe_property('XLim', xlim_vals, 'YLim', ops.plot.example_yaxis);

    % Plot 3: Population SDF plot
    single_unit_fig(3, 1) = gramm('x', ops.sound_sdf_window, 'y', pop_sdf_in, 'color', pop_sound_cond);
    single_unit_fig(3, 1).stat_summary();
    single_unit_fig(3, 1).axe_property('XLim', xlim_vals, 'YLim',  ops.plot.pop_yaxis);

    % Plot 4: Population beta weight plot
    single_unit_fig(4, 1) = gramm('x', -150:10:750, 'y', pop_beta_in, 'color', pop_sound_cond);
    single_unit_fig(4, 1).stat_summary();
    single_unit_fig(4, 1).axe_property('XLim', xlim_vals, 'YLim', ops.plot.beta_yaxis);

    % Figure setup ////////////////////////////////////////////////////////
    single_unit_fig(1, 1).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(2, 1).axe_property('XTick', [], 'XColor', [1 1 1]);
    single_unit_fig(3, 1).axe_property('XTick', [], 'XColor', [1 1 1]);

    single_unit_fig(1, 1).geom_vline('xintercept', 0, 'style', 'k-');
    single_unit_fig(2, 1).geom_vline('xintercept', 0, 'style', 'k-');
    single_unit_fig(3, 1).geom_vline('xintercept', 0, 'style', 'k-');
    single_unit_fig(4, 1).geom_vline('xintercept', 0, 'style', 'k-');

    single_unit_fig(1, 1).set_names('y', 'Trials');
    single_unit_fig(2, 1).set_names('x', '', 'y', 'FR (spk/sec)');
    single_unit_fig(3, 1).set_names('x', '', 'y', 'FR (norm)');
    single_unit_fig(4, 1).set_names('x', 'Time from event (ms)', 'y', 'Beta weight');

    % Figure layout ////////////////////////////////////////////////////////////
    single_unit_fig(1, 1).set_layout_options('Position', [0.2 0.85 0.75 0.1], 'legend', false, 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(2, 1).set_layout_options('Position', [0.2 0.6 0.75 0.2], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(3, 1).set_layout_options('Position', [0.2 0.35 0.75 0.2], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);
    single_unit_fig(4, 1).set_layout_options('Position', [0.2 0.1 0.75 0.2], 'margin_height', [0.00 0.00], 'margin_width', [0.0 0.00], 'redraw', false);

    % Create the figure window with specific size and renderer settings
    figure('Renderer', 'painters', 'Position', [100 100 300 600]);
    single_unit_fig.draw();  % Draw all plots in the figure
end