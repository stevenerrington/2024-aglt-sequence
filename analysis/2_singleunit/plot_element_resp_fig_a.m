function plot_element_resp_fig_a(neuron_example, sdf_soundAlign_data, ops, pop_neurons_in, normFR_in)
    % Define list of sound stimuli
    sound_list = {'A', 'C', 'D', 'F', 'G'};

    % Set axis limits for the plots
    xlim_vals = [-100 600];
    ylim_vals = [0 60];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Filter example neuron data for non-baseline and non-viol data points
    example_valid_idx = strcmp(string(sdf_soundAlign_data{neuron_example}(:,5)) ,'nonviol') & ~strcmp(string(sdf_soundAlign_data{neuron_example}(:,3)) ,'Baseline') & cell2mat(sdf_soundAlign_data{neuron_example}(:,9)) == 1;
    example_sdf_in = cell2mat(sdf_soundAlign_data{neuron_example}(example_valid_idx, 1));  % Extract SDFs for valid trials

    % Smooth each row (trial) of SDF data with a smoothing window of 50
    for i = 1:size(example_sdf_in, 1)
        example_sdf_in(i, :) = smooth(example_sdf_in(i, :), 50);
    end

    % Extract raster, sound condition, and order condition data for example neuron
    example_raster_in = sdf_soundAlign_data{neuron_example}(example_valid_idx, 2);
    example_sound_cond = sdf_soundAlign_data{neuron_example}(example_valid_idx, 3);
    example_order_cond = sdf_soundAlign_data{neuron_example}(example_valid_idx, 4);

    % Plot 2: Example neuron SDF plot
    single_unit_fig(1, 1) = gramm('x', ops.sound_sdf_window, 'y', example_sdf_in, 'color', example_order_cond);
    single_unit_fig(1, 1).stat_summary();
    single_unit_fig(1, 1).axe_property('XLim', xlim_vals, 'YLim', ops.plot.example_yaxis);
    single_unit_fig(1, 1).geom_vline('xintercept',[0 413 563]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Concatenate normalized SDF data for population neurons across different sounds
    pop_sdf_in = [normFR_in.norm_fr_soundA(pop_neurons_in, :);
        normFR_in.norm_fr_soundC(pop_neurons_in, :);
        normFR_in.norm_fr_soundD(pop_neurons_in, :);
        normFR_in.norm_fr_soundF(pop_neurons_in, :);
        normFR_in.norm_fr_soundG(pop_neurons_in, :)];

    % Create a condition label vector for population data
    pop_sound_cond = [repmat({'A'}, length(pop_neurons_in), 1);
        repmat({'C'}, length(pop_neurons_in), 1);
        repmat({'D'}, length(pop_neurons_in), 1);
        repmat({'F'}, length(pop_neurons_in), 1);
        repmat({'G'}, length(pop_neurons_in), 1)];


    % Plot 3: Population SDF plot
    single_unit_fig(1, 2) = gramm('x', ops.sound_sdf_window, 'y', pop_sdf_in, 'color', pop_sound_cond);
    single_unit_fig(1, 2).stat_summary();
    single_unit_fig(1, 2).axe_property('XLim', xlim_vals, 'YLim',  ops.plot.pop_yaxis);
    single_unit_fig(1, 2).geom_vline('xintercept',[0 413 563]);

    % Create the figure and draw the plots
    figure('Renderer', 'painters', 'Position', [100 100 900 300]);
    single_unit_fig.draw();
end