ft_defaults

ops = struct();
ops.rootZ = 'H:\Desktop\acute_test\ks\marmo-acute-meeko-2025-04-17';
ops.bin_file = 'H:\Desktop\acute_test\bin\marmo-acute-meeko-2025-04-17.dat';
ops.nCh = 16;
ops.fs = 32000;

[spikes] = phy2mat(ops);
[spk_info] = phyinfo2mat(ops);


%% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = 'H:\Desktop\acute_test\raw\'; ops.filename = '2025-04-16_13-30-59'; ops.session_n = '';
clear event_table_raw event_table
event_table_raw = get_event_table(ops);
ops.event_port = 2;
event_table = event_table_raw(event_table_raw.port == ops.event_port, :);
event_table = event_table(event_table.value > 150, :);

letters = char('A' + (0:24-1)); % Creates A to X
map = containers.Map(179:202, cellstr(letters')); 

for i = 1:size(event_table,1)
    event_table.sound_cond{i} = map(event_table.value(i));
end

%% LFP data 
ops.timewin = [-500:1000];
% Local field potential data -------------------------------------------------------
filelabels_lfp = get_ncs_filelabel(fullfile(ops.dirs.raw_data,[ops.filename '\']), ['LFP1' '.ncs'],16);
lfp = ft_read_neuralynx_interp(filelabels_lfp);
lfp = lfp.trial{1};

lfp_aligned = get_lfp_aligned(lfp,event_table.timestamp_ms,ops);

trials = 1:length(event_table.timestamp_ms);


figuren('Renderer', 'painters', 'Position', [100 100 900 800]); hold on

for ch_i = 1:16
    lfp_in = [];
    lfp_in = lfp_aligned.(['lfp_' int2str(ch_i)]);


    if ismember(ch_i,1:16) ...
            color_line_value = 'r';
    else
        color_line_value = 'b';
    end

    plot(ops.timewin, nanmean(lfp_in(trials,:))+(20*ch_i),'color',color_line_value)


    lfp_out(ch_i,:) =  nanmean(lfp_in(trials,:));
end

set(gca,'YDir','Reverse')
hline(16.5*20,'k')
vline(0,'k')
xlabel('Time from sound onset (ms)')
clear electrode_mean

figuren;
plot(ops.timewin, nanmean(lfp_out([2:4,7:10,12:16],:)))

%% Spikes data
ops.timewin = -500:1000;
ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(spikes,event_table.timestamp_ms+200,ops);

% Get list of neurons
names = fieldnames( spikes.time );

% Define times for figures
xlim_vals = [-100 500];
ylim_vals = [0 10];

for ch_i = 1:length(names)
    clear spk_figure_idv fig

    % Get channel label
    ch = names{ch_i}(4:end);

    % Plot raster
    spk_figure_idv(1,1)=gramm('x',raster.(['DSP' ch]));%,'color',event_table.sound_cond);
    spk_figure_idv(1,1).geom_raster('geom',{'line'});
    spk_figure_idv(1,1).axe_property('XLim',xlim_vals);

    % Plot SDF
    spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]));%,'color',event_table.sound_cond);
    spk_figure_idv(2,1).stat_summary();
    spk_figure_idv(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);

    % Plot waveform average
    spk_figure_idv(3,1)=gramm('x',[1:length(spikes.waveform.(['WAV' ch]))],'y',spikes.waveform.(['WAV' ch]));
    spk_figure_idv(3,1).stat_summary();
    spk_figure_idv(3,1).axe_property('XLim',[-41 40],'YLim',[-35 35]);

    % Figure properties
    spk_figure_idv(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
    spk_figure_idv(3,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
    spk_figure_idv(1,1).geom_vline('xintercept',0,'style','k-');
    spk_figure_idv(2,1).geom_vline('xintercept',0,'style','k-');
    spk_figure_idv(1,1).set_names('y','Trials');
    spk_figure_idv(2,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');

    % Figure layout
    spk_figure_idv(1,1).set_layout_options...
        ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv(2,1).set_layout_options...
        ('Position',[0.1 0.1 0.8 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv(3,1).set_layout_options...
        ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend',false,...
        'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
        'margin_width',[0.0 0.00],...
        'redraw',false);

    spk_figure_idv.set_title(['DSP' ch]);

    % Draw figure
    fig = figure('Renderer', 'painters', 'Position', [100 100 700 600]);
    spk_figure_idv.draw();


end


