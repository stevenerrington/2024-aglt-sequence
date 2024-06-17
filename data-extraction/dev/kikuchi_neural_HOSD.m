
%% Setup environment
% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Admin --------------------------------------------------------------
monkey = 'troy'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2021-11-05-AGLt'; % Experimental raw data
task = 'agl_t'; % Experiment type [agl, opto]
session_n = '0001'; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name

%% Extract broadband signal

for ch_n = 1:n_channels
    clear filepart_name NCSpath spk_ncs_in

    filepart_name = ['CSC' int2str(ch_n) '_' session_n];
    NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

    spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
    spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
end


%% Run HOSD Spike Detection Algorithm

% Define algorithm parameters
params.lowpass = 4000; % lowpass cutoff ( Hz )
params.highpass = 200; % highpass cutoff ( Hz )
params.limit_memory = 5e8; 
params.check_features = false;
params.randseed = 1;
saveAllRaw = 0; % boolean to keep spikes data (rather than just the clustered data)
saveFigFile = 0;

% Setup data for function
clear data
data.dat = spk_ncs_out;
data.fs = 32000;


% - 2024-06-14: Ran this code and it took an age (stopped after 45 minutes
% to access MATLAB for other tasks). Will run again over the weekend?
mat_savefile = fullfile(dirs.mat_data,['troy-agl_t-2021-11-05-HOSD.mat']);

if exist(mat_savefile) == 2
    load(fullfile(dirs.mat_data,['troy-agl_t-2021-11-05-HOSD.mat']))
    fprintf('Loading %s ... \n', 'troy-agl_t-2021-11-05-HOSD.mat')
else
    [spikes,hos,znrm] = HOSD_spike_detection(data,params); % Run extraction
    spike_cluster = sort_spikes(spikes); % Sort spikes
    save(mat_savefile,'spikes','hos','znrm','spike_cluster','-v7.3')
end


%% Plot summary sorting figure
fig = plot_clusters_2(spike_cluster);
set(fig,'renderer','painters','Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,['C:\KIKUCHI-LOCAL\script\kikuchi-data\data-extraction\doc\troy-agl_t-2021-11-05-HOSD','.pdf'],'-r400','-bestfit','-dpdf')

%% Restructure spike data for future analysis
n_clusters = spike_cluster.Nclust; % Find number of spikes

clear spikes_out
for cluster_i = 1:n_clusters

    % Get spike times
    spikes_out.time.(['DSP' int2str(cluster_i)]) = ...
        spike_cluster.spike_times(find(spike_cluster.cl == cluster_i))*1000;

    % Get spike waveforms
     spikes_out.waveform.(['WAV' int2str(cluster_i)]) = ...
        spike_cluster.avg_waves(:,cluster_i)';
end

%% Align spikes on events

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table events
event_table = get_event_table(ops);

% Align data to event
% Define parameters
ops.event_port = 2;
events = get_agl_t_trials(event_table, ops);
aligntime = events.stimulusOnset_ms;

ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes_out,aligntime,ops);

%% Produce spike figures for each neuron

% Get list of neurons
names = fieldnames( spikes_out.time );

% Define times for figures
xlim_vals = [-1000 5000];
ylim_vals = [0 150];

for ch_i = 1:length(names)
    clear spk_figure_idv fig

    % Get channel label
    ch = names{ch_i}(4:end);

    % Plot raster
    spk_figure_idv(1,1)=gramm('x',raster.(['DSP' ch]),'color',events.cond_label);
    spk_figure_idv(1,1).geom_raster('geom',{'line'});
    spk_figure_idv(1,1).axe_property('XLim',xlim_vals,'YLim',[0 size(raster.(['DSP' ch]),1)]);

    % Plot SDF
    spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',events.cond_label);
    spk_figure_idv(2,1).stat_summary();
    spk_figure_idv(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);

    % Plot waveform average
    spk_figure_idv(3,1)=gramm('x',[1:length(spikes_out.waveform.(['WAV' ch]))],'y',spikes_out.waveform.(['WAV' ch]));
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

    % Save figure
    set(fig,'renderer','painters','Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,['C:\KIKUCHI-LOCAL\script\kikuchi-data\data-extraction\doc\troy-agl_t-2021-11-05-HOSD-',['DSP' ch],'.pdf'],'-r400','-bestfit','-dpdf')
    close all
end
