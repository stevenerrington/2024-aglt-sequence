
%% Workspace configuration and setup
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'beans'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2023-11-03-Marmo'; % Experimental raw data
task = 'marmo'; % Experiment type [agl, opto]
session_n = ''; % Experimental file tag

% Define experimental/data dirphyectories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

spike_status = 'phy';

switch spike_status
    
    % Convert .ncs files into a binary file for use in Kilosort --------/
    case 'bin' 
        for ch_n = 1:n_channels
            clear filepart_name NCSpath spk_ncs_in

            filepart_name = ['CSC' int2str(ch_n)];
            NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

            spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
            spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
        end

        % Create a binary file and export the restructure broadband data
        clear bin_out_file
        bin_out_file = fopen([dirs.bin_data outfile_name '.dat'],'wb');
        fwrite(bin_out_file,spk_ncs_out,'int16');
        fclose(bin_out_file);

        % Run kilosort
        mkdir(fullfile(dirs.kilosort,outfile_name));
        % - Run in python: to be integrated here.

    % Import data phy-curated Kilosort data ---------------------------/   
    case 'phy'  
        % Import phy curated data
        ops = struct();
        ops.rootZ = fullfile(dirs.kilosort,outfile_name);
        ops.bin_file = [dirs.bin_data outfile_name '.dat'];
        ops.nCh = n_channels;
        ops.fs = 32000;

        [spikes] = phy2mat(ops);
        [spk_info] = phyinfo2mat(ops);
end

% Local field potential data -------------------------------------------------------
filelabels_lfp = get_ncs_filelabel(fullfile(dirs.raw_data,[exp_filename '\']), ['LFP1' '.ncs'],32);
lfp = ft_read_neuralynx_interp(filelabels_lfp);
lfp = lfp.trial{1};

% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table agl_t_event
event_table = get_event_table(ops);

% Event extraction
clear agl_t_event aligntime
ops.event_port = 2;
agl_t_event = get_agl_t_trials(event_table, ops);
aligntime = agl_t_event.stimulusOnset_ms ;

% Get aligned neural activity ------------------------------------------
ops.timewin = -2000:7000;
ops.sdf_filter = 'PSP';
ops.freq = [2 6];

lfp_aligned = get_lfp_aligned(lfp,aligntime,ops);


%%
figuren('Renderer', 'painters', 'Position', [100 100 900 800]); hold on

nsubplot(2,2,[1 2],1)
for ch_i = 1:32
    if ismember(ch_i,1:16) ...
            color_line_value = 'r';
    else
        color_line_value = 'b';
    end

    plot(ops.timewin, nanmean(lfp_aligned.(['lfp_' int2str(ch_i)]))+(10*ch_i),'color',color_line_value)

end
set(gca,'YDir','Reverse')
hline(16.5*10,'k')
vline(0,'k')

for ch_i = 1:32
    electrode_mean(ch_i,:) = nanmean(lfp_aligned.(['lfp_' int2str(ch_i)]));
end

nsubplot(2,2,1,2)
plot(ops.timewin, nanmean(electrode_mean(1:16,:)),'color','r')
vline(0,'k')

nsubplot(2,2,2,2)
plot(ops.timewin, nanmean(electrode_mean(17:32,:)),'color','b')
vline(0,'k')



%% Align spikes

[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);







%% HOSD test

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

% Run sorting algorithm
% % [spikes,hos,znrm] = HOSD_spike_detection(data,params); % Run extraction
% % spike_cluster = sort_spikes(spikes); % Sort spikes

mat_savefile = fullfile(dirs.mat_data,['beans-test.mat']);
save(mat_savefile,'spikes','hos','znrm','spike_cluster','-v7.3')

% Processing HOSD spikes
% Plot summary sorting figure
fig = plot_clusters_2(spike_cluster);
set(fig,'renderer','painters','Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


% Restructure spike data for future analysis
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

% Align spikes on events
% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table events
event_table = get_event_table(ops);

% Align data to event
% Define parameters
ops.event_port = 2;
events = get_agl_t_trials(event_table, ops);
events = events(~strcmp(events.cond_label,'error'),:);
aligntime = events.stimulusOnset_ms;

ops.timewin = -1000:5000;
ops.sdf_filter = 'PSP';

[sdf, raster] = get_spikes_aligned(spikes_out,aligntime,ops);

% Produce spike figures for each neuron

% Get list of neurons
names = fieldnames( spikes_out.time );

% Define times for figures
xlim_vals = [-1000 5000];
ylim_vals = [0 60];

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
    % print(fig,['C:\KIKUCHI-LOCAL\script\kikuchi-data\data-extraction\doc\troy-agl_t-2021-11-05-HOSD-',['DSP' ch],'.pdf'],'-r400','-bestfit','-dpdf')
    % close all
end


