
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

%% Plot stimuli
sound_dir = 'S:\DATA\MARMO\ACUTE\Beans\Stimuli\marmoMvoc\BEST_renamed';


wavFiles = dir(fullfile(sound_dir,'*.wav'));
wavFiles = {wavFiles.name};


n_stimuli = size(wavFiles,2);
figuren('Renderer', 'painters', 'Position', [100 100 400 1200]); hold on;;

threshold = 0.05;

for file_i = 1:n_stimuli

    audio_file = wavFiles{file_i};

    [sounddata,fs_check] = audioread...
        ([sound_dir ,'\' audio_file]);

    sounddata = sounddata(:,1);

    time = (1:length(sounddata))/fs_check;

    clip_duration(file_i,1) = max(time);

    subplot(n_stimuli,1,file_i)
    plot(time, sounddata(:,1))
    xlim([0 5])
    set(gca,'XColor',[1 1 1], 'YColor', [1 1 1])
end



%% Adjust event table and align
% Behavioral data -------------------------------------------------------
% Read in events
ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n; 
clear event_table agl_t_event
event_table = get_event_table(ops);

% Event extraction
clear agl_t_event aligntime
ops.event_port = 2;
agl_t_event = get_agl_t_trials_nlx(event_table, ops);

wav_offset = [314 219 231 1367 1485 1434 784 884 672 617 536 596 1587 0 0 774 701 859];


for trial_i = 1:size(agl_t_event,1)
    agl_t_event.stimulusOnset_ms(trial_i) = agl_t_event.stimulusOnset_ms(trial_i) + wav_offset(agl_t_event.cond_value(trial_i));

    switch agl_t_event.cond_value(trial_i)
        case { 1 2 3 }
            agl_t_event.cond_label{trial_i} = 'nonpreferred';
        case { 4 5 6 }
            agl_t_event.cond_label{trial_i} = 'nonpreferred';
        case { 7 8 9 }
            agl_t_event.cond_label{trial_i} = 'preferred';
        case { 10 11 12 }
            agl_t_event.cond_label{trial_i} = 'nonpreferred';
        case { 13 14 15 }
            agl_t_event.cond_label{trial_i} = 'nonpreferred';
        case { 16 17 18 }
            agl_t_event.cond_label{trial_i} = 'nonpreferred';
    end
end


aligntime = agl_t_event.stimulusOnset_ms;


%% Plot waveforms
waves = fieldnames( spikes.waveform );

figuren('Renderer', 'painters', 'Position', [100 100 900 400]); hold on;
for wave_i = 1:length(waves)
    subplot(3,5,wave_i)
    if spk_info.site(wave_i) < 17
        color = 'r';
    else
        color = 'b';
    end

    plot(nanmean(spikes.waveform.(waves{wave_i})),'color',color,'LineWidth',1.5)
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1])
end

%% Plot LFP
% Setup LFP data
ops.timewin = [-1000:5000];
ops.freq = [1 30];
lfp = patch_fault_ch(lfp,4);
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);

trial_average_lfp = nanmean(lfp_aligned,3);

figuren('Renderer', 'painters', 'Position', [100 100 600 500]);
channels = 1:32;
fig_plot_exp = 20;

subplot(1,2,1); hold on
for ch_i = 1:16
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+fig_plot_exp*(ch_i-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-15 (length(1:16)*fig_plot_exp)]); yticks([fig_plot_exp*([1:16]-1)]); yticklabels(num2cell([1:16]))
xlim([-250 1000])
vline(0,'k')

subplot(1,2,2); hold on
for ch_i = 17:32
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+fig_plot_exp*(ch_i-16-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-15 (length(1:16)*fig_plot_exp)]); yticks([fig_plot_exp*([1:16]-1)]); yticklabels(num2cell([1:16]))
xlim([-250 1000])
vline(0,'k')


%% Plot theta LFP

ops.timewin = [-1000:5000];
ops.freq = [4 9];
lfp = patch_fault_ch(lfp,4);
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);

trial_average_lfp = nanmean(lfp_aligned,3);

figuren('Renderer', 'painters', 'Position', [100 100 600 500]);
channels = 1:32;
fig_plot_exp = 5;

subplot(1,2,1); hold on
for ch_i = 1:16
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+fig_plot_exp*(ch_i-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-15 (length(1:16)*fig_plot_exp)]); yticks([fig_plot_exp*([1:16]-1)]); yticklabels(num2cell([1:16]))
xlim([-250 1000])
vline(0,'k')

subplot(1,2,2); hold on
for ch_i = 17:32
    if ismember(ch_i,[1:16])
        color_line = [204 0 102]./255;
    else
        color_line = [0 153 153]./255;
    end

    plot(ops.timewin,  trial_average_lfp(ch_i,:)+fig_plot_exp*(ch_i-16-1),'color',color_line)
end
set(gca,'ydir', 'reverse')
ylim([-15 (length(1:16)*fig_plot_exp)]); yticks([fig_plot_exp*([1:16]-1)]); yticklabels(num2cell([1:16]))
xlim([-250 1000])
vline(0,'k')

%% Plot event-related LFP

ops.timewin = [-1000:5000];
ops.freq = [5 12];
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);

trial_average_lfp_amy_pref = nanmean(lfp_aligned(1:16,:,strcmp(agl_t_event.cond_label,'preferred')),3);
trial_average_lfp_amy_nonpref = nanmean(lfp_aligned(1:16,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);

trial_average_lfp_mpfc_pref = nanmean(lfp_aligned(17:32,:,strcmp(agl_t_event.cond_label,'preferred')),3);
trial_average_lfp_mpfc_nonpref = nanmean(lfp_aligned(17:32,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);


xlim_vals = [-250 1000];
ylim_vals = [-20 20];
lfp_figure_idv(1,1)=gramm('x',ops.timewin,'y',[trial_average_lfp_amy_pref; trial_average_lfp_amy_nonpref],'color',[repmat({'preferred'},16,1);repmat({'nonpreferred'},16,1)]);
lfp_figure_idv(1,1).stat_summary();
lfp_figure_idv(1,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);

lfp_figure_idv(2,1)=gramm('x',ops.timewin,'y',[trial_average_lfp_mpfc_pref; trial_average_lfp_mpfc_nonpref],'color',[repmat({'preferred'},16,1);repmat({'nonpreferred'},16,1)]);
lfp_figure_idv(2,1).stat_summary();
lfp_figure_idv(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);


fig = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
lfp_figure_idv.draw();

%% Plot spike activity
ops.timewin = -1000:5000;
ops.sdf_filter = 'Gauss';
[sdf, raster] = get_spikes_aligned(spikes,aligntime,ops);

% Get list of neurons
names = fieldnames( spikes.time );

% Define times for figures
xlim_vals = [-250 600];
ylim_vals = [0 30];
stimuli_id = 2;

for ch_i = 10
    clear spk_figure_idv fig

    % Get channel label
    ch = names{ch_i}(4:end);

    % Plot raster
    spk_figure_idv(1,1)=gramm('x',raster.(['DSP' ch]),'color',agl_t_event.cond_label);
    spk_figure_idv(1,1).geom_raster('geom',{'line'});
    spk_figure_idv(1,1).axe_property('XLim',xlim_vals);

    % Plot SDF
    spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',agl_t_event.cond_label);
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

%% Power spectrum

% ops.timewin = [-1000:5000];
% ops.freq = [4 30];
% [~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);
% trial_average_lfp_preferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')),3);
% trial_average_lfp_nonpreferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);

ops.freq = [4 30];
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);
%trial_average_lfp_preferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')),3);
%trial_average_lfp_nonpreferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);

lfp_in_pref = []; lfp_in_pref = squeeze(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')));
lfp_in_nonpref = []; lfp_in_nonpref = squeeze(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')));

clear trial_average_lfp*
trial_average_lfp_preferred = reshape(lfp_in_pref,[32, size(lfp_in_pref,2)*size(lfp_in_pref,3)]);
trial_average_lfp_nonpreferred = reshape(lfp_in_nonpref,[32, size(lfp_in_nonpref,2)*size(lfp_in_nonpref,3)]);



% Parameters for pwelch
window = 500; % Length of each segment
noverlap = 250; % Number of overlapping samples
nfft = 5000; % Number of FFT points

clear power_preferred f
for channel_i = 1:32
    [power_preferred(channel_i,:), f(channel_i,:)] = pwelch(trial_average_lfp_preferred(channel_i,:), window, noverlap, nfft, 1000, 'power');
    [power_nonpreferred(channel_i,:), f(channel_i,:)] = pwelch(trial_average_lfp_nonpreferred(channel_i,:), window, noverlap, nfft, 1000, 'power');

    if channel_i < 17
        channel_label{channel_i,1} = 'aHPC';
    else
        channel_label{channel_i,1} = 'mPFC';
   end
end


power_spectrum_figure(1,1)=gramm('x',f(1,:),'y',[power_preferred;power_nonpreferred],'color',[repmat({'preferred'},32,1); repmat({'nonpreferred'},32,1)],...
    'column',repmat(channel_label,2,1));
power_spectrum_figure(1,1).stat_summary();
power_spectrum_figure(1,1).axe_property('XLim',[0 30],'YLim',[0 250]);
fig = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
power_spectrum_figure.draw();

%% Coherence
ops.freq = [4 30];
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);
%trial_average_lfp_preferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')),3);
%trial_average_lfp_nonpreferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);


lfp_in_pref = []; lfp_in_pref = squeeze(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')));
lfp_in_nonpref = []; lfp_in_nonpref = squeeze(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')));

clear trial_average_lfp*
trial_average_lfp_preferred = reshape(lfp_in_pref,[32, size(lfp_in_pref,2)*size(lfp_in_pref,3)]);
trial_average_lfp_nonpreferred = reshape(lfp_in_nonpref,[32, size(lfp_in_nonpref,2)*size(lfp_in_nonpref,3)]);

clear Cxy f
for pair_i = 1:16
    % Compute the magnitude-squared coherence
    [coherence_pref(pair_i,:), f(pair_i,:)] = mscohere(trial_average_lfp_preferred(pair_i,:), trial_average_lfp_preferred(pair_i+16,:), window, noverlap, nfft, 1000);
    [coherence_nonpref(pair_i,:), f(pair_i,:)] = mscohere(trial_average_lfp_nonpreferred(pair_i,:), trial_average_lfp_nonpreferred(pair_i+16,:), window, noverlap, nfft, 1000);

end

clear coherence_figure
coherence_figure(1,1)=gramm('x',f(1,:),'y',[coherence_pref;coherence_nonpref],'color',[repmat({'preferred'},16,1); repmat({'nonpreferred'},16,1)]);
coherence_figure(1,1).stat_summary();
coherence_figure(1,1).axe_property('XLim',[0 30],'YLim',[0 0.8]);
fig = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
coherence_figure.draw();

%% Cross-correlation
ops.freq = [4 30];
[~, lfp_aligned] = get_lfp_aligned(lfp,aligntime,ops);
trial_average_lfp_preferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'preferred')),3);
trial_average_lfp_nonpreferred = nanmean(lfp_aligned(:,:,strcmp(agl_t_event.cond_label,'nonpreferred')),3);

clear x_corr_* lags
for pair_i = 1:16
    % Compute the crosscorrelation between signals
    [x_corr_pref(pair_i,:),lags(pair_i,:)] = xcorr(trial_average_lfp_preferred(pair_i,:), trial_average_lfp_preferred(pair_i+16,:),'normalized');
    [x_corr_nonpref(pair_i,:),~] = xcorr(trial_average_lfp_nonpreferred(pair_i,:), trial_average_lfp_nonpreferred(pair_i+16,:),'normalized');

end


clear xcorr_figure
xcorr_figure(1,1)=gramm('x',lags(1,:),'y',[x_corr_pref;x_corr_nonpref],'color',[repmat({'preferred'},16,1); repmat({'nonpreferred'},16,1)]);
xcorr_figure(1,1).stat_summary();
xcorr_figure(1,1).axe_property('XLim',[-300 300],'YLim',[-0.6 0.6]);
fig = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
xcorr_figure.draw();

%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% %
% %
% % %% HOSD test
% %
% % % Define algorithm parameters
% % params.lowpass = 4000; % lowpass cutoff ( Hz )
% % params.highpass = 200; % highpass cutoff ( Hz )
% % params.limit_memory = 5e8;
% % params.check_features = false;
% % params.randseed = 1;
% % saveAllRaw = 0; % boolean to keep spikes data (rather than just the clustered data)
% % saveFigFile = 0;
% %
% % % Setup data for function
% % clear data
% % data.dat = spk_ncs_out;
% % data.fs = 32000;
% %
% % % Run sorting algorithm
% % % % [spikes,hos,znrm] = HOSD_spike_detection(data,params); % Run extraction
% % % % spike_cluster = sort_spikes(spikes); % Sort spikes
% %
% % mat_savefile = fullfile(dirs.mat_data,['beans-test.mat']);
% % save(mat_savefile,'spikes','hos','znrm','spike_cluster','-v7.3')
% %
% % % Processing HOSD spikes
% % % Plot summary sorting figure
% % fig = plot_clusters_2(spike_cluster);
% % set(fig,'renderer','painters','Units','Inches');
% % pos = get(fig,'Position');
% % set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% %
% %
% % % Restructure spike data for future analysis
% % n_clusters = spike_cluster.Nclust; % Find number of spikes
% %
% % clear spikes_out
% % for cluster_i = 1:n_clusters
% %
% %     % Get spike times
% %     spikes_out.time.(['DSP' int2str(cluster_i)]) = ...
% %         spike_cluster.spike_times(find(spike_cluster.cl == cluster_i))*1000;
% %
% %     % Get spike waveforms
% %      spikes_out.waveform.(['WAV' int2str(cluster_i)]) = ...
% %         spike_cluster.avg_waves(:,cluster_i)';
% % end
% %
% % % Align spikes on events
% % % Behavioral data -------------------------------------------------------
% % % Read in events
% % ops.dirs.raw_data = dirs.raw_data; ops.filename = exp_filename; ops.session_n = session_n;
% % clear event_table events
% % event_table = get_event_table(ops);
% %
% % % Align data to event
% % % Define parameters
% % ops.event_port = 2;
% % events = get_agl_t_trials(event_table, ops);
% % events = events(~strcmp(events.cond_label,'error'),:);
% % aligntime = events.stimulusOnset_ms;
% %
% % ops.timewin = -1000:5000;
% % ops.sdf_filter = 'PSP';
% %
% % [sdf, raster] = get_spikes_aligned(spikes_out,aligntime,ops);
% % 
% % % Produce spike figures for each neuron
% % 
% % % Get list of neurons
% % names = fieldnames( spikes_out.time );
% % 
% % % Define times for figures
% % xlim_vals = [-1000 5000];
% % ylim_vals = [0 60];
% % 
% % for ch_i = 1:length(names)
% %     clear spk_figure_idv fig
% % 
% %     % Get channel label
% %     ch = names{ch_i}(4:end);
% % 
% %     % Plot raster
% %     spk_figure_idv(1,1)=gramm('x',raster.(['DSP' ch]),'color',events.cond_label);
% %     spk_figure_idv(1,1).geom_raster('geom',{'line'});
% %     spk_figure_idv(1,1).axe_property('XLim',xlim_vals,'YLim',[0 size(raster.(['DSP' ch]),1)]);
% % 
% %     % Plot SDF
% %     spk_figure_idv(2,1)=gramm('x',ops.timewin,'y',sdf.(['DSP' ch]),'color',events.cond_label);
% %     spk_figure_idv(2,1).stat_summary();
% %     spk_figure_idv(2,1).axe_property('XLim',xlim_vals,'YLim',ylim_vals);
% % 
% %     % Plot waveform average
% %     spk_figure_idv(3,1)=gramm('x',[1:length(spikes_out.waveform.(['WAV' ch]))],'y',spikes_out.waveform.(['WAV' ch]));
% %     spk_figure_idv(3,1).stat_summary();
% %     spk_figure_idv(3,1).axe_property('XLim',[-41 40],'YLim',[-35 35]);
% % 
% %     % Figure properties
% %     spk_figure_idv(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
% %     spk_figure_idv(3,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
% %     spk_figure_idv(1,1).geom_vline('xintercept',0,'style','k-');
% %     spk_figure_idv(2,1).geom_vline('xintercept',0,'style','k-');
% %     spk_figure_idv(1,1).set_names('y','Trials');
% %     spk_figure_idv(2,1).set_names('x','Time from event (ms)','y','FR (spk/sec)');
% % 
% %     % Figure layout
% %     spk_figure_idv(1,1).set_layout_options...
% %         ('Position',[0.1 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
% %         'legend',false,...
% %         'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
% %         'margin_width',[0.0 0.00],...
% %         'redraw',false);
% % 
% %     spk_figure_idv(2,1).set_layout_options...
% %         ('Position',[0.1 0.1 0.8 0.6],... %Set the position in the figure (as in standard 'Position' axe property)
% %         'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
% %         'margin_width',[0.0 0.00],...
% %         'redraw',false);
% % 
% %     spk_figure_idv(3,1).set_layout_options...
% %         ('Position',[0.7 0.5 0.2 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
% %         'legend',false,...
% %         'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
% %         'margin_width',[0.0 0.00],...
% %         'redraw',false);
% % 
% %     spk_figure_idv.set_title(['DSP' ch]);
% % 
% %     % Draw figure
% %     fig = figure('Renderer', 'painters', 'Position', [100 100 700 600]);
% %     spk_figure_idv.draw();
% % 
% %     % Save figure
% %     set(fig,'renderer','painters','Units','Inches');
% %     pos = get(fig,'Position');
% %     set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% %     % print(fig,['C:\KIKUCHI-LOCAL\script\kikuchi-data\data-extraction\doc\troy-agl_t-2021-11-05-HOSD-',['DSP' ch],'.pdf'],'-r400','-bestfit','-dpdf')
% %     % close all
% % end
% % 
% % 
