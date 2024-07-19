
%{ 
///////////////////////////////////////////////////////////////////////////
Plot neurophysiology signals for methods figure
- 16 channels of broadband signal for 5 seconds
- 1 channel of broadband signal with spikes notes
- Example waveform
///////////////////////////////////////////////////////////////////////////
%} 

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
session_i = 44;
monkey = ephysLog.monkey{session_i}; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = ephysLog.data_folder{session_i}; % Experimental raw data
task = ephysLog.task{session_i}; % Experiment type [agl, opto]
session_n = ephysLog.file_n{session_i}; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = ephysLog.session{session_i}; % Processed file name

dirs.raw_data = ephysLog.data_dir{session_i};

%% Load processed data
load(fullfile(dirs.mat_data,outfile_name))

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array

clear tdt_data 
tdtFun = @TDTbin2mat;
tdt_data = tdtFun(fullfile(dirs.raw_data,exp_filename));

%% Produce figures
close all
plot_time_ms = [200000:210000];
plot_time_samples = round((plot_time_ms/1000)*tdt_data.streams.Raw1.fs);

figuren('Renderer', 'painters', 'Position', [100 100 300 500]); hold on;
for ch_i = 1:16
    data_in = [];
    data_in = lfp_filter(tdt_data.streams.Raw1.data(ch_i,plot_time_samples)*10000,600, 8000, round(tdt_data.streams.Raws.fs));
    plot(data_in + 1.8*(ch_i-1),'k','LineWidth',0.1);

end
set(gca,'YDir','Reverse','XColor',[ 1 1 1 ],'YColor',[ 1 1 1 ])

%% 
ch_i = 16;
figuren('Renderer', 'painters', 'Position', [100 100 300 100]); hold on;
data_in = lfp_filter(tdt_data.streams.Raw1.data(ch_i,plot_time_samples)*15000, 600, 8000, round(tdt_data.streams.Raw1.fs));

plot(plot_time_ms/1000,data_in,'color',[0.5 0.5 0.5],'LineWidth',0.1);
vline(spikes.time.DSP16a(spikes.time.DSP16a > plot_time_ms(1) & spikes.time.DSP16a < plot_time_ms(end) )/1000,'r-')

xlim([200 210])