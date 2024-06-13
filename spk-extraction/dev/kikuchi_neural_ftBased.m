clear all; clc; warning off

%% Configuration & setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Admin --------------------------------------------------------------
monkey = 'troy'; % Monkey name [troy, chief]

% Experimental parameters -------------------------------------------
n_channels = 32; % Number of channels recorded in session

% Key setup variables
exp_filename = '2021-07-05_09-56-06_AGL2'; % Experimental raw data
task = 'agl'; % Experiment type [agl, opto]
session_n = '0004'; % Experimental file tag

% Define experimental/data directories -------------------------------
outfile_name = [monkey '-' task '-' exp_filename(1:10)]; % Processed file name
set_extract_dirs % Set experimental directories (i.e. data, scripts, etc...)


filelabels_spk = get_ncs_filelabel(fullfile(raw_dir,[exp_filename '\']), ['CSC1_' session_n '.ncs'],32);
data = ft_read_neuralynx_interp(filelabels_spk);

for ch_n = 1:n_channels
    clear filepart_name NCSpath spk_ncs_in

    filepart_name = ['CSC' int2str(ch_n) '_' session_n];
    NCSpath = [fullfile(raw_dir,exp_filename,filepart_name) '.ncs'];

    spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(raw_dir,exp_filename));
    spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
end

% Create a binary file and export the restructure broadband data
clear bin_out_file
bin_out_file = fopen([bin_dir outfile_name '.dat'],'wb');
fwrite(bin_out_file,spk_ncs_out,'int16');
fclose(bin_out_file);