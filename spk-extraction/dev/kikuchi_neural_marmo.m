clear all; warning off

% Configuration -----------------------------------------------------
% Admin
monkey = 'beans';
task = 'vocal';

% File naming
raw_filename = '2023-11-03_14-32-27';
outfile_name = [monkey '-' task '-' raw_filename(1:10)];

% Define experiment parameters
session_n = 1;
n_channels = 32;

% Setup directories
raw_dir  = 'C:\KIKUCHI-LOCAL\data\ephys\raw\';
bin_dir  = 'C:\KIKUCHI-LOCAL\data\ephys\bin\';
ks_dir   = 'C:\KIKUCHI-LOCAL\data\ephys\ks\';
proc_dir = 'C:\KIKUCHI-LOCAL\data\ephys\mat';

% Behavioral data -------------------------------------------------------
% Read in events
evnt = readnev(['Events.nev'],fullfile(raw_dir,raw_filename));

% Map events


% Spike data -------------------------------------------------------
% Loop through all recorded channels, get the ncs, and process it.
% - Restructure data as nCh x nSample array
for ch_n = 1:n_channels
    clear filepart_name NCSpath spk_ncs_in 
    filepart_name = ['CSC' int2str(ch_n)];
    NCSpath = [fullfile(raw_dir,raw_filename,filepart_name) '.ncs'];
    fprintf('Extracting data from file %s ... \n', filepart_name);

    spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(raw_dir,raw_filename));
    spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
end

% Create a binary file and export the restructure broadband data
clear bin_out_file
bin_out_file = fopen([bin_dir outfile_name '.dat'],'wb');
fwrite(bin_out_file,spk_ncs_out,'int16');
fclose(bin_out_file);

% Run kilosort
% - Run in python: to be integrated here.
% conda activate kilosort
% python -m kilosort


% Import phy curated data
ops = struct();
ops.rootZ = fullfile(ks_dir,outfile_name);
ops.bin_file = [bin_dir outfile_name '.dat'];
ops.nCh = n_channels;
ops.fs = 32000;

[spikes] = phy2mat(ops);
[spk_info] = phyinfo2mat(ops);



% Local field potential data -------------------------------------------------------
for ch_n = 1:n_channels
    clear filepart_name NCSpath lfp_ncs_in 

    filepart_name = ['LFP' int2str(ch_n) '_000' int2str(session_n)];
    NCSpath = [fullfile(raw_dir,raw_filename,filepart_name) '.ncs'];

    lfp_ncs_in = readncs([filepart_name '.ncs'],fullfile(raw_dir,raw_filename));
    lfp_ncs_out(ch_n,:) = lfp_ncs_in.dat';
end





