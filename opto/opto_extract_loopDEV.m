
% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
optoLog = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'opto'));

unique_days = unique(optoLog.date);

for day_i = 1:length(unique_days)
    day_sessions = [];
    day_sessions = find(optoLog.date == unique_days(day_i));
    clc; fprintf('Extracting broadband signal from day %i of %i... \n', day_i, length(unique_days))

    concat_spk_ncs = [];

    for session_n = 1:length(day_sessions)
        session_i = day_sessions(session_n);
        fprintf(' - merging session %i of %i... \n', session_n, length(day_sessions))
        session_tag = optoLog.file_n{session_i}; % Experimental file tag

        spk_ncs_out = [];

        for ch_n = 1:n_channels
            clear filepart_name NCSpath spk_ncs_in

            filepart_name = ['CSC' int2str(ch_n) session_tag];
            NCSpath = [fullfile(dirs.raw_data,exp_filename,filepart_name) '.ncs'];

            spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(dirs.raw_data,exp_filename));
            spk_ncs_out(ch_n,:) = spk_ncs_in.dat';
        end

        concat_spk_cut(session_n) = size(spk_ncs_out,2);
        concat_spk_ncs = [concat_spk_ncs, spk_ncs_out];


    end

    % Create a binary file and export the restructure broadband data
    clear bin_out_file
    bin_out_file = fopen([dirs.bin_data outfile_name(1:end-1) '.dat'],'wb');
    fwrite(bin_out_file,concat_spk_ncs,'int16');
    fclose(bin_out_file);
    mkdir(fullfile(dirs.kilosort,outfile_name(1:end-1)));

end