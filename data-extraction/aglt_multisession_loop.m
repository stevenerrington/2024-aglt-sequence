%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- AGLt multi-session extraction script  -------------------------------
      S P Errington, 2024

This script was developed with the intention to loop through raw data files
listed in the ephysLog, determine what system is was recorded with, and
then extract the relevant signals and information. This should allow for a
common data format (spikes, spk info, event times, and lfps) to be
extracted for data recorded on both NeuraLynx and TDT.
///////////////////////////////////////////////////////////////////////////
%} 

%% Run initial extraction
% 

% Get session number
for session_i = 37
    % Clear workspace to reduce clutter
    clearvars -except C session_i dirs ephysLog

    % Detemine recording system
    exp_system = ephysLog.sys{session_i}; % Experiment type [agl, opto]

    % Provide user update
    h = waitbar(session_i/size(ephysLog,1),{['Processing session ' int2str(session_i) ' of ' int2str(size(ephysLog,1))], ephysLog.session{session_i}});

    % Run extraction functions
    switch exp_system
        case 'plex' % For Neuralynx
            kikuchi_nlx_extractFunction(ephysLog, session_i, dirs)
        case 'tdt' % For TDT
            kikuchi_tdt_extractFunction(ephysLog, session_i, dirs)
    end
    
    close(h)
end

%% Import spike sorted data and add to processed data file
% Once Kilosort has been run

for session_i = 1:33 % size(ephysLog,1)
    tic
    outfile_name = ephysLog.session{session_i}; % Processed file name
    exp_system = ephysLog.sys{session_i}; % Experiment type [agl, opto]

    % Provide user update
    %h = waitbar(session_i/size(ephysLog,1),{['Processing session ' int2str(session_i) ' of ' int2str(size(ephysLog,1))], ephysLog.session{session_i}});

    switch exp_system
        case 'plex' % For Neuralynx (32000Hz sampling)
            kikuchi_phy_import(outfile_name, dirs, 32000)
        case 'tdt' % For TDT  (24414.0625Hz sampling)
            kikuchi_phy_import(outfile_name, dirs, 24414.0625)
    end

    %close(h)
    toc
end
