
%% 
%{ 
///////////////////////////////////////////////////////////////////////////
----- Kikuchi lab main script -----------------------------------------
      S P Errington, 2024
///////////////////////////////////////////////////////////////////////////
%} 

%% Workspace configuration and setup //////////////////////////////////////
% This series of commands and scripts must be ran prior to any other
% scripts, as they serve as dependencies.

% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
ephysLog_all = import_exp_map();
ephysLog = clean_exp_map(ephysLog_all);


spk_info_troy = [];

for session_i = 1:33 % size(ephysLog,1)
    outfile_name = ephysLog.session{session_i}; % Processed file name
    load(fullfile(dirs.mat_data,[outfile_name '.mat']),'spk_info');

    % Admin information
    n_units_session = size(spk_info,1);
    n_units_probe1 = sum(ismember(spk_info.site,[1:16]));
    n_units_probe2 = sum(ismember(spk_info.site,[17:32]));
    
    session = repmat(ephysLog.session(session_i),n_units_session,1);

    session_allidx = find(strcmp(ephysLog.session{session_i},ephysLog_all.session));
        
    area = [repmat(ephysLog_all.area_label(session_allidx(1)),n_units_probe1,1);...
        repmat(ephysLog_all.area_label(session_allidx(2)),n_units_probe2,1)];

    ml = repmat(ephysLog.nan_x(session_i),n_units_session,1);
    ap = repmat(ephysLog.nan_y(session_i),n_units_session,1);

    admin_table = table(session, area, ml, ap);
    spk_info_troy = [spk_info_troy; [admin_table, spk_info]];
end


%%