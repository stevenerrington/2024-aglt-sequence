
n_trls = length(dat.T.start);
eye_win = [0.45 0.55];
clear eye_x eye_y time
for trial_i = 1:n_trls
    raw_timepoint = round(dat.eye.dat{trial_i}.timepoint*1000);
    raw_eyex = dat.eye.dat{trial_i}.gx;
    raw_eyey = dat.eye.dat{trial_i}.gy;

    for timestamp = 1:6000
        eye_x{trial_i}(timestamp,1) = raw_eyex(find(raw_timepoint <= timestamp,1,'last'));
        eye_y{trial_i}(timestamp,2) = raw_eyey(find(raw_timepoint <= timestamp,1,'last'));
        time{trial_i}(timestamp,3) = timestamp;

        outwin_flag(trial_i,timestamp) = eye_x{trial_i}(timestamp,1) > eye_win(2) | eye_x{trial_i}(timestamp,1) < eye_win(1);

    end
end


%% 

viol_trials = []; viol_trials = find(strcmp(event_table.cond_label,'viol'));
nonviol_trials = []; nonviol_trials = find(strcmp(event_table.cond_label,'nonviol'));

outwin_flag(viol_trials,:)