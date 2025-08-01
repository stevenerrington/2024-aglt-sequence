% Clear environment
clear all; clc; warning off

% Setup data directories for use throughout scripts
dirs = set_directories();

% Import and curate experimental log
[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = clean_exp_map(ephysLog);
spike_log = clean_spike_map(spike_log);


% Load session data
session_i = 1;
datafile = ephysLog.session{session_i};
load(fullfile(dirs.mat_data,datafile),'event_table')

% Initialize transition probability matrix with learned prior probabilities
% Assuming we have a prior transition probability matrix based on compatible sequences
% Example: Replace with actual prior probabilities if available
sounds = {'Start', 'YAG', 'KEM', 'LEK', 'PAV', 'RAP', 'X'};  % Define your sounds

for seq_i = 1:16
    observedSequenceList{seq_i}{1,1} = 'Start';
    for ele_i = 2:6
    observedSequenceList{seq_i}{1,ele_i} =...
        stimulusLog.(['sound_' int2str(ele_i-1) '_label']){seq_i};
    end
    observedSequenceList{seq_i}{1,7} = 'X';
end

for trial_i = 1:size(event_table,1)
    trial_seq = event_table.cond_value(trial_i);
    observedSequence_trial{trial_i,1} = observedSequenceList{trial_seq};
end

alpha = 1;  % Decay factor for adaptive updating, set between 0 and 1

learned_transitionProb = transitional_prob_array;  % Prior distribution from AGLe

learned_transitionProb_bck = learned_transitionProb./sum(learned_transitionProb,1);  % Prior distribution from AGLe
learned_transitionProb_fwd = learned_transitionProb./sum(learned_transitionProb,2);  % Prior distribution from AGLe

bayesian_transitionProb = learned_transitionProb;
bayesian_transitionProb_bck = [];
bayesian_transitionProb_fwd = [];

delta_bck_prob = [];
delta_fwd_prob = [];


% Iterate over each trial sequence
for trial = 1:length(observedSequence_trial)
    sequence = observedSequence_trial{trial};
    
    % Update transition probabilities based on observed transitions in sequence
    for t = 2:length(sequence)
        % Get current and next sound in the sequence
        prevElement = find(strcmp(sounds, sequence{t-1}));
        nextElement = find(strcmp(sounds, sequence{t}));
        
        % Calculate observed transition probability for this transition
        observedProb = 1;  % Since we observed this transition, probability is 1

        bayesian_transitionProb(prevElement, nextElement) = ...
            bayesian_transitionProb(prevElement, nextElement) + observedProb;
        
        % Normalize row to ensure probabilities sum to 1
        bayesian_transitionProb_bck(:,:,trial) = bayesian_transitionProb./sum(bayesian_transitionProb,1);
        bayesian_transitionProb_fwd(:,:,trial) = bayesian_transitionProb./sum(bayesian_transitionProb,2);

        if trial == 1
            delta_bck_prob(trial, t-1) = bayesian_transitionProb_bck(prevElement,nextElement,trial) - learned_transitionProb_bck(prevElement, nextElement);
            delta_fwd_prob(trial, t-1) = bayesian_transitionProb_fwd(prevElement,nextElement,trial) - learned_transitionProb_fwd(prevElement, nextElement);
        else
            delta_bck_prob(trial, t-1) = bayesian_transitionProb_bck(prevElement,nextElement,trial) - bayesian_transitionProb_bck(prevElement,nextElement,trial-1);
            delta_fwd_prob(trial, t-1) = bayesian_transitionProb_fwd(prevElement,nextElement,trial) - bayesian_transitionProb_fwd(prevElement,nextElement,trial-1);
        end

    end
end



violation_trial = 1;
event_table.cond_value(violation_trial)

figuren; hold on
plot(1:6,delta_bck_prob(1,:)*100)
plot(1:6,delta_bck_prob(115,:)*100)
xticks(1:6); xticklabels(observedSequence_trial{115}(2:end))
vline(5)


figuren; hold on
plot(smooth(delta_bck_prob(:,1)*100,20))
plot(smooth(delta_bck_prob(:,2)*100,20))
plot(smooth(delta_bck_prob(:,3)*100,20))
plot(smooth(delta_bck_prob(:,4)*100,20))
plot(smooth(delta_bck_prob(:,5)*100,20))
plot(smooth(delta_bck_prob(:,6)*100,20))
legend({'1','2','3','4','5'})