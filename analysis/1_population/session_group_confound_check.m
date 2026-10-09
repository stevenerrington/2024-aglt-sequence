%% Session/group confound check: is sequence identity confounded with which session it came from?
%
% If certain sequence groups are disproportionately represented in
% particular sessions (different days, different neuron populations,
% different monkeys), then "decoding sequence group" from population
% activity could trivially reduce to "decoding which session this trial
% came from" -- since sessions differ in which neurons were recorded at
% all, with no requirement for any real-time computation. This would be a
% much more mundane explanation for the frontal 89% / auditory 67%
% identical-stimulus confound than genuine anticipatory coding, and it is
% a different failure mode from the trial-order check already run (that
% tested WITHIN-session drift; this tests BETWEEN-session imbalance).
%
% No neural data needed -- only each session's event table. Uses
% chi2_independence.m (already in analysis/0_functions) for a formal test
% of the session x group contingency table, plus per-monkey tallies since
% Monkey T and Monkey W were recorded on different rigs/days entirely.
%
% Assumes: dirs, spike_log already in the workspace.

seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};

unique_sessions = unique(spike_log.session);
n_sessions = length(unique_sessions);

session_group_counts = zeros(n_sessions, 4);
session_monkey = cell(n_sessions,1);

for session_i = 1:n_sessions

    try
        event_table_in = load(fullfile(dirs.mat_data, [unique_sessions{session_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);
    cond_value_session = event_table.cond_value(nonviol_mask);

    for g = 1:4
        session_group_counts(session_i,g) = sum(ismember(cond_value_session, seq_cond_values{g}));
    end

    monkey_rows = strcmp(spike_log.session, unique_sessions{session_i});
    monkey_here = unique(spike_log.monkey(monkey_rows));
    if ~isempty(monkey_here)
        session_monkey{session_i} = monkey_here{1};
    end

end

valid_sessions = sum(session_group_counts,2) > 0;
session_group_counts = session_group_counts(valid_sessions,:);
session_monkey = session_monkey(valid_sessions);
unique_sessions_valid = unique_sessions(valid_sessions);

fprintf('Sessions with usable trials: %d\n', size(session_group_counts,1));
fprintf('Per-group trial totals across all sessions: %s\n', mat2str(sum(session_group_counts,1)));

%% Formal test: session x group contingency table

[chi2_stat, df, p_value, expected] = chi2_independence(session_group_counts);
fprintf('\nSession x sequence-group chi-square test of independence:\n');
fprintf('chi2(%d) = %.2f, p = %.4g\n', df, chi2_stat, p_value);
if p_value < 0.05
    fprintf(['--> Sequence group is NOT evenly distributed across sessions.\n' ...
             '    This is a plausible confound: "decoding sequence group" could partly or wholly\n' ...
             '    reduce to "decoding which session/neuron-population this trial came from."\n']);
else
    fprintf(['--> No evidence sequence group is unevenly distributed across sessions.\n' ...
             '    Session-level imbalance does not obviously explain the identical-stimulus confound.\n']);
end

%% Same test restricted to each monkey separately (in case it only shows up within one animal)

for m = {'troy','walt'}
    m_idx = strcmp(session_monkey, m{1});
    if sum(m_idx) < 2
        continue
    end
    counts_m = session_group_counts(m_idx,:);
    [chi2_m, df_m, p_m, ~] = chi2_independence(counts_m);
    fprintf('\nMonkey %s only (%d sessions): chi2(%d) = %.2f, p = %.4g\n', m{1}, sum(m_idx), df_m, chi2_m, p_m);
end

%% Print the full table so you can eyeball it directly

fprintf('\nSession-by-session group counts (columns = group 1-4, i.e. cond_value pairs {1,5} {2,6} {3,7} {4,8}):\n');
for s = 1:length(unique_sessions_valid)
    fprintf('%-30s [%s]  monkey: %s\n', unique_sessions_valid{s}, ...
        num2str(session_group_counts(s,:)), session_monkey{s});
end
