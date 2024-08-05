import numpy as np

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def create_violation_alignment_event(event_table):
    # Ensure the DataFrame has the necessary columns
    if 'violation_ms' not in event_table.columns:
        event_table['violation_ms'] = np.nan

    for trial_i in range(len(event_table)):
        cond_value = event_table.at[trial_i, 'cond_value']
        stimulus_onset = event_table.at[trial_i, 'stimulusOnset_ms']

        if cond_value in {3, 7, 14, 1, 5, 13}:
            event_table.at[trial_i, 'violation_ms'] = stimulus_onset + 1127
        elif cond_value in {4, 8, 15, 2, 6, 16}:
            event_table.at[trial_i, 'violation_ms'] = stimulus_onset + 2253
        elif cond_value in {9, 10, 11, 12}:
            event_table.at[trial_i, 'violation_ms'] = np.nan

    return event_table
