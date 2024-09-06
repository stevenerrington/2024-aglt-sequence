import numpy as np
import scipy

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

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def adjust_audio_onset_ms(session_i, event_table):
    # Load the .mat file
    session_audio_latency = scipy.io.loadmat(r'data-extraction\doc\session_audio_latency.mat')
    session_audio_latency = session_audio_latency['session_audio_latency']

    for trial_i in range(len(event_table)):
        event_table.at[trial_i, 'stimulusOnset_ms'] = event_table.at[trial_i, 'stimulusOnset_ms'] + session_audio_latency[session_i][0][trial_i][0]

    return event_table