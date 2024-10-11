# Import dependencies -------------------------------------------------------------
import pandas as pd
import sys


# Setup directory structure -------------------------------------------------------------
def set_directories():
    """
    Setup directories
    """

    if sys.platform == "win32":
        dirs = {
            "raw_data": r"C:\KIKUCHI-LOCAL\data\ephys\raw",
            "bin_data": r"C:\KIKUCHI-LOCAL\data\ephys\bin",
            "kilosort": r"C:\KIKUCHI-LOCAL\data\ephys\ks",
            "mat_data": r"C:\KIKUCHI-LOCAL\data\ephys\mat",
            "doc_data": r"C:\KIKUCHI-LOCAL\script\2024-aglt-sequence\data-extraction\doc",
            "py_data": r"C:\KIKUCHI-LOCAL\data\ephys\py_conv",
        }
    elif sys.platform == "darwin":
        dirs = {
            "raw_data": r"",  # Analysis not valid on mac
            "bin_data": r"",  # Analysis not valid on mac
            "kilosort": r"",  # Analysis not valid on mac
            "mat_data": r"/Users/stevenerrington/Library/CloudStorage/OneDrive-NewcastleUniversity/Research/agl_t/mat",
            "doc_data": r"/Users/stevenerrington/Desktop/Projects/2024-aglt-sequence/data-extraction/doc",
            "py_data": r"/Users/stevenerrington/Library/CloudStorage/OneDrive-NewcastleUniversity/Research/agl_t/py",
        }

    return dirs


# Import experiment details -------------------------------------------------------------
def import_exp_map():
    # Define the base URL and sheet names
    base_url = "https://docs.google.com/spreadsheets/d/{}/gviz/tq?tqx=out:csv&sheet={}"
    sheet_id = "1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg"

    # Read the ephysLog sheet
    ephysLog_url = base_url.format(sheet_id, "agl_t")
    ephysLog = pd.read_csv(ephysLog_url)

    # Read the stimulusLog sheet
    stimulusLog_url = base_url.format(sheet_id, "agl_t_stimuli")
    stimulusLog = pd.read_csv(stimulusLog_url)

    # Read the spike_log sheet
    spike_log_url = base_url.format(sheet_id, "agl_t_spikes")
    spike_log = pd.read_csv(spike_log_url)

    return ephysLog, stimulusLog, spike_log


# Clean and curate ephys log -------------------------------------------------------------
def clean_exp_map(ephysLog):
    """
    Clean the ephysLog DataFrame by removing the double session entry.

    Parameters:
    ephysLog (pd.DataFrame): The input DataFrame containing spike log data.

    Returns:
    pd.DataFrame: The cleaned DataFrame.
    """

    # Select every other row
    ephysLog_out = ephysLog.iloc[::2].reset_index(drop=True)

    return ephysLog_out


# Clean and curate spike log -------------------------------------------------------------
def clean_spike_map(spike_log):
    """
    Clean the spike_log DataFrame by removing neurons that did not meet criteria.

    Parameters:
    spike_log (pd.DataFrame): The input DataFrame containing spike log data.

    Returns:
    pd.DataFrame: The cleaned DataFrame.
    """
    # Remove neurons that did not meet criteria post-examination
    spike_log_out = spike_log[spike_log["useNeuron"] == 1]
    spike_log_out = spike_log_out[spike_log_out['area'] != 'hpc']

    return spike_log_out
