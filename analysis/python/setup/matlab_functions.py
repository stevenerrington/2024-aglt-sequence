import scipy.io
import pandas as pd
import h5py
import os

def load_mat_session(dirs, session_name):
    full_path = os.path.join(dirs['mat_data'], session_name + '.mat')

    with h5py.File(full_path, 'r') as mat_file:
        
        event_table =  mat_file['event_table']
        lfp =  mat_file['lfp']
        spikes =  mat_file['spikes']
        spk_info =  mat_file['spk_info']

        # Convert 'event_table' and 'spk_info' to pandas DataFrame
        # Assuming they are 2D arrays
        event_table = pd.DataFrame(event_table[:])
        spk_info = pd.DataFrame(spk_info[:])

        # Convert 'lfp' and 'spikes' to dictionaries
        # Assuming they are 1D arrays or can be converted to dictionaries
        lfp = dict(enumerate(lfp[:]))  # Convert to dict with indices as keys
        spikes = dict(enumerate(spikes[:]))  # Convert to dict with indices as keys


        return event_table, lfp, spikes, spk_info

