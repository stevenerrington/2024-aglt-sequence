import numpy as np
from scipy.signal import convolve


# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def spk_convolve(spk_data, session_end_time, conv_type):
    """
    Convolve spike data with a specified kernel to obtain a SessionSDF (Spike Density Function).

    Parameters:
    spk_data (array-like): Array of spike times (in samples).
    session_end_time (int): The length of the session (number of time points).
    conv_type (str): Type of convolution kernel to use ('PSP' or 'Gauss').

    Returns:
    numpy.ndarray: Convolved spike density function.
    """
    if conv_type == 'PSP':
        xrange = np.arange(101)  # 0 to 100
        R_non_normalized = (1 - np.exp(-xrange / 1)) * np.exp(-xrange / 100)
        R = R_non_normalized / np.sum(R_non_normalized)
        R2use = np.concatenate([np.zeros(len(xrange)), R])
    
    elif conv_type == 'Gauss':
        tg = 1  # Not used in this context
        td = 20  # From Pouget et al 2005 (in turn from Sayer et al 1990)
        normalize = 1  # Not used in this context
        
        mu = 0
        sd = td
        N = int(sd * 5)
        t = np.arange(-N, N + 1)
        R2use = (1 / (np.sqrt(2 * np.pi * sd**2))) * np.exp(-t**2 / (2 * sd**2))
    
    else:
        raise ValueError("Invalid ConvType. Use 'PSP' or 'Gauss'.")

    spk_data = np.array(spk_data)
    spk_data = spk_data[spk_data > 0]
    S2 = np.zeros(session_end_time)
    S2[spk_data] = 1

    # Perform convolution
    session_sdf = convolve(S2, R2use, mode='same') * 1000
    
    return session_sdf

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def spk_align(sdf_session, spk_times, align_times, time_win):
    """
    Aligns spike density functions and spike times based on alignment times and time window.

    Parameters:
    sdf_session (numpy.ndarray): The spike density function session data.
    spk_times (numpy.ndarray): The spike times.
    align_times (numpy.ndarray): Times to align the data.
    time_win (tuple): Time window for alignment (start, end).

    Returns:
    tuple: (sdf_aligned, raster_aligned)
    - sdf_aligned: Array of aligned spike density functions.
    - raster_aligned: List of arrays with aligned spike times.
    """
    # Initialize the output arrays
    sdf_aligned = np.full((len(align_times), len(time_win)), np.nan)
    raster_aligned = [None] * len(align_times)

    # Time window range
    time_win_start = time_win[0]
    time_win_end = time_win[-1]
    time_win_range = len(time_win)
    
    # Align the SDF and raster data
    for ii in range(len(align_times)):
        try:
            # Define the start and end indices for slicing
            start_idx = align_times[ii] + time_win_start
            end_idx = align_times[ii] + time_win_end
            
            # Extract the SDF segment
            sdf_aligned[ii, :min(time_win_range, len(sdf_session[start_idx:end_idx]))] = \
                sdf_session[start_idx:end_idx]
            
            # Extract and align spike times
            raster_aligned[ii] = spk_times[(spk_times > start_idx) & (spk_times < end_idx)] - align_times[ii]
        
        except IndexError:
            # Handle case where indices are out of bounds
            sdf_aligned[ii, :] = np.nan
            raster_aligned[ii] = []

    return sdf_aligned, raster_aligned

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def get_sound_aligned_sdf(sdf_aligned, raster_aligned, stimulusLog, event_table):
    '''
    The get_sound_aligned_sdf function aligns and extracts data for each sound presented in a trial, 
    creating a DataFrame that includes aligned SDF values, raster data, sound codes, sound positions, 
    and condition labels. This is useful for analyzing how neural responses are aligned with sound events 
    across different experimental conditions.
    '''
    
    # Define the onset times of sounds in milliseconds.
    sound_onset_ms = [0, 563, 1126, 1689, 2252]

    # Initialize a counter and create an empty DataFrame to store results.
    count = 0
    sound_sdf = pd.DataFrame(columns=['SDF_Value', 'Raster_Value', 'Sound_Code', 'Sound_Position', 'Cond_Label'])

    # Iterate over each trial in the event table.
    for trial_i in range(len(event_table)):
        # Process only trials that do not have 'error' as their condition label.
        if event_table['cond_label'][trial_i] != 'error':
            # For each sound onset time, align the SDF and raster data.
            for index, sound_i in enumerate(sound_onset_ms):
                count += 1
                
                # Extract the SDF values around the sound onset time.
                sdf_value = sdf_aligned[trial_i][1000 + np.arange(-500, 500, 1) + sound_i]
                
                # Compute the raster values relative to the sound onset time.
                raster_value = [int(r - sound_i) for r in raster_aligned[trial_i]]
                
                # Retrieve the sound code corresponding to the current trial's condition.
                sound_code = stimulusLog['sound_' + str(index+1) + '_code'][int(event_table['cond_value'][trial_i])-1]

                # Determine the position of the sound (e.g., 'position_0', 'position_1', etc.).
                sound_position = f'position_{index}'
                
                # Get the condition label for the trial.
                cond_label = event_table['cond_label'][trial_i]     
                    
                # Create a new DataFrame row with the extracted values.
                new_row = pd.DataFrame({
                    'SDF_Value': [sdf_value],
                    'Raster_Value': [raster_value],
                    'Sound_Code': [sound_code],
                    'Sound_Position': [sound_position],
                    'Cond_Label': [cond_label]
                })
        
                # Append the new row to the existing DataFrame.
                sound_sdf = pd.concat([sound_sdf, new_row], ignore_index=True)

    # Return the DataFrame containing aligned SDF and raster data.
    return sound_sdf

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------

def zscore_sdf(data, baseline_period):
    """
    Z-scores each trial in the data array using a specified baseline period.

    Parameters:
    - data: 2D numpy array of shape (trials, time).
    - baseline_period: Tuple specifying the start and end indices of the baseline period.

    Returns:
    - zscored_data: 2D numpy array of z-scored data with the same shape as input data.
    """
    start_idx, end_idx = baseline_period
    
    # Initialize the array for z-scored data
    zscored_data = np.zeros_like(data)
    
    # Iterate over each trial to compute the z-scored values
    for trial in range(data.shape[0]):
        # Extract the baseline period data for the current trial
        baseline_data = data[trial, start_idx:end_idx]
        
        # Compute the mean and standard deviation of the baseline period
        baseline_mean = np.nanmean(baseline_data)
        baseline_std = np.nanstd(baseline_data)
        
        # Z-score the entire trial based on the baseline period statistics
        zscored_data[trial] = (data[trial,:] - baseline_mean) / baseline_std
    
    return zscored_data

# ----------------------------------------------------------------------------------------------
#  /////////////////////////////////////////////////////////////////////////////////////////////
# ----------------------------------------------------------------------------------------------
def find_sig_modulation(data, nstds, nsamples, baseline_timewin, timewin):
    """
    Identifies periods in which the values of the array go above or below nstds (standard deviations)
    from the baseline for at least nsamples consecutive samples.

    Parameters:
    data (numpy array): Input array of shape (1, ntimes)
    nstds (float): Number of standard deviations from the baseline to consider
    nsamples (int): Minimum number of consecutive samples to consider as a period
    baseline_timewin (tuple): Start and end times (in seconds) for baseline period
    timewin (tuple): Start and end times (in seconds) for the whole data period

    Returns:
    List of tuples: Each tuple contains the start and end indices of a period
    """
    # Example usage
    # nstds = 2  # Number of standard deviations
    # nsamples = 50  # Minimum number of consecutive samples
    
    baseline_start_time, baseline_end_time = baseline_timewin
    
    baseline_indices = np.where((timewin >= baseline_start_time) & (timewin <= baseline_end_time))[0]

    # Calculate the baseline and standard deviation
    baseline = np.nanmean(data[baseline_indices])
    std_dev = np.nanstd(data[baseline_indices])
   
    # Calculate the upper and lower thresholds
    upper_threshold = baseline + (nstds * std_dev)
    lower_threshold = baseline - (nstds * std_dev)
    
    # Find the periods above or below the thresholds
    above_or_below_threshold = (data > upper_threshold) | (data < lower_threshold)
    
    # Find periods with at least nsamples consecutive True values
    mod_epochs = []
    start = None

    for i in range(len(above_or_below_threshold)):
        if above_or_below_threshold[i]:
            if start is None:
                start = i
        else:
            if start is not None:
                if i - start >= nsamples:
                    mod_epochs.append((start, i - 1))
                start = None

    # Check if the last period extends to the end of the array
    if start is not None and len(above_or_below_threshold) - start >= nsamples:
        mod_epochs.append((start, len(above_or_below_threshold) - 1))
    
    return mod_epochs
