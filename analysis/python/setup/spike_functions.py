import numpy as np
from scipy.signal import convolve

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