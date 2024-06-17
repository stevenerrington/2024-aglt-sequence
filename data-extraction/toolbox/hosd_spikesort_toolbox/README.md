# HOSD Spike Sorting

The code in this repository implements higher-order spectral decomposition (HOSD) based spike sorting for Matlab.

For general details about the algorithm, see

C. K. Kovach and M. A. Howard III, [Decomposition of higher-order spectra for blind multiple-input deconvolution, pattern identification and separation](https://www.sciencedirect.com/science/article/pii/S0165168419302609), Signal Processing, 165 (2019), pp. 357–379.
 
A manuscript describing the application to spike sorting is in preparation. 

# How to use it
The full pipeline is contained in 3 scripts, to be run in sequence 
#### 1. Detect spikes 
        
```
    detected_spikes = HOSD_spike_detection(data)
    % where data is a structure with the following fields
    %    data.dat = data array as Samples x Channels
    %    data.fs  = Sample rate
```
#### 2. Sort spikes 


    sorted_spikes = sort_spikes(detected_spikes)

#### 3. Create summary plots and statistics for review

    plot_clusters(sorted_spikes)

## Additional notes
Default parameter values can be reviewed and adjusted in [HOSD_default_params.m](https://research-git.uiowa.edu/kovachc/hosd_spikesort/-/blob/master/HOSD_default_params.m).

To modify parameters outside of HOSD_default_params, do the following:
```
params.([some_field]) = [some_value]; % Set the value(s) you want to change
params = HOSD_default_params(params); %Set all other values to their defaults
detected_spikes = HOSD_spike_detection(data,params) % Run the algorithm
    ...
```

Performance metrics comparing output to ground truth, when availabe, can be computed with   

    perf = cluster_performance(sorted_spikes,ground_truth_spike_times)






