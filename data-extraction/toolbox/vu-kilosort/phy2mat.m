function [spikes] = phy2mat(ops)
%% Import spike information from Phy
ops.rootZ = ops.rootZ;
sp = loadKSdir(ops.rootZ);
[spikeTimes, spikeAmps, ~, spikeSites, spikeTemplates] = ksDriftmap(ops.rootZ);

%% Waveform extraction
%  Get Parameters
gwfparams.dataDir = ops.rootZ;           % KiloSort/Phy output folder
gwfparams.fileName = ops.bin_file;       % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    spikeTimes;    % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu;        % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

%  Run main extraction
wf = getWaveForms(gwfparams);


%% Setup variable space for spike and waveform data
% Note: I've checked and it seems as though noise clusters are dropped from
% the import, and we don't need to include/exclude them manually.

unitList = double(unique(sp.clu));
nUnits = length(unitList);
% Find site for each ID'd cluster
for unitIdx = 1:nUnits
    unit = unitList(unitIdx);
    cluster(unitIdx,1) = unit;
    site(unitIdx,1) = double(mode(spikeSites(sp.clu == unit)));
end

% Get labels for the output (e.g. DSP01a, WAV01a = first unit on ch 1)
for unitIdx = 1:nUnits
    unit_site = find(unitIdx == find(site == site(unitIdx,1)));
    dspString = num2str(site(unitIdx,1),'DSP%02i');
    wavString = num2str(site(unitIdx,1),'WAV%02i');
    clustLetter = char(unit_site+96); % 97 is the char code for 'a"
    unitDSP{unitIdx,1} = [dspString clustLetter];
    unitWAV{unitIdx,1} = [wavString clustLetter];
end

spkTable = table(cluster,site,unitDSP,unitWAV);
spkTable = sortrows(spkTable,'site');
%% Get spike information for each unit
% NOTE: Sampling rate adjustment for spike time is made here. After this,
% spikes are in ms, and not in samples.

for unitIdx = 1:nUnits
    unit = spkTable.cluster(unitIdx);
    spikes.time.(spkTable.unitDSP{unitIdx}) = (spikeTimes(sp.clu == unit)/ops.fs)*1000;
    spikes.amplitudes.(spkTable.unitDSP{unitIdx}) = spikeAmps(sp.clu == unit);
end


%% Get spike waveform for each unit
% NOTE: Due to array size issues, we have subsampled 2000 spikes.
% NOTE: Sampling rate adjustment for spike time is made here. After this,
% spikes are in ms, and not in samples.

for unitIdx = 1:nUnits
    unit = spkTable.cluster(unitIdx);
    spikes.waveform.(spkTable.unitWAV{unitIdx}) =...
        squeeze(wf.waveForms(wf.unitIDs == unit,:,spkTable.site(unitIdx),:));
    
    spikes.waveform_spkTime.(spkTable.unitWAV{unitIdx}) =...
        (wf.spikeTimeKeeps(wf.unitIDs == unit,:)/ops.fs)*1000;
end
