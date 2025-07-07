
function kikuchi_phy_import(outfile_name, dirs, fs_in, n_channels)

ops = struct();
ops.rootZ = fullfile(dirs.kilosort,outfile_name);
ops.bin_file = [dirs.bin_data outfile_name '.dat'];
ops.nCh = n_channels;
ops.fs = fs_in;

if exist(fullfile(ops.rootZ,'params.py')) == 2
    try
        [spikes] = phy2mat(ops);
        [spikes.spk_info] = phyinfo2mat(ops);
        fprintf(['- phy import successful!: ' outfile_name ' \n'])
    catch
        fprintf(['- error importing phy curation: ' outfile_name ' \n'])
    end

    save(fullfile(dirs.mat_data,[outfile_name '.mat']),'-append','spikes');
else
    disp(['no kilosort output detected: ' outfile_name])
end


end