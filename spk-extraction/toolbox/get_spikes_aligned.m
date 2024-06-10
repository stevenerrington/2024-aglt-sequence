function [sdf, raster] = get_spikes_aligned(spikes,alignTimes,ops)

timeWin = ops.timewin;


% Get spike times
names = fieldnames( spikes.time );
subStr = 'DSP';
DSPstruct = rmfield( spikes.time, names( find( cellfun( @isempty, strfind( names , subStr ) ) ) ) );
DSPnames = fieldnames(DSPstruct);

for ch_i = 1:size(DSPnames,1)
    DSPidx = ch_i;
    DSPlabel = DSPnames{DSPidx};
    spkTimes = []; spkTimes = round(spikes.time.(DSPlabel));

    sdf_session = [];
    sdf_session = SpkConvolver (spkTimes, round(max(alignTimes)+5000), ops.sdf_filter);
    sdf_aligned = nan(length(alignTimes),range(timeWin)+1);

    for ii = 1:length(alignTimes)
        try
            sdf_aligned(ii,:) = sdf_session(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));
            raster_aligned{ii,1} = spkTimes(spkTimes>alignTimes(ii)+timeWin(1) & spkTimes<alignTimes(ii)+timeWin(end))-alignTimes(ii);
        catch
            sdf_aligned(ii,:) = nan(1,length(timeWin));
            raster_aligned{ii,1} = [];
        end
    end

    sdf.(DSPlabel) = sdf_aligned;
    raster.(DSPlabel) = raster_aligned;

end