
spkTimes = spikes.time;
timeWin = [-1000:2000];

names = fieldnames( spkTimes );
subStr = 'DSP';
DSPstruct = rmfield( spkTimes, names( find( cellfun( @isempty, strfind( names , subStr ) ) ) ) );
DSPnames = fieldnames(DSPstruct);

DSPidx = 1;
DSPlabel = DSPnames{DSPidx};

marmo_test_timestamps
alignTimes = timestamps;

sdf_session = SpkConvolver (spkTimes.(DSPlabel), round(max(spkTimes.(DSPlabel))+5000000), 'Gauss');



sdf_aligned = nan(length(alignTimes),range(timeWin)+1);

for ii = 2:length(alignTimes)
    if isnan(alignTimes(ii)) | alignTimes(ii) == 0
        continue
    else
        sdf_aligned(ii,:) = sdf_session(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));
    end
end

figure;
plot(timeWin,nanmean(sdf_aligned))
vline(0,'k')