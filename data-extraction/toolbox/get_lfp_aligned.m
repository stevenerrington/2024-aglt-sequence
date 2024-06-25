function lfp = get_lfp_aligned(lfp_ncs_out,alignTimes,ops)

if isfield(ops,'freq')
    lowFreq = ops.freq(1);
    highFreq = ops.freq(2);
else
    lowFreq = 3;
    highFreq = 30;
end

timeWin = ops.timewin;

for ch_i = 1:size(lfp_ncs_out,1)
    ch_lfp = [];
    ch_lfp = lfp_ncs_out(ch_i,:);
    filt_ch_lfp = [];
    filt_ch_lfp = lfp_filter(ch_lfp,lowFreq, highFreq, 1000);


    lfp_aligned = nan(length(alignTimes),range(timeWin)+1);

    for ii = 1:length(alignTimes)
        try
            lfp_aligned(ii,:) = filt_ch_lfp(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));
        end
    end

    lfp.(['lfp_' int2str(ch_i)]) = lfp_aligned;

end