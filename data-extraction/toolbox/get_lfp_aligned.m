function lfp = get_lfp_aligned(lfp_ncs_out,alignTimes,ops)

timeWin = ops.timewin;

for ch_i = 1:size(lfp_ncs_out,1)
    ch_lfp = [];
    ch_lfp = lfp_ncs_out(ch_i,:);
    filt_ch_lfp = [];
    filt_ch_lfp = lfp_filter(ch_lfp,3, 180, 1000);


    lfp_aligned = nan(length(alignTimes),range(timeWin)+1);

    for ii = 1:length(alignTimes)
        try
            lfp_aligned(ii,:) = filt_ch_lfp(alignTimes(ii)+timeWin(1):alignTimes(ii)+timeWin(end));
        end
    end

    lfp.(['lfp_' int2str(ch_i)]) = lfp_aligned;

end