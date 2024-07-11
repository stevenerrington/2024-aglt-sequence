function lfp = patch_fault_ch(lfp,fault_ch_idx)

% Note: this function is currently not configured to deal with instances in
% which the faulty contact does not have another contact adjacent on both
% sides (i.e. the top most and bottom most electrode contacts).

if ismember(fault_ch_idx,[1,16,17,32])
    fprintf('!Check function before use')
end

if any(fault_ch_idx ~= -999)
    for fault_i = 1:length(fault_ch_idx)
        fault_ch = fault_ch_idx(fault_i);
        lfp(fault_ch,:) = mean([lfp(fault_ch-1,:); lfp(fault_ch+1,:)]);
    end
end

end