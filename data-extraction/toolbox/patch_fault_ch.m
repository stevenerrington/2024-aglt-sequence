function lfp = patch_fault_ch(lfp,fault_ch_idx)

% Note: this function is currently not configured to deal with instances in
% which the faulty contact does not have another contact adjacent on both
% sides (i.e. the top most and bottom most electrode contacts).

if ismember(fault_ch_idx,[1,16,17,32])
    fprint('!Check function before use')
end

lfp(fault_ch_idx,:) = mean([lfp(fault_ch_idx-1,:); lfp(fault_ch_idx+1,:)]);


end