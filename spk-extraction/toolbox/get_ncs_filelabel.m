function filelabels = get_ncs_filelabel(directory, file_convention, n_channels)


ch = [1:n_channels];
temp = repmat({file_convention(1:3)},length(ch),1);
ext = repmat({'.ncs'},length(ch),1);
ch = cellfun(@num2str, num2cell(ch), 'UniformOutput', 0)';

file_num = file_convention(strfind(file_convention,'_'):strfind(file_convention,'.ncs')-1);
file_num = repmat({file_num},length(ch),1);
filelabels = cellfun(@horzcat, repmat({directory},length(ch),1), temp,ch,file_num,ext,'UniformOutput',false);
clearvars ch ext temp file_num


end