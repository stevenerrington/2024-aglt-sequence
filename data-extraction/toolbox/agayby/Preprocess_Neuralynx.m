%%20210409 Finished first version of the script
clc
clear all
%to make file selection with a gui
% [file,path,indx] = uigetfile('*.ncs','Select NCS file');

%For more script function
path = pwd;
file = 'CSC1_0001.ncs';

% tic

%% Getting individual file names for channels of interest
tic
cfg =[];
ch = [1:16];
temp = repmat({file(1:3)},length(ch),1);
ext = repmat({'.ncs'},length(ch),1);
ch = cellfun(@num2str, num2cell(ch), 'UniformOutput', 0)';
if isempty( strfind(file,'_') )
    cfg.fname = cellfun(@horzcat, temp,ch,ext,'UniformOutput',false);
    clearvars ch ext temp
else
    file_num = file(strfind(file,'_'):strfind(file,'.ncs')-1);
    file_num = repmat({file_num},length(ch),1);
    cfg.fname = cellfun(@horzcat, repmat({'C:\KIKUCHI-LOCAL\data\ephys\raw\2021-07-05_09-56-06_AGL2\'},length(ch),1), temp,ch,file_num,ext,'UniformOutput',false);
    clearvars ch ext temp file_num
end

% cfg.latency = [0 300];
%%
data = ft_read_neuralynx_interp(cfg.fname);
fprintf('\nChannel files loaded and processed, took %2.2f seconds\n', toc);
%%
tic
hdr   = ft_read_header(fullfile(path,file));

if isempty(strfind(file,'_'))    %% detecting if the file is the original one events.nev or the numbered ones events_xxxx.nev
    event_file = 'Events.nev';
else
    file_num = file(strfind(file,'_'):end-4);
    event_file = ['Events' file_num '.nev'];
end

event = ft_read_event_BA(event_file);

for i=1:length(event)
  % the first sample in the datafile is 1
  event(i).sample = (event(i).timestamp-double(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1;
end
fprintf('\nEvents loaded and processed, took %2.2f seconds\n', toc);

%%


%%
tic
% parts = strsplit(path, '\');
% DirPart = parts{end-1};
[~,fileNew] = fileparts(file); 
fileNew = fullfile(path,[fileNew '.mat']); %% end-4 to remove the extension
save(fileNew,'data','event','hdr')
fprintf('\nOutput saved, took %2.2f seconds\n', toc);
