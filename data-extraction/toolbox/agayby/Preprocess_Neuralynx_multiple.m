%%20210409 Finished first version of the script
clc
clear all
%to make file selection with a gui
% [file,path,indx] = uigetfile('*.ncs','Select NCS file');

%For more script function
path = pwd;
% file = 'CSC1_0005.ncs';
files{1} =  'CSC1.ncs';
for ind =1:27
    files{ind+1} = [files{1}(1:end-4) '_' num2str(ind,'%.4i') '.ncs'];
end
l = length(files);
for ind =l+1:2*l
    files{ind} = ['LFP' files{ind-l}(4:end)];
end
% tic 
total_time = tic;
%% Getting individual file names for channels of interest
ch_all = [1:16;17:32];

f = waitbar(0, 'Starting');
for ind = 1:length(files)
    waitbar(ind/length(files), f, sprintf('File number: %d out of %d', ind,length(files)));
    pause(0.1);
    ylim([0.8 1.5])
    for ind2 = 1:2
        file = files{ind};
        tic
        cfg =[];
        %         ch = [1:16];
        ch = ch_all(ind2,:);
        ch_store = ch;
        
        temp = repmat({file(1:3)},length(ch),1);
        ext = repmat({'.ncs'},length(ch),1);
        ch = cellfun(@num2str, num2cell(ch), 'UniformOutput', 0)';
        if isempty( strfind(file,'_') )
            cfg.fname = cellfun(@horzcat, temp,ch,ext,'UniformOutput',false);
            clearvars ch ext temp
        else
            file_num = file(strfind(file,'_'):strfind(file,'.ncs')-1);
            file_num = repmat({file_num},length(ch),1);
            cfg.fname = cellfun(@horzcat, temp,ch,file_num,ext,'UniformOutput',false);
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
        if ch_store(1) == 1
            fileNew = [fileNew '_1'];
        elseif ch_store(1) == 17
            fileNew = [fileNew '_2'];
        end
        if ~isempty(strfind(file,'CSC'))
            cfg = [];
            cfg.resamplefs = 1000;
            data = ft_resampledata(cfg,data);
        end
        fileNew = fullfile(path,[fileNew '.mat']); %% end-4 to remove the extension
        save(fileNew,'data','event','hdr','-v7.3')
        fprintf('\nOutput saved, took %2.2f seconds\n', toc);
    end
end
total_time = toc(total_time)/60;
close(f)

fprintf('\nTotal Output saved, took %2.2f minutes\n', total_time);