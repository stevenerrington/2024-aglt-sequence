% Subject 522:
% ddir = fullfile('~/hbrl/HBRL_Analyses','Anal17_NLX-TDT_505-','Subject.522A.DBS.STN_TDT','Greenlee_Rohl');
% % File name is 522-002_SPECIALevents.mat. We will look at spike from Inpt_RZ2_chn001-Inpt_RZ2_chn003. The data is already high passed filtered.  Attached are the event time stamps for Speech production. 
% datafile = fullfile(ddir,'522-002_SPECIALevents.mat');
% outputdir = 'output';

% Subject 535:
% ddir = fullfile('~/hbrl/HBRL_Analyses','Anal17_NLX-TDT_505-550','Subject.535A.DBS.STN_TDT','Greenlee_Johari_OR_Timing');
% file name  for STN is 535-001_SPECIALevents.mat , We will look at spikes from Inpt_RZ2_chn001-Inpt_RZ2_chn003. Input 3 is in STN
% datafile = fullfile(ddir,'535-001_SPECIALevents.mat');
% and for Nigra is 535-002_SPECIALevents.mat  we will look at Inpt_RZ2_chn001-Inpt_RZ2_chn003, 
% ddir = fullfile('~/hbrl','HBRL_Upload','for_JBerger','JBData','STN','SpikeData')

if ~exist('jobindex','var')
    ddir = fullfile('~/hbrl','HBRL_Upload','for_JBerger','JBData','ClickMapping','SpikeData')
    subdirs = dir(fullfile(ddir,'*-*'));
    % datafiles = {fullfile(ddir,'273-002 ','273-002_ch001-003.mat'),fullfile(ddir,'269-001 ','269-001_ch001-003.mat'),fullfile(ddir,'253-002 ','253-002_ch001-003.mat')};
    datafiles = arrayfun(@(x)dir(fullfile(x.folder,x.name,sprintf('%s_*.mat',x.name))),subdirs,'uniformoutput',false);
    datafiles = cellfun(@(x)fullfile(x(1).folder,{x(~contains({x.name},'filt')).name}),datafiles,'uniformoutput',false);
    datafiles = [datafiles{:}];
    outputdir = fullfile(ddir,'hosd_sorting');
else
    datafiles = dir(fullfile(inputdir,'*.mat'));
    ddir = outputdir;
    datafiles = fullfile(datafiles(1).folder,{datafiles(jobindex).name});
end
%Create rarameter structure (see also HOSD_default paramters 
params = struct(...
        'windur', .015,...
        'online_plot', true,...
        'normalize_power', false,...
        'highpass', 300,...
        'ncomp', 6);

    %%% load the data

for di = 1:length(datafiles)%:-1:1
    datafile = datafiles{di};

    %%% PARAMETERS FOR HOSD
    params = HOSD_default_params(params);
    
    
    %%
    [pth,fn]= fileparts(datafile);
    blk = regexp(fn,'^\d\d\d-\d\d\d','match','once');
    try
        wh = whos('-file',datafile,'PDesMat*');
        if ~isempty(wh)
            ld = load(datafile);
            clear dat
            for k = 1:size(ld.PDesMat,1)
                dat(k).dat = double(ld.PDesMat(k,:)');
                dat(k).fs = ld.srate;
                dat(k).block = blk;
                dat(k).chan = k; %#ok<*SAGROW>
            end
        else
            ld = load(datafile,'Inpt*');
            fldn = fieldnames(ld);
            ch = regexp(fldn,'chn(\d*)','tokens','once');
            ch = str2double([ch{:}]);
            clear dat
            for k = 1:length(ch)
                dat(k).dat = ld.(fldn{k}).dat(:);
                dat(k).fs = ld.(fldn{k}).fs(1);
                dat(k).block = blk;
                dat(k).chan = ch(k);
            end
        end
    catch err
        warning(err.message)
        continue
    end
%     if isfield(dat,vname)
%         dat = dat.(vname);
%     else
%         ld = load(datafile);
%         dat.dat = double(ld.PDesMat(chan,:)');
%         dat.fs = ld.srate;
    
    try
        tmst = fullfile(ddir,'Timestamps',sprintf('%s.mat',regexprep(blk,'-','_')));
        ldts = load(tmst);
    catch
        ldts.FIDX = ld.FIDX;
    end
    for k = 1:length(dat)
        dat(k).FIDX = ldts.FIDX;
    end
%         vname = sprintf('chan00%i',chan);
%          dat.chan = chan;
%          dat.block = blk;
%     end

    for chan = 1:length(dat)
    %     chan = 1;

%         vname = sprintf('Inpt_RZ2_chn%03i',chan);
%         try
%             dat = load(datafile,vname);
%         catch
%             continue
%         end
%         if isfield(dat,vname)
%             dat = dat.(vname);
%         else
%             ld = load(datafile);
%             dat.dat = double(ld.PDesMat(chan,:)');
%             dat.fs = ld.srate;
%             tmst = fullfile(ddir,'Timestamps',sprintf('%s.mat',regexprep(blk,'-','_')));
%             ldts = load(tmst);
%             dat.FIDX = ldts.FIDX;
%             vname = sprintf('chan00%i',chan);
%              dat.chan = chan;
%              dat.block = blk;
%         end
        
        spike_detect = HOSD_spike_detection(dat(chan),params);
        
        spike_cluster = sort_spikes(spike_detect);
        
         vname = sprintf('chan00%i',dat(chan).chan);
 
        %%
        close all
        
       fig = plot_clusters(spike_cluster,params);
       
       nclust = spike_cluster.Nclust;
       cols = cmap(nclust);
        if isfield(dat,'FIDX') && ~isempty(ld.FIDX.evnt)

            [Tstim,ttstim] = chopper([-1 7],dat(chan).FIDX.time,dat(chan).fs(1));
            [unq,~,unqi] =unique(dat(chan).FIDX.evnt);
            Tstim(Tstim<1)=1;
            Tstim(Tstim>length(dat(chan).dat))=length(dat(chan).dat);
            spkimp = zeros(size(dat(chan).dat));
            spkimp(spike_detect.spike_indices) = spike_cluster.cl;
            smwin=ones(round(.005*dat(chan).fs),1);
            ximp = convn(spkimp==(1:nclust),smwin,'same')>0;
            smwin2=hann(round(.2*dat(chan).fs));
            ximp2 = convn(spkimp==(1:nclust),smwin2,'same')./sum(smwin2);
           
            fig2 = figure('Units','normalized','position',[0 0 1 1],'PaperOrientation','landscape');
     
            for kk = 1:nclust
               subplot(2,nclust,kk)
               Ximp =ximp(Tstim+(kk-1)*length(dat(chan).dat))'.*permute(cols(kk,:),[1 3 2]);
%                Ximp =ximp(Tstim+(kk-1)*length(dat(chan).dat))';
               Ximp2 =ximp2(Tstim+(kk-1)*length(dat(chan).dat))';
               imagesc(ttstim,[],Ximp)
               hold on
               plot(ttstim([1 end]),[1 1]'*find(diff([0;unqi])>0)','r--')
               axis xy
               ylabel trial
               title(sprintf('Cluster %i',kk));%,'Color',cols(kk,:))
                 subplot(2,nclust,kk+nclust)
                 hold on
                 clear plh
               for stk = 1:length(unq)

                   plh(stk) = plot(ttstim,mean(Ximp2(unqi==stk,:)).*dat(1).fs(1));
               end
               xlim(ttstim([1 end]));
               if kk==1
                 ylabel 'Mean firing rate (s^{-1})'
               xlabel 'time (s)'
               end
                 grid on
            end
            legend(plh,strcat({'Event '},num2str((unq))))
            colormap(flipud(gray))
        end
    %%
%         out.spike_times = spike_detect.spike_indices/dat(chan).fs(1);
%         out.components_active= spike;
%         out.waveform_cluster = cl; %%% This is clustering based on PCA applied to aligned waveforms.
%         out.skewness_cluster = cl3; %%% This is clustering based on peri-spike time skewness within the respec
         out.sampling_rate = dat(chan).fs(1);
        out.HOSD_waveforms = [spike_detect.hos.feature];
        out.HOSD_detection_filters = [spike_detect.hos.filterfun];
        
        out.source_file = datafile;
        out.source_variable = vname;
        out.spike_detect = spike_detect;
        out.spike_cluster = spike_cluster;
        
        fid = fopen([mfilename,'.m'],'r');
        out.COM = fread(fid,'uchar=>char')';
        fclose(fid)
        [pth,fn] = fileparts(datafile);
        outd = fullfile(outputdir,fn);
        if ~exist(outd,'dir')
            mkdir(outd)
        end

        

        outputfile = fullfile(outd,sprintf('HOSD_spike_sort_%s',vname));
        save(outputfile,'-struct','out')
        set(fig,'renderer','painters')
        print(fig,[outputfile,'.pdf'],'-r400','-bestfit','-dpdf')
        if ishandle(fig2)
            set(fig2,'renderer','painters')
            print(fig2,[outputfile,'_raster','.pdf'],'-r400','-bestfit','-dpdf')
        end
        try %Matlab is crashing on saving the fig files, for some reason.
            savefig(fig,[outputfile,'.fig']);
            savefig(fig2,[outputfile,'_raster','.fig']);
        catch
            warning('Failed to save figure')
        end
      %         close all

    end
end
