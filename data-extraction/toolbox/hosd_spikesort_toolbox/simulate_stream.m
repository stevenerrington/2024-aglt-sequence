
% Git revision $Id: 9bc91694af5c882afdd9a62f2e141b7b808fd12a $


clear all

% ncascade = 3;
ncascades = [3 3 3];
nhist = 1000;
lowpass = 4000; 
windur = .006;
burnintime = 20;
hos_learning_rate = 1e-5;
filter_adaptation_rate = .01;

n_recent_spikes = 100;

use_residual_xfilt = true; %For clustering, apply xfilt the normal way on the residual. This might cause some splitting of clusters. 

make_movie=false;
max_cycle = 1+burnintime/100; %%% How many laps to run through the data

plot_sep = 20;
run_model = 'offline';
% run_model = 'online';

spike_sorting = 'hosd';
spike_detection = 'hosd';

window_edge = .25;
wintype = @(x)tukeywin(x,window_edge);
%  wintype = 'sasaki';
% wintype = 'rectwin';

outdir = fullfile('hosd_output','online');
%%
for simulation_set = 3:-1:1
    %%
    ncascade = ncascades(simulation_set);
    for lev = 4:-1:1
        close all

    simulate_dataset = sprintf('osort/sim%i/simulation%i.mat',simulation_set*[1 1]);

    stream_data = ~exist('simulate_dataset','var')||~exist(simulate_dataset,'file');

    cols = {'b','g','r','c','m'};
    ld = load(simulate_dataset);

    [~,fn] = fileparts(simulate_dataset);
    Fs = 25000;
    streamwin = Fs*windur;
    povlp = (1-window_edge/2);
    x = zscore(ld.spiketrains{lev}');
    % x = x+randn(size(x))*1;
    [T,tt] = chopper([0 streamwin-1],1:streamwin*povlp:length(x)-streamwin,1);

    truvals = ld.spiketimes;
    nspikes = length(truvals);
    truval = zeros(size(x,1),nspikes);
    for k = 1:length(truvals)
        truval(round(truvals{k}),k)=1;
    %     TRV(:,:,k) = truval(T+(k-1)*length(x));
    end

    smwin = ones(.004*Fs,1);
    truvalsm = convn(truval,smwin,'same')>0;
    TRV = truvalsm(T(:,:,ones(1,size(truvalsm,2)))+permute(0:size(truvalsm,2)-1,[1 3 2])*length(x));
    % truespike = squeeze(any(TRV));

    H = zeros(ncascade+1,nspikes+1);

    N = round(Fs*windur);


    fig = figure('units','inches','position',[0.9583 1.2639 19.0417 9.9028]);

    btext = {'Start','Stop'};
    bh = uicontrol(fig,'units','normalized','position',[.333 .025 .095 .04],'string','Stop','Style','togglebutton','callback',@(h,varargin)set(h,'String',btext{1+h.Value}),'Value',1);

    clear hos

    for k = 1:ncascade
        hos(k) = hosobject(3);
    end
    if isa(hos,'hosobject')
    hos.initialize(N,Fs,lowpass,[],[],'window',wintype,'poverlap',povlp,'glowpass',lowpass,'use_adaptive_threshold',true,...
                  'hos_learning_rate',hos_learning_rate,'filter_adaptation_rate',filter_adaptation_rate,'filter_burnin',1./filter_adaptation_rate,'hos_burnin',1./hos_learning_rate);
    else
    hos.initialize(N,Fs,lowpass,[],[],'window',wintype,'poverlap',povlp,'glowpass',lowpass,...
                  'hos_learning_rate',hos_learning_rate,'filter_adaptation_rate',filter_adaptation_rate,'burnin',1./hos_learning_rate);
    end    
    tbuf =fftshift(hos(1).sampt)/hos(1).sampling_rate;

    wb = fftshift(hos(1).freqindx.Bfreqs{1});

    % subplot(1,3,2)
    ax1=axes('units','normalized','position',[0.1 .1 .2 .55]);
    hold on
    histy = 10*(0:ncascade-1)'+zeros(1,n_recent_spikes);
    plh0 = reshape(plot(tbuf*1e3,nan(1,length(tbuf))+histy(:)),ncascade,n_recent_spikes);
    plh1(2,:) = plot(tbuf*1e3,[hos.feature]+plot_sep*(0:ncascade-1),'b','linewidth',2);
    set(ax1,'ytick',plot_sep*(0:ncascade-1),'yticklabel',1:ncascade);

    for k = 1:length(plh1(2,:)), plh1(2,k).Color=cols{k}; set(plh0(k,:),'color',(plh1(2,k).Color+3)/4,'linewidth',.5); end
    plh1(1,:) = plot(tbuf*1e3,[hos.shiftbuffer]+10*(0:ncascade-1),'k');
%     plh1(3,:) = plot(tbuf*1e3,[hos.feature]+10*(0:ncascade-1),'k','linewidth',2);
    ylim([-plot_sep ncascade*plot_sep])
    legend(plh1(:,1),{'data','template','reconstructed signal'})
    xlabel 'time (ms)'
    ylabel('Component')
    title('Input data and HOSD templates')
    % th = title(' ','position',[.01 43]);
    % th = annotation('textbox','position',[0.02    0.9    0.1    0.1],'linestyle','none');
    th3 = annotation('textbox','position',[  0.78    0.1    0.21    0.2],'linestyle','none','String','','FontName','FixedWidth','FontSize',8,'FontWeight','Bold');

    % subplot(1,3,3)
    ax2=axes('units','normalized','position',[0.35 .1 .2 .55]);
    plh2(1,:) = plot(tbuf,[hos.shiftbuffer],'k');
    hold on
    plh2(2,:) = plot(tbuf,[hos.feature],'r');
    plh2(3,:) = plot(tbuf([1 end]),[0 0]'*[1:ncascade],'g--');
    set(ax2,'ytick',4*(0:ncascade-1),'yticklabel',1:ncascade);
    xlabel 'time (s)'
    legend(plh2(:,1),{'Filter output','Thresholded output'})
    title('Filter output and thresholding')

    for k = 1:length(hos)
        ax(k)=subplot(3,ncascade,k);
        imh(k) = imagesc(wb,wb,mean(fftshift(abs(hos(k).bicoh)),3));
        th2(k) = title('EDF: 0');
        axis xy
        axlim = axis;axlim([1 3])=0;
        axis(axlim)
        xlabel 'Freq. (Hz)'
        ylabel 'Freq. (Hz)'
    end
    colorbar(ax(end),'EastOutside')

    mscatter = @(x,varargin)scatter(x(:,1),x(:,2),varargin{:});
    % peakhist = nan(nhist,ncascade);
    event_binary = zeros(2^(ncascade),nspikes+1);
    nevnt = zeros(1,ncascade);
    confusion = zeros(nspikes+1);
    % ax3 = axes('units','normalized','position',[0.6 0.1 0.3 0.55]);
    ax3 = axes('units','normalized','position',[0.6 0.1 0.1865 0.55]);
    imh3 = imagesc(event_binary);
    caxis([0 1]);
    colorbar
    set(ax3,'ytick',1:2^ncascade,'yticklabel',dec2bin(0:2^ncascade-1),'xtick',1:nspikes+1,'xticklabel',0:nspikes);
    xlabel('True spike type')
    ylabel 'HOSD components active'
    title(sprintf('Component activation vs. true spike type (Sim = %i, Lev. = %i)',simulation_set,lev))
    for k = 1:2^ncascade
        for kk = 1:nspikes+1
            txth(k,kk) = text(kk,k,'0','horizontalalignment','center','Color',[.25 0 0],'FontSize',10);
        end
    end


    ax4 = axes('position',[ 0.83    0.4    0.15    0.25]);
    imh4 = imagesc(zeros(nspikes));
    ax4.CLim = [0 1];
    set(ax4,'ytick',1:nspikes+1,'yticklabel',0:nspikes,'xtick',1:nspikes+1,'xticklabel',1:nspikes);
    for k = 1:nspikes+1
        for kk = k:nspikes
            %if k~=kk

            txth2(k,kk) = text(kk,k,'(0%)','horizontalalignment','center','Color',[.25 0 0],'FontSize',10);
            %end
        end
    end
    axis xy
    xlabel('True spike type')
    ylabel('True spike type')
    title(sprintf('Pairwise discriminability\nKL dist (Theoretical confusion rate)'))
    cb = colorbar('SouthOutside','position',[0.8339    0.3361    0.1500    0.0133])
    xlabel(cb,'KL Divergence (bits)')
    % cb = colorbar;
    % ylabel(cb,'Information (bit)')

    Hfun = @(p)-nansum(p.*log2(p+eps));
    hfun = @(p)-p.*log2(p+eps)-(1-p).*log2(1-p+eps);
    hd = @(p)-log2(p)+log2(1-p);
    hd2 = @(p)-1./(p.*(1-p))/log(2);
    hsolv = @(H,p) p+(-hd(p) + sqrt(hd(p).^2-2*hd2(p).*(hfun(p)-H)))./hd2(p);    

    indx = 1;
    go = true;
    t0=tic;
    burninN = burnintime./windur/povlp;
    nsamp=0;
    xrecon = zeros(size(x));
    xwght = xrecon;
    win = window(wintype,size(T,1));
    fr = getframe(fig);
    %%
    while go

        if bh.Value==0
            pause(.1)
            keyboard
            continue
        end
    %    hos.get_input(x(T(:,indx)));

    %         win = window(wintype,size(T,1));
        xrem = x(T(:,indx)).*win;%.*win;
        h = zeros(size(H)-[1 0]);
        evb=0;
        evbb = zeros(size(event_binary));
        xr = 0;
        dt = zeros(size(hos));
        for k =1:length(hos)
            [xf,xfshift] = hos(k).apply_filter(xrem,false,true);
            xf = fftshift(xf);
            sdxf = std(xf);
            xthr = hos(k).filter_threshold(xf);
    %         xf = fftshift(hos(k).xfilt(xrem));
    %         xthr = fftshift(hos(k).xthresh(xrem));
            if isa(hos,'hosobject')
                xthr = xthr./sqrt(hos(k).running_var);
                xf = xf./sqrt(hos(k).running_var);
            end
            currthr = nthroot(hos(k).current_threshold,hos(k).threshold_order);
      
    %      
    %        xthr = xthr./std(xf);%./win;
    %        xf = xf./std(xf);
    %         featfilt = hos(k).xfilt(hos(k).feature);
    %         xthr = xthr./max(featfilt);%./win;
    %         xf = xf./max(featfilt);
        %         xthr = xthr + 0./xthr;
            if  any(xthr>0) || nsamp<burninN %|| hos(k).current_threshold < 1
                hos(k).do_filter_update = true;
                hos(k).do_wave_update = true;
                hos(k).do_bsp_update=true;
              
               
    %             win = circshift(win,hos(k).delay);
            else
                hos(k).do_filter_update = false;%|true;
                hos(k).do_wave_update = false;%|true;
                 hos(k).do_bsp_update=true;
                 hos(k).do_CDF_update=true;
            end
    %           hos(k).use_adaptive_threshold=true;
            hos(k).get_input(xrem);
    %         xr0 = hos(k).xrec(xrem.*win);
            xr0 = hos(k).xrec(xrem);
            xrem = xrem-xr0;
            xr=xr0+xr;

              if  nsamp>burninN && any(xthr>0)
                 plh0(k,mod(nevnt(k),n_recent_spikes)+1).YData =  hos(k).xrec(hos(k).shiftbuffer) + plot_sep*(k-1);
                 nevnt(k) = nevnt(k)+1;
            end
            evb = evb+2^(k-1)*any(xthr>0);

    %         sdxfs(k) = sdxf;
             xfs(:,k) = xf./currthr;
    %         mxff(k) = max(featfilt);
            xthrs(:,k) = xthr./currthr;%./(win+eps);
    %         h(k+1,1) = ~any(tru)*any(xthr>0);
    %         h(k+1,2:end)= any(xthr>0)*truespike(indx,:);
    %         h = h*(hos(k).EDF>500);
    %      
           xfx = hos(k).xfilt(hos(1).feature); % Align component features on peak of the detection filter applied to the first feature
%            dt(k) = hos(1).sampt*(abs(xfx)==max(abs(xfx)));
%            hos(k).lag= hos(k).lag.*exp(-1i*2*pi*dt(k)./hos(1).buffersize);
        end

        xrecon(T(:,indx)) = (xrecon(T(:,indx)).*xwght(T(:,indx))+xr.*win)./(xwght(T(:,indx))+win);
        xwght(T(:,indx)) = xwght(T(:,indx))+win;

         hosdelay = [hos.delay].*any(xthrs>0);
         Tsh = T(Fs*.002:end-Fs*.002,indx)+permute(hosdelay,[1 3 2]) + (0:nspikes-1)*length(x);
         Tsh(Tsh<1)=1;Tsh(Tsh>size(truvalsm,2)*length(x)) = size(truvalsm,2)*length(x);
         tru = any(any(truvalsm(Tsh)),3);
         tru = squeeze([~any(tru),tru]);
    %      tru = [~any(truespike(indx,:)),truespike(indx,:)];
         evbb(evb+1,tru) =   1./sum(tru);

    %       if evbb(1,2)>0 %&&bh.Value==0
    %           keyboard
    %      end
    %     if any(evbb(1,2:end)) && nsamp>burninN
    %             keyboard
    %     end
    %      event_binary =   event_binary*(1-hos(1).filter_adaptation_rate)+hos(1).filter_adaptation_rate*evbb;
         event_binary =   event_binary + (nsamp>burninN)*evbb;
         Pjoint = (event_binary+eps)./sum(event_binary(:)+eps);
         Pevent = Pjoint./sum(Pjoint); %Use Prob(respone | spike type)
         Pclass = Pjoint./sum(Pjoint,2); %Prob(spike type | response)
         [~,mlclass]  = max(Pclass(evb+1,:));
         confusion(mlclass,tru)=confusion(mlclass,tru)+(nsamp>burninN);

         if mod(indx,100)==0

    %%
            for k = 1:length(hos)
                imh(k).CData = mean(abs(fftshift(abs(hos(k).bicoh))),3);
                if k==1
                    cax=ax(1).CLim;
                else
                    ax(k).CLim=cax;
                end

    %             xrec = circshift(hos(k).xrec(xrem),hos(k).delay);
                if any(xthrs(:))
%                     xrec = hos(k).xrec(hos(k).shiftbuffer);
                else
                    xrec(:)=0;
                end
                plh1(1,k).YData = hos(k).shiftbuffer+plot_sep*(k-1);
                plh1(2,k).YData = hos(k).feature+plot_sep*(k-1);
%                 plh1(3,k).YData = xrec+plot_sep*(k-1);
                plh2(1,k).YData = xfs(:,k)+4*(k-1);
                plh2(2,k).YData = xthrs(:,k)+4*(k-1);
    %            plh2(3,k).YData = nthroot(hos(k).current_threshold,3)*[1 1]+plot_sep*(k-1);  
               plh2(3,k).YData = [1 1]+4*(k-1);  
    %             plh2(3,k).YData = nthroot(hos(k).current_threshold,3)*sdxfs(k)./mxff(k)*[1 1]+plot_sep*(k-1);  

                set(th2(k), 'string',sprintf('Resid. Bicoherence %i EDF: %0.1f',k,hos(k).EDF));



            end
            set(th2(1), 'string',sprintf('Full Bicoherence EDF: %0.1f',hos(1).EDF));

          set(th3, 'string',sprintf('Simulation set %i\nNoise level: %i:\nData time: %0.2f\nReal time: %0.2f',simulation_set,lev,    mean(T([1 end],indx))./Fs,toc(t0)));

           if nsamp>burninN

            txth(1).String = sprintf('TN\n%i\n(%0.0f%%)',round(event_binary(1)),100*Pevent(1));
            txth(1).Color=[1 1 1]-(Pevent(1)>.25);
            for k = 2:2.^ncascade
                txth(k,1).String = sprintf('FP\n%i\n(%0.0f%%)',round(event_binary(k,1)),100*Pevent(k,1));
                 txth(k,1).Color=[1 1 1]-(Pevent(k,1)>.25);
                 for kk = 2:nspikes+1
                    txth(k,kk).String = sprintf('%i\n(%0.0f%%)',round(event_binary(k,kk)),100*Pevent(k,kk));
                   txth(k,kk).Color=[1 1 1]-(Pevent(k,kk)>.25);
                end
            end
            for k = 2:size(truvals,2)+1
                txth(1,k).String = sprintf('FN\n%i\n(%0.0f%%)',round(event_binary(1,k)),100*Pevent(1,k));
                txth(1,k).Color=[1 1 1]-(Pevent(1,k)>.25);

            end
                 misclass = confusion.*(1-eye(size(confusion)));
                str = round([0:nspikes;
                             confusion(1,:);
                             100*confusion(1,:)./sum(confusion);
                             sum(misclass(2:end,:));
                             100*sum(misclass(2:end,:))./sum(confusion);
                             sum(confusion(2:end,:));
                             100*(sum(confusion(2:end,:)./sum(confusion)));
                             sum(confusion-misclass);
                             100*sum(confusion-misclass)./sum(confusion)]);
    %             str = round([0:nspikes;event_binary(1,:);100*event_binary(1,:)./sum(event_binary);sum(event_binary(2:end,:)); 100*(1-event_binary(1,:)./sum(event_binary))]);
    %             th3.String = sprintf('%s\n\n        Overall Rates (for unsorted spikes)\n\n               FP              TN\nSpike 0:    %i (%i%%)       %i (%i%%)\n\n              misses (FN)   detections (TP)\n%s-----------------------------------------\nAll  :%s',...
    %                 sprintf('Simulation set %i\nNoise level: %i:\nData time: %0.2f\nReal time: %0.2f',simulation_set,lev,    mean(T([1 end],indx))./Fs,toc(t0)),...
    %                 str([4 5 2 3],1),...
    %                 sprintf('Spike %i:    %5i (%3i%%)    %5i (%3i%%)\n',str(:,2:end)),...  
    %                 sprintf('      %5i (%3i%%)    %5i (%3i%%)',round(sum(event_binary(1,2:end))),round(100*sum(sum(Pjoint(1,2:end)))./sum(sum(Pjoint(1:end,2:end)))),round(sum(sum(event_binary(2:end,2:end)))),round(100*sum(sum(Pjoint(2:end,2:end)))./sum(sum(Pjoint(1:end,2:end))))));
                th3.String = sprintf('%s\n\n        Classification based on threshod crossing\n\n               FP              TN\nSpike 0:    %i (%i%%)       %i (%i%%)\n\n            Misses (FN)         |        Hits(TP)\n           noise      misclass. |   unsorted    class \n%s----------------------------------------------------------\nAll  :%s',...
                    sprintf('Simulation set %i\nNoise level: %i:\nData time: %0.2f\nReal time: %0.2f',simulation_set,lev,    mean(T([1 end],indx))./Fs,toc(t0)),...
                    str([4 5 2 3],1),...
                    sprintf('Spike %i:%5i(%3i%%) %5i(%3i%%) |%5i(%3i%%) %5i(%3i%%)\n',str(:,2:end)),...  
                    sprintf('  %5i(%3i%%) %5i(%3i%%) |%5i(%3i%%) %5i(%3i%%)',sum(confusion(1,2:end)),round(100*sum(confusion(1,2:end))./sum(sum(confusion(1:end,2:end)))),...
                                                                                 round(sum(sum(misclass(2:end,2:end)))),round(100*sum(sum(misclass(2:end,2:end)))./sum(sum(confusion(2:end,2:end)))),...
                                                                                 sum(sum(confusion(2:end,2:end))),round(100*sum(sum(confusion(2:end,2:end)))/sum(sum(confusion(1:end,2:end)))),...
                                                                                 sum(sum(confusion(1:end,2:end)-misclass(:,2:end))),round(100*sum(sum(confusion(1:end,2:end)-misclass(:,2:end)))/sum(sum(confusion(1:end,2:end))))));
    %         shg
               imh3.CData = event_binary./sum(event_binary);
        %         ax3.CLim = [0 max(event_binary(2:end))+eps];
                set(plh2(1:2,:),'Xdata',(tt+T(1,indx))/Fs);
                set(plh2(3,:),'Xdata',(tt([1 end])+T(1,indx))/Fs);
                ax2.XLim = (tt([1 end])+T(1,indx))/Fs;


                Pmarginal = Pevent+permute(Pevent,[1 3 2]);
                Ppair = (Pevent)./sum(Pmarginal+eps);
                Ppaircond = (Pevent)./(Pmarginal);
                Pmarginal = Pmarginal./sum(Pmarginal);
                Hpaircond = squeeze(nansum(hfun(Ppaircond).*Pmarginal));
                Hpairmarginal = squeeze(hfun(nansum(Ppair)));
                KLdistpair = Hpairmarginal-Hpaircond;
                KLdistpairnorm = KLdistpair./Hpairmarginal;
                imh4.CData = KLdistpairnorm(1:end-1,2:end);
                Pconfusion=.3;
                for k = 1:10
                    Pconfusion = real(hsolv(Hpaircond,Pconfusion));
                end

                for k = 1:nspikes
                    for kk = k:nspikes
    %                     if k~=kk
        %                 txth2(k,kk).String = sprintf('(%0.0f%%)',round(KLdistpairnorm(k,kk)*100));
                            imh4.CData(kk+1,k) =nan;
                          txth2(k,kk).String = sprintf('(%0.0f%%)',round(Pconfusion(k,kk+1)*100));
                          txth2(k,kk).Color = [1 1 1]-(KLdistpairnorm(k,kk+1)>.25);
    %                     end
                    end
                end

                Poutput = sum(Pjoint,2);
                Ptype_cond_output = Pjoint./Poutput;
           end
    %               set(th, 'string',sprintf('Simulation set %i\nNoise level: %i:\nData time: %0.2f\nReal time: %0.2f',simulation_set,lev,    mean(T([1 end],indx))./Fs,toc(t0)),'fontsize',12);
            if make_movie
                fr(end+1) = getframe(fig);
            end   
              drawnow
               hos(1).current_threshold
         end

    %     h(1,:) = ~any(h);
    %     H = H + h;
        indx = mod(indx,size(T,2))+1;

        nsamp=nsamp+1;

        go = nsamp<max_cycle*size(T,2);
    %     vars(nsamp,:) = [hos.running_var];

    end

     if make_movie

         avi = VideoWriter(fullfile(outdir,sprintf('online_sim%i_lev%i_ncomp%i.mp4',simulation_set,lev,ncascade)),'MPEG-4');
         avi.FrameRate=4;
         avi.Quality = 25;
         avi.open;
         avi.writeVideo(fr([end 1:end]));
         avi.close();

         set(fig,'Units','Inches');
        pos = get(fig,'Position');
        set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
         print(fig,'-dpdf',fullfile(outdir,sprintf('online_sim%i_lev%i_ncomp%i.pdf',simulation_set,lev,ncascade)))
     end


    %%

        hosoffline = hosobject(hos);
        hosoffline.reset;
        for k = 1:length(hosoffline)
            hosoffline(k).do_update = true;
            hosoffline(k).do_filter_update=true;
            hosoffline(k).do_wave_update=true;
            hosoffline(k).use_adaptive_threshold=false;
%             hosoffline(k).poverlap=.5;
%             hosoffline(k).window='sasaki';

        end

        hosoffline.get_block(x(T));

        %%% Align component features to the first feature according to the peak 
        %%% of the detection filter applied to the 1st feature 
        %%% using the window circular shift parameter to the first feature
        %%% This should reduce the likelihood of spike times splitting into
        %%% into multiple times across the different components.
        for k = 1:length(hosoffline)   
            xfx = hosoffline(k).xfilt(hosoffline(1).feature); 
             dt(k) = hosoffline(1).sampt*(abs(xfx)==max(abs(xfx)));
               hosoffline(k).lag= hosoffline(k).lag.*exp(-1i*2*pi*dt(k)./hos(1).buffersize);
         end 
        % segs = [hosoffline.segment];
        % T2 = chopper(segs(1).Trange,[segs(1).wint],segs(1).fs);
        T2 = T(:,:,ones(1,ncascade))+ permute(0:length(hosoffline)-1,[1 3 2])*length(x);



        T2adj = T2 + cat(3,hosoffline.delay) ;
        T2adj(T2adj<1) = 1;
        T2adj(T2adj>ncascade*length(x))=ncascade*length(x);
        for k = 1:length(hos)
            hos(k).do_update=false;
            hos(k).use_adaptive_threshold=false;
        end
        xfilt = hosoffline.xfilt(x);
         xthr = hosoffline.xthresh(x);
        % xthr = hos.xthresh(x);
        xthrsm = sum(convn(full(xthr),hann(Fs*.001),'same'),2);
        Xfiltadj =  reshape(permute(xfilt(T2adj),[1 3 2]),size(T2adj,1)*size(T2adj,3),size(T2adj,2));
        [~,pk] = getpeak2(xthr);
        [~,pk2] = getpeak2(sum(xthrsm,2));
        pksm = convn(full(pk==1),ones(round(Fs*.004),1),'same')>0;
        pkcat = (pk==1)*2.^(0:ncascade-1)';
        pkt = find(pk2==1);
        
        % Use realignHOSD to improve alignment
        [Xrl,newTimes,tshift,clusterOut] = realignSpikesHOSD(x,pkt,struct('samplingFreq',Fs,'detectionParams',struct('hos',hosobject(hosoffline))));
        pkt = round(newTimes(~clusterOut.discarded));
        cl = clusterOut.clust(~clusterOut.discarded);
        [Tpk,ttpk] = chopper(minmax(hos(1).sampt)/Fs,pkt/Fs,Fs);
        Tpk(Tpk<1)=1;
        Tpk(Tpk>length(x))=length(x);
       
%         %Align to the center of gravity of the half-wave rectified signal
%         %to reduce splitting that results from threshodling and peak detection.
%         sxthr = sum(xfilt.^2.*(xfilt>0),2);
%         sxthr = sxthr(Tpk(abs(ttpk)<1e-3,:,1)).*hann(sum(abs(ttpk)<1e-3));
%         tshift = round((ttpk(abs(ttpk)<1e-3)'*sxthr)./sum(sxthr)*hosoffline(1).sampling_rate);        
%             Tpk = round(Tpk+tshift);
        
        pktype = pksm(pkt,:)*2.^(0:ncascade-1)';
        [srtpk,srtipk] = sort(pktype);

        % [Tpk,ttpk] = chopper([-1 1]*windur/2,pkt/Fs,Fs);
%         xfx = hos(1).xfilt([hos.feature]);
%         dts = hos(1).sampt*(abs(xfx)==max(abs(xfx)))./sum(abs(xfx)==max(abs(xfx)));
        Tpk = Tpk + length(x)*permute(0:max(nspikes,ncascade)-1,[1 3 2]);

        %         hos2 = hosobject(hosoffline);
%         A = hos2.get_block(x(Tpk(:,:,1)));
%         Tpk2 = Tpk(:,:,1)+cat(3,hos2.delay);
%         AA = reshape(permute(A,[1 3 2]),ncascade*size(Tpk,1),size(Tpk,2));
        if use_residual_xfilt
%             xfilt = sum(hosoffline.xfilt(x),2);
            xfilt = zscore(hosoffline.xfilt(x));
        else
            for k = 1:length(hos)
                xfilt(:,k) = zscore(hosoffline(k).xfilt(x));
            end
        end
%         B = xfilt(Tpk(:,:,1)+cat(3,hos2.delay)).*sasaki(size(Tpk,1));
%     %     B = xfilt(Tpkx).*hamming(size(Tpk,1));
%         BB = reshape(permute(B,[1 3 2]),ncascade*size(Tpk,1),size(Tpk,2));
%         Xpkfilt = fftshift(hosoffline.xfilt(x(Tpk(:,:,1)).*hamming(size(Tpk,1))),1);
        Xpkfilt = xfilt(Tpk(:,:,1:ncascade)).*hamming(size(Tpk,1));
        XXpkf = reshape(permute((Xpkfilt),[1 3 2]),ncascade*size(Tpk,1),size(Tpk,2));
        Xpk = x(Tpk(:,:,1)).*hamming(size(Tpk,1));

        trupk = squeeze(truvalsm(Tpk(round(end/2),:,1:nspikes)))*2.^(0:nspikes-1)';
        [srttrupk,srtitrupk] = sort(trupk);
        [a,b,c,d] =crosstab(pktype,trupk);

        % q = crossn(2*ones(1,ncascade))-1;
        % val = q*2.^(0:ncascade-1)'; [srt,srti] = sort(val); q=q(srti,:);
        % % q = q-q(2,:);
        % r = splitreg(regressor(repmat(q,length(pktype),1),'noptions',size(q,1)));
        % for k = 1:length(r)
        %     r(k).label = dec2bin(2^(k-1),ncascade);
        % end
        % rint = interaction(r,'intxnord',length(r),'codeincr',r(end).code);
        % rr = [r,rint];
        % out = modelFit(pktype+1,rr);

        %Xpk = x(Tpk(:,:,1));
%          [u,l,v] = svd(XXpkf);
    %      [u,l,v] = svd(Xfiltadj);
          [u,l,v] = svd(Xpk);

        xfilt =fftshift(hosoffline.xfilt(x(Tpk(:,:,1))),1);
        % xsk = squeeze(skewness(xfilt(T2)));f
        % sigsk = mean(xsk)./std(xsk)*sqrt(size(T,2)*min(povlp,1));

        hos.get_input(x(T2(:,:,1)));
        T3adj = T2(:,:,ones(1,ncascade)) + cat(3,hos.delay) + permute(0:length(hos)-1,[1 3 2])*length(x);
        T3adj(T3adj<1)=1;
        T3adj(T3adj>ncascade*length(x))=ncascade*length(x);
        ximp = full(hosoffline.ximp(x(:)));
        ximp2 = full(hos.ximp(x(:)));
        [srt,srti] = sort(squeeze(ximp(T2adj(round(size(T,1)/2),:,:)))*2.^(0:ncascade-1)');
        [srt2,srti2] = sort(squeeze(ximp2(T3adj(round(size(T,1)/2),:,:)))*2.^(0:ncascade-1)');
        [Ttru,ttrue] = chopper([0 streamwin-1]-(streamwin-1)/2,find(any(truval,2)),1);
        Ttru(Ttru<1)=1;Ttru(Ttru>length(x))=length(x);
        Ttru = Ttru(:,:,ones(1,max(nspikes,ncascade)))  + permute(0:max(nspikes,ncascade)-1,[1 3 2])*length(x);
        % [srttru,srtitru] = sort(squeeze(mode(truvalsm(T2adj([-10:10]+round(size(T,1)/2),:,:))))*2.^(0:2)');

        trv = truval(any(truval,2),:)*(1:nspikes)';
        [srttv,srtitv] = sort(trv);

        misses = sum(pksm(Ttru(end/2,:,1:ncascade)),3)==0;
        misstab = hist3([misses',trv],{0:1,1:nspikes});
    %     [srttru,srtitru] = sort(squeeze((truvalsm(T2adj(round(end/2),:,1:nspikes))))*2.^(0:nspikes-1)');
        %%
        fig2 = figure('units','normalize','position',[0.1365 0.5267 0.5557 0.3875]);

        cm =[0     0     0
             0     0     1     
             0     1     0
             1     0     0
             0     1     1
             1     0     1];
        for k = 0:nspikes-1
            for kk = k+1:nspikes

                if k>0

                    subtp = (trupk==2^(k-1))+2*(trupk==2^(kk-1));
                else
                    subtp = (trupk==0)+2*(trupk==2^(kk-1));
                end
                vv = v(subtp>0,1:nspikes);
                vv=vv-mean(vv);
                tp = subtp(subtp>0);
                ms = [mean(vv(tp==1,:),1);mean(vv(tp==2,:),1)];
                S = cov(vv-ms(tp,:));
                dproj = S^-.5*diff(ms)';

                D(k+1,kk) = dproj'*dproj;
                vproj = vv*dproj./sqrt(sum(diff(ms).^2));

                ax=subplot(nspikes,2*(nspikes),k+1+(kk-1)*2*(nspikes),'ytick',[]);
                hh1(kk,k+1)=histogram(vproj(tp==1),50,'FaceColor',cm(k+1,:),'EdgeColor','none');
                hold on
                hh2(kk,k+1)=histogram(vproj(tp==2),50,'FaceColor',cm(kk+1,:),'EdgeColor','none');
                title(sprintf('$d=%0.1f$',sqrt(D(k+1,kk))),'interpreter','latex');
                if k==0
                    yl = ylabel(sprintf('Spike %i\n%i%% det.',kk,round(100*misstab(1,kk)/sum(misstab(:,kk)))),'Rotation',90)
                    yl.Position(1) = yl.Position(1)-diff(xlim)*.1;
                end
                if kk==nspikes
                    if k==0
                    xlabel(sprintf('Spike %i (noise)',k))
                    else                
                    xlabel(sprintf('Spike %i',k))
                    end
                end
        %         ax.YTickLabel='';
            grid on
            end
        end
        allms = [];
        for k = [0,1:2^nspikes]
            allms(k+1,:) = mean(v(trupk==k,1:nspikes),1);
        end
        allS = cov(v(trupk>0,1:nspikes)-allms(trupk(trupk>0)+1,:));

        H = allS^-1*cov(allms,'omitrows');
        [mu,ml] = svd(H);
        mscatter3 = @(x,varargin)scatter3(x(:,1),x(:,2),x(:,3),varargin{:});
        mscatter = @(x,varargin)scatter(x(:,1),x(:,2),varargin{:});
      
        for subi =  1:2
        subplot(2,2,2*subi)
        q=floor(log2(trupk))+1;
        q(isinf(q))=0;

        mscatter(v(:,[1 1+subi]),5,q);
        colormap(cm(1:nspikes+1,:))
        caxis([-.5 nspikes+.5])
        grid on
        if subi==1
            title('Offline performance vs. ground truth using PCA of peri-detection filter outputs')
        end
        xlabel(sprintf('PC%i',1));
        ylabel(sprintf('PC%i',1+subi));
        fig2.Units = 'inches';
        ch = colorbar;
         set(ch,'ytick',0:5,'yticklabel',{'Spike 0 (noise)','Spike 1','Spike 2','Spike 3','Spike 4','Spike 5'})
        end
          pos = get(fig2,'Position');
       set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
       th4 = annotation(fig2,'textbox','position',[ 0.2848    0.7021    0.21    0.2],'linestyle','none','String','','FontName','FixedWidth','FontSize',12,'FontWeight','Bold');
       th4.String = sprintf('Simulation set %i\nNoise level:%i\nNumber of spikes: %i\nNumber of HOSD components: %i',simulation_set,lev,nspikes, ncascade) 
       th5 = annotation(fig2,'textbox','position',[  0.381  0.427    0.21    0.2],'linestyle','none','String','','FontName','FixedWidth','FontSize',8,'FontWeight','Bold');
       strtab = [0:nspikes;sum(trupk==0),misstab(2,:);length(trupk),sum(misstab);round(100*sum(trupk==0)/length(trupk)),round(100*misstab(2,:)./sum(misstab))];
        th5.String = sprintf('--False pos.--\n%3i/%3i (%i%%)\n\n   --Misses--\n%s',strtab(2:end,1),sprintf('spike %i: %3i/%3i (%i%%)\n',strtab(:,2:end)));
         print(fig2,'-dpdf',fullfile(outdir,sprintf('offline_sim%i_lev%i_ncomp%i.pdf',simulation_set,lev,ncascade)))
%     try
        %%
        if true
            prfigs = [13 14 15 8888];
            fldn = {'hosd','wvlt'};
            tmpdir1 = tempname;

            vnames = {'\#Spikes','\#Detected','\#algn','\#Det-algn','TP','FP','FPnoise','FPothclus','det-TP','nrAssgn'};
            pdfnames = {};
            mdtype = [3 1];
            for deti = 1:length(fldn)
                [osort(simulation_set,lev).(fldn{deti}).perf,osort(simulation_set,lev).(fldn{deti}).nrAssigned,osort(simulation_set,lev).(fldn{deti}).assigned,osort(simulation_set,lev).(fldn{deti}).params] = runSimulatedEval(simulation_set, lev, mdtype(deti),false, true);   

    %        [osort(simulation_set,lev).wvlt.perf,osort(simulation_set,lev).wvlt.nrAssigned,osort(simulation_set,lev).wvlt.assigned,osort(simulation_set,lev).wvlt.params] = runSimulatedEval(simulation_set, lev, 1,false, true);   
                mkdir(fullfile(tmpdir1,fldn{deti}));
                fg=figure(81218)
                cla
                text(.5,.5,sprintf('Performance for %s',fldn{deti}),'Fontsize',20,'horizontalalignment','center')
                axis off
                pdfnames{end+1}=fullfile(tmpdir1,fldn{deti},sprintf('fig%i.eps',0));
                print(fg,'-depsc', pdfnames{end});
                for k = 1:length(prfigs)
                 pdfnames{end+1}=fullfile(tmpdir1,fldn{deti},sprintf('fig%i.eps',k));
                 print(prfigs(k),'-depsc', pdfnames{end});
                end
            end
            fg = figure(12321);

            x1=[1:nspikes;osort(simulation_set,lev).hosd.perf'];
            x2=[1:nspikes;osort(simulation_set,lev).wvlt.perf'];

                 txt1 = sprintf('Sim. set:%i\nNoise lev: %i\n\n\\begin{tabular}{l%s}\\hline\\multicolumn{%i}{c}{\\textbf{HOSD Detection}} \\\\ \\hline &%s \\\\ \\hline %s \\end{tabular}',...
              simulation_set,lev,repmat('|l',1,length(vnames)+1),length(vnames)+1,sprintf(' %s & ',vnames{:}),...
              sprintf([' spike %i & ',repmat(' %i &',1,length(vnames)),' \\\\ '],x1(:)));
          txt2 = sprintf('\\begin{tabular}{l%s}\\hline\\multicolumn{%i}{c}{\\textbf{Wavelet Detection}} \\\\ \\hline &%s \\\\ \\hline %s \\end{tabular}',...
              repmat('|l',1,length(vnames)+1),length(vnames)+1,sprintf(' %s & ',vnames{:}),...
              sprintf([' spike %i & ',repmat(' %i &',1,length(vnames)),' \\\\ '],x2(:)));
       

          th2(1) = text(-.1,.66,txt1,'interpreter','latex')  ; 
          th2(2) = text(-.1,.33,txt2,'interpreter','latex')  ; 
              axis off
              pdfnames{end+1}=fullfile(tmpdir1,'perf.eps');
              print(fg,'-depsc', pdfnames{end});
              outfn = fullfile(outdir,sprintf('summary_set%i_lev%i.pdf',simulation_set,lev));
             com1 = sprintf('/usr/local/bin/gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s %s',...
                               outfn,sprintf(' %s',pdfnames{[end,1:end-1]}));

               com2=sprintf('/Library/TeX/texbin/pdfcrop %s %s',outfn,outfn);
    %  [err2,out2]=system(com2);

            fid = fopen('temp.sh','w');
            fprintf(fid,'LD_LIBRARY_PATH=\n%s\n\n%s',com1,com2);
            fclose(fid);
            system('chmod 755 temp.sh');
    %         [err,out] = system('/opt/X11/bin/xterm -e "./temp.sh;read;"');
            [err,out] = system('./temp.sh');
        end
%          system(com1);
        
        
%       catch
%         
%     end
    end
end
%%

Pevent = event_binary./sum(event_binary(:));
Pmarginal = Pevent+permute(Pevent,[1 3 2]);
Ppair = Pevent./sum(Pmarginal);
Ppaircond = Pevent./Pmarginal;
Hpaircond = squeeze(sum(hfun(Ppaircond).*Pmarginal));
Hpairmarginal = squeeze(hfun(sum(Ppair)));
KLdistpair = Hpairmarginal-Hpaircond;
KLdistpairnorm = KLdistpair./Hpairmarginal;Pmarginal = Pmarginal./sum(Pmarginal);



Pconfusion=.3;
for k = 1:10
    Pconfusion = real(hsolv(Hpaircond,Pconfusion));
end

PPair2 = Pevent./(Pevent+permute(Pevent,[1 3 2]));
PPairmarginal = sum(Pevent)./(sum(Pevent)+sum(Pevent)');
Hpair2 = -PPair2.*log2(PPair2)-(1-PPair2).*log2(1-PPair2);

Hmarginal = -squeeze(sum(Pmarginal.*log2(Pmarginal+eps)));
Hconditional = -squeeze(sum(Ppair.*(log2(Ppair+eps)-log2(sum(Ppair)+eps))));
KLdist =squeeze( Hmarginal-Hconditional-Hconditional');
KLdistNorm = KLdist./Hpair2;

figure, imagesc(100*KLdistpairnorm)
title('Normalized KL distance (% entropy decrease in classifier output)')
set(gca,'xtick',1:nspikes,'xticklabel',0:nspikes-1);
set(gca,'ytick',1:nspikes,'yticklabel',0:nspikes-1);
caxis([0 100]);
