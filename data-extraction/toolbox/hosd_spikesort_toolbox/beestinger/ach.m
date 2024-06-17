
sps =[];
ncl = 0;
for k = 1:length(clust)
    
    sp = clust(k).spike_sep;
    
    for kk = 1:length(sp)
        sp(kk).probe = k;
        sp(kk).cluster_orig = sp(kk).cluster;
        sp(kk).cluster = sp(kk).cluster_orig+ncl;
       
    end
    ncl = ncl+length(sp);
    sps = [sps,sp];
end

%%
fig = figure('Units','inches');
fig.Position = [0 0 8.5 11];


nperpage=5;

outdir = 'figs';
nfiles = 0;
checkfs = 1e4;

for k = 1:length(sps)
    sp = sps(srti(k));
    spts = sp.spike_times;
    q = zeros(round(checkfs*max(spts*2)),1);
    q(round(checkfs*spts))=1;
    ach = ifft(abs(fft(q)).^2);
    ach(1) = 0;
    t = ((0:length(q)-1)-floor(length(q)/2))/checkfs;
    ach = fftshift(ach);
    ach(abs(t)>1) = [];
    t(abs(t)>1) = [];
    smw = zeros(size(t(:)));
    smw(abs(t)<=.001) = 1;
    achsm = ifft(fft(ach).*fft(ifftshift(smw)));
    achsm = achsm(1:10:end);
    tsm = t(1:10:end);
    
    
    
    spw = double(sp.spike_waves);
    m =squeeze( mean(spw,2));
    [mx,mxi ] = max(diff(minmax(m)));
    
    qt = quantile(spw(:,:,mxi)',[.005 .995]);
    hrange = [-1 1]*max(abs(qt(:)));

    H = hist(spw(:,:,mxi)',hrange(1):hrange(end));
    H([1 end],:) = nan;
    
    subplot(nperpage,2,mod(k-1,nperpage)*2+1)
    imagesc(wint*1e3,hrange*0.195,1-H(:,:,[1 1 1])./max(H(:)));
    hold on
    plot(wint*1e3,m(:,mxi)*0.195,'r','linewidth',2);
    xlim([-5 5]);
    ylabel '{\mu}V'
    xlabel 'ms'
    hold off
    axis xy
    title(sprintf('Cluster %i (N=%i)\n{\\color{red} Pcnt. KS missed: %0.1f%%}',sp.cluster+(sp.probe-1)*16,sp.cluster_count,100*pnull(srti(k))));
    
%     DT = pdist(spts');
%     
%     [h,b]= hist(DT,0:1e-3:1);
    sp.ach = achsm;
    sp.achb=tsm;
    
    
     subplot(nperpage,2,mod(k-1,nperpage)*2+2)
   
%     bar([-b(end:-1:2),b]*1e3,[h(end:-1:2),h])
    bar(tsm*1e3,achsm)
    xlim([-25 25])
    xlabel 'ms'
    title('ACH')
    drawnow
    
    if mod(k,nperpage)==0
        print(fig,'-dpdf',fullfile(outdir,sprintf('clusterplot%i.pdf',nfiles)))
        nfiles = nfiles+1;
    end
    
end