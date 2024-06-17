
function figs = plot_clusters(clust,params)


figs = figure('Color','w');


nclust = clust.Nclust;

cols = cmap(nclust);

% max_traces = 100; %Choose at most this many traces for plotting
max_points = 4000; %Choose at most this many points for plotting the distribution

a=almostSquare(nclust);



for k = 1:nclust
    subplot(a(1),a(2)*3,k + 2*a(2)*floor((k-1)/a(2)))
    sp = clust.spike_sep(k);
    M  = squeeze(mean(sp.spike_waves,2));
    [~,mxi] = max(max(M)-min(M));
    
    spw = sp.spike_waves(:,:,mxi)';
    qt = quantile(spw,[.001 .999]);
    qt = minmax(qt(:));
    [H,bb] = hist(spw,linspace(qt(1),qt(2),100)); 
%     g1 = gausswin(11);
%     g2 = gausswin(ceil(.0005./mean(diff(clust.wint))));
%     G = g1*g2';
%     G = G./sum(G(:));
%     H = round(convn(H,G,'same'));
    mxH = max(H(:));
%     Wgt = (H./mxH);
    cm = flipud(gray(mxH));% + .5*(1-gray(mxH)).*(cols(k,:)+1));
    H = ind2rgb(H,cm);
%     H = Wgt.*H + (1-Wgt);
%     rp = randperm(size(sp.spike_waves,2));
%     rp = rp(1:min(max_traces,size(sp.spike_waves,2)));
%     plot(clust.wint*1e3,sp.spike_waves(:,rp),'color',[1 1 1]*.5);
    imagesc(clust.wint*1e3,bb,H)   
    axis xy
    hold on, 
    plot(clust.wint*1e3,mean(sp.spike_waves(:,:,mxi),2),'linewidth',2,'color',cols(k,:));
    if size(M,2)>1       
        title(sprintf('%i (N=%i; max ch%i)',k,sp.cluster_count,mxi))
    else
        title(sprintf('%i (N=%i)',k,sp.cluster_count))        
    end
    
    if mod(k,a(2))==1 && ceil(k/a(2))==a(1)
        xlabel('ms')
    end
    axis tight 
    grid on
    xlim([-1 1]*min([abs(xlim),5])); %% USE THIS NORMALLY
%     xlim([-1 2]);
    ylim([-1 1]*max(abs(ylim))) %% USE THIS NORMALLY
%     ylim([-150 150])
%     ylim([-0.4 0.4])
%     ylim([-1 1]*6*std(sp.spike_waves(:)));
end


%Projection to maximize separation.
MM = clust.Mx'*clust.Mx;
S = mean(clust.SS,3);
[u,~] = svd(S^-1*MM);
xsep = clust.x*u(:,1:min(3,clust.params.keepc));

if size(xsep,2)==1
    xsep(:,2) = randn(size(xsep))*std(xsep)/10;
end

subplot(2,3,2)

rp = randperm(size(xsep,1));
rp = rp(1:min(max_points,size(xsep,1)));

if size(xsep,2) ==2
     scatter(xsep(rp,1),xsep(rp,2),4,cols(clust.cl(rp),:),'markerfacecolor','flat')
else
    scatter3(xsep(rp,1),xsep(rp,2),xsep(rp,3),4,cols(clust.cl(rp),:),'markerfacecolor','flat')
end
% colormap(cols)
caxis([.5 nclust+.5]);
view(2)
set(gca,'xticklabel',[],'yticklabel',[]);
grid on
xlabel 'dim 1'
ylabel 'dim 2'
zlabel 'dim 3'
title('Clustering space (rotated PCA)')
% subplot(2,4,4)
% 
% scatter(xsep(rp,1),xsep(rp,3),4,cols(clust.cl(rp),:),'markerfacecolor','flat')
% %colormap(cols)
% %caxis([.5 nclust+.5]);
% view(2)
% set(gca,'xticklabel',[],'yticklabel',[]);
% grid on
% xlabel 'dim 1'
% ylabel 'dim 3'

% c = colorbar;
% set(c,'ytick',1:nclust);

[srt,~] = sort(clust.cl);

ax1=subplot(2,3,3);

plot_d_projection(clust,ax1,cols)

if isfield(clust.spike_sep,'active_feature')
    ax2 = subplot(2,6,11);

    A = full(cat(2,clust.spike_sep.active_feature));
    imagesc(A);
    caxis([0 quantile(A(A>0),.95)])
    hold on, 
    for k = 1:nclust
        plot(.5 +0./(srt==k),'linewidth',4,'color',cols(k,:));
        plot(size(A,1)+.5 +0./(srt==k),'linewidth',4,'color',cols(k,:));
    end
    ff = find(diff([0;srt(:);size(A,1)+1]));
    set(gca,'xtick',(ff(1:end-1)+ff(2:end))/2,'xticklabel',1:nclust,'ytick',1:size(A,1))
    axis 
    title 'HOSD Feature Scores' 
    ylabel('HOSD Feature')
    xlabel('cluster')

    
    ax4 =subplot(2,6,12);
    hos = clust.hos;

    if (isa(hos,'mvhosd')||isa(hos,'hosminimal')) && hos(1).subspace_dim>0
        for k = 1:length(hos)
            F(:,k,:) = squeeze(hos(k).feature)*hos(k).projection';
        end
       
    else
        F = [hos.feature];
    end
    F = F./nanstd(F(:));
    if size(F,3)>1
        for k = 1:size(F,2)
            imagesc(fftshift(hos(1).sampt)./hos(1).sampling_rate*1e3,.9*((0:size(F,3)-1))/size(F,3)+k-1,squeeze(F(:,k,:))');
            hold on
        end
    else    
        plot(fftshift(hos(1).sampt)./hos(1).sampling_rate*1e3,-F/20+(1:length(hos)),'k');
    end
    ylim([.5 length(hos)+.5])
    xlim([-1 1]*5)
    xlabel ms
    axis ij
    title('Feature Waveforms')
    if size(F,3)>1
        yt = .9*(0:1/size(F,3):1-1/size(F,3))'+(0:size(F,2)-1);
        yyaxis left
        set(gca,'ytick', mean(yt),'yticklabel',[])
        axis tight
        yl = ylim;
        yyaxis right
        yh = get(gca,'YAxis');

        ytl = num2str((1:size(F,3))');
        ytl(mod((1:size(F,3))',ceil(size(F,3)/4))~=0,:) = ' ';
        set(gca,'ytick',yt(:),'yticklabel',ytl)
        yh(2).Color = 'k';
        yh(2).Limits = yl;
        ylabel channel
    else
        set(gca,'ytick',[])
    end
    set(ax2,'Units','Normalized','Position',[0.6682 0.1100 0.1254 0.3412])
end
ax3 = subplot(2,3,5);
plot_isi(clust,ax3,cols)

set(ax1,'Units','Normalized','Position',[0.6625 0.5347 0.2425 0.3903])

%%% Print figure
set(figs,'unit','inches','Position',[0 0 11 8.5]*2,...
    'PaperOrientation','Landscape','PaperSize', [ 11 8.5]*2,'PaperPositionMode','auto',...
    'renderer','painter');
ttl = '';
if isfield(clust,'block')
    ttl = sprintf('Block: %s',clust.block);
end
if isfield(clust,'chan')
   ttl = sprintf('%s      Channel: %03i',ttl,clust.chan(1));
   if length(clust.chan)>1
       ttl = sprintf('%s%s',ttl,sprintf(',%03i',clust.chan(2:end)));
   end
end
if ~isempty(ttl)
    abx=annotation(figs,'textbox',...
        [0.1323 0.9580 0.1827 0.03649],...
        'String',ttl,...
        'FitBoxToText','on',...
        'FontSize',20,'edgecolor','none'); 
end

function plot_isi(clust,ax,cols)

axes(ax)
hold on
nclust = clust.Nclust;
a = almostSquare(nclust);
nx = a(2);
ny=a(1);
dx = 10;
dy = 10;
pltick = [1e-3 3e-3 10e-3 1e-1 1 10 100 1000];
qtile = [0  .1 .25 .5 .75 .9  1];
for k = 1:nclust
  
    sp = clust.spike_sep(k);
    isi = diff(sp.spike_times);
    logisi = log10(isi);
%     qt = quantile(logisi,qtile);
    
    if sp.cluster_count<10
        continue
    end
    [h,b] = hist(logisi,50);
    
     pltt = pltick(log10(pltick)>=min(b) & log10(pltick)<=max(b));
    if isempty(pltt)
        continue
    end
    plnorm = @(x)(x-min(b))./(max(b)-min(b))*.9*dx + mod(k-1,nx)*dx;
    ynorm = @(x) x*dy*.75 - floor((k-1)./nx)*dy;
    bpl = plnorm(b);
    b03 = plnorm(log10(.003));
    btick = plnorm(log10(pltt));
    
    cs=cumsum(h)./sum(h);
    plot(plnorm(b([1 end]))'*ones(size(qtile)),[1 1]'*ynorm(qtile),'Color',[1 1 1]*.5)
    plot([1 1]'*btick,ynorm([0 1]),'-','Color',[1 1 1]*.5)
    hb = bar(bpl,h./max(h),'hist');
    hb.Vertices(:,2) = ynorm(hb.Vertices(:,2));
    set(hb,'FaceColor',cols(k,:),'EdgeColor','none','FaceAlpha',1)
    plh = plot(bpl,ynorm(cs),'k','linewidth',2);
    
    intp = interp1(cs(diff(cs)>0),b(diff(cs)>0),.03);
    plot(plnorm(intp),ynorm(.03),'r.','markersize',25)
    pl1 = plot([1 1]'*plnorm(intp),ynorm([0 1]),'r--','linewidth',1);
    plot([1 1]*b03,ynorm([0 1]),'k--','linewidth',2)
    if ceil(k./nx) == ny
        for kk = 1:length(pltt)
            text(btick(kk),ynorm(-.075),sprintf('%i',pltt(kk)*1e3),'HorizontalAlignment','center','rotation',0)
        end
        if mod(k-1,nx)==0
            text(mean(bpl),ynorm(-.15),'ms','HorizontalAlignment','center','FontWeight','bold')
        end
    end
    txh = text(plnorm(intp),ynorm(1.1),sprintf('P_{.03}=%0.1f ms',10^(intp+3)),'HorizontalAlignment','center','Color','k');
    if intp>log10(.003)
        txh.FontWeight = 'Bold';
        txh.Color = 'r';
        pl1.LineWidth = 2;
    end
    if mod(k-1,nx)==0
        for kk = 1:length(qtile)
            text(plnorm(b(1))-.1,ynorm(qtile(kk)),sprintf('%i',100*qtile(kk)),'HorizontalAlignment','Center')
        end
        if ceil(k./nx) == ny
           text(plnorm(b(1))-1,ynorm(.5),'quantile (%)','HorizontalAlignment','center','FontWeight','bold','Rotation',90)
        end
    end
    axis tight
    set(gca,'xtick',[],'ytick',[])
    title('inter-spike intervals')
    ylim(ylim.*[1 1.2])
    axis off
end


% figs(2) = figure;

function plot_d_projection(clust,ax,cols)

nclust = clust.Nclust;
if nclust>1
    P = clust.x*clust.project;
else
    P = [];
end

if nargin < 3
    cols = cmap(nclust);
end

sqidx = squareform(1:nclust*(nclust-1)/2);
% subplot(1,2,1)
axes(ax)

dx = 10;
dy = 10;

hold on   
for k = 1:nclust
    for kk=k:nclust
        
%         subplot(nclust,2*nclust,2*nclust*(k-1) + kk)
%       
   
        if k~=kk
            [~,b] = hist(P(clust.cl==k|clust.cl==kk,sqidx(k,kk)),50); %#ok<*HIST>
            [h1,b1] = hist(P(clust.cl==k,sqidx(k,kk)),b);
            b1=(b1-min(b1))./(max(b1)-min(b1))*dx*.9;
            [h2,b2] = hist(P(clust.cl==kk,sqidx(k,kk)),b);
            b2=(b2-min(b2))./(max(b2)-min(b2))*dx*.9;

            m1 = mean(P(clust.cl==k,sqidx(k,kk)));
            m2 = mean(P(clust.cl==kk,sqidx(k,kk)));
            s1 = std(P(clust.cl==k,sqidx(k,kk)));
            s2 = std(P(clust.cl==kk,sqidx(k,kk)));
            dm = (m1*s2 + m2*s1)./(s1+s2);
            dm = (dm-min(b))./(max(b)-min(b))*dx*.9;
            hb1 = bar(b1 + (k-1)*dx ,h1/max(h1)*dy/2  ,'hist');
            hb1.Vertices(:,2) = hb1.Vertices(:,2)+(kk-1)*dy;
            set(hb1,'FaceColor',cols(k,:),'EdgeColor','none','FaceAlpha',1)
            hb2 = bar(b2 + (k-1)*dx ,h2/max(h2)*dy/2  ,'hist');
            hb2.Vertices(:,2) = hb2.Vertices(:,2)+(kk-1)*dy;
            plot(b2([1 end])+(k-1)*dx,(kk-1)*dy*[1 1],'k')
            plot(dm+(k-1)*dx*[1 1],[0 dy/2]+(kk-1)*dy,'k')
            set(hb2,'FaceColor',cols(kk,:),'EdgeColor','none','FaceAlpha',.8)

            txh = text(dm+(k-1)*dx*[1],(kk-.4)*dy ,sprintf('%0.1f',sqrt(clust.D(sqidx(k,kk)))),'HorizontalAlignment','Center','FontWeight','Bold');
    
        else
            [~,cdists] = clusterdists(clust.x,clust.cl);
%             u=chi2cdf(clust.d2centers(clust.cl==k),size(clust.x,2));
            u=chi2cdf(cdists(clust.cl==k),size(clust.x,2));
            uexp =(.5:length(u))'./length(u);
            udev = sort(u)-uexp;
            [KS,mxi] = max(udev*sqrt(sum(clust.cl==k)));
            [h,b] = hist(u,[0:.05:1]);
             b=b./(max(b)-min(b))*dx*.9;
             hb = bar(b + (k-1)*dx ,h/max(h)*dy/4  ,'hist');
             hb.Vertices(:,2) = hb.Vertices(:,2)+(kk-1)*dy;
             set(hb,'FaceColor',cols(k,:),'EdgeColor','none','FaceAlpha',1)
%             plot(uexp*dx*.9+ (k-1)*dx,(udev+1)*dy/2+(kk-1)*dy,'linewidth',2,'color',cols(k,:))
            plot(b([1 end])+(k-1)*dx,(kk-1)*dy*[1 1],'k')
            txh = text((k-.5)*dx*1,(kk-.4)*dy ,sprintf('KS_{\\chi^2}: %0.1f',KS),'HorizontalAlignment','Center');
        end
    end
%     if k >1
        text(-dx*.125,(k-.75)*dy,sprintf('%i',k),'FontWeight','bold','FontSize',12,'Color',cols(k,:))
%     end
%     if k < nclust
        text(dx*(k-.5),(nclust)*dy,sprintf('%i',k),'FontWeight','bold','FontSize',12,'Color',cols(k,:))
%     end  
%     axis equal
    axis([-.5*dx nclust*dx 0 (nclust+.5)*dy])
    set(gca,'xtick',[],'ytick',[])
    axis off
    title(sprintf('Cluster sep. \n & intracluster d.'))
end
        
   



function out =  almostSquare(N)

rg = 0:sqrt(N);
divs = arrayfun(@(x)find(rem(x,1:x)==0),N+rg,'uniformoutput',false);
divr = arrayfun(@(x,y)y{1}-x./y{1},N+rg,divs,'uniformoutput',false);

[mn,mni] = cellfun(@(x,y)min(abs(x)),divr,divs);
[~,mni2] = min(mn);
dm = divs{mni2}(mni(mni2));
out = [(N+rg(mni2))/dm  dm];