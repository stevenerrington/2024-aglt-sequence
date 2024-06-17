order = 3;

load selectedfeats 
load tempinput
load K3
tref = inputdat.tref;
quality = [selectedfeats.featfiltmax];
[srt,srti] = sort(quality,'descend');
afsrt = selectedfeats(srti);

hospar = [afsrt.hosparams];
%%% Use only features obtained with highpass under 1 hz.
keepfeat = [hospar.highpass]<.9;

afsrt = afsrt(keepfeat);


switch order
    case 3
        getu = [1 2 3 4 5];
        x = u3(getu,keepfeat);
    case 4
        getu = 2:4;
        x = u4(getu,keepfeat);
end        
nclust=2048;

F = [afsrt.feature];
F(isnan(F))=0;
PDF = [afsrt.pdfilt];
PDF(isnan(PDF))=0;
Fnorm = real(ifft(fft(F).*abs(fft(PDF))));
% 
% net = selforgmap([20 20],100,3,'hextop','linkdist');
% 
% [net,tr] = train(net,x);

[knnclass0,knncm0] = kmeans(x',nclust);


% y = net(x);

knnclassn0 = hist(knnclass0,1:nclust);

dd = sqrt(sum((permute(knncm0,[1 3 2])-permute(knncm0,[3 1 2])).^2,3));
sc = median(min(dd + 1e9*(dd==0)));

[knncl,hoss] = clust_algn(Fnorm,knnclass0,order);


%% Validation clustering and sign alignment for even orders
Dst = momdist(full(knncl.mFalgn),4);
[uu,ll] = svd(Dst);
% if mod(order,2)~=0
%     Dst = (Dst+1)/2;
% end
nclust2 = 256;
[scclass,V,E]  = spectralcluster(Dst,nclust2,'distance','precomputed');
jt = jet(nclust2);
% [scclass8,V8,E8]  = spectralcluster(Dst8,12,'distance','precomputed');

[sccl,shoss] = clust_algn(full(knncl.mFalgn),scclass(knnclass0),4);

if mod(order,2)==0
    Dst2 = momdist(full(sccl.mFalgn),order);
    [scclass2,V2,E2]  = spectralcluster(Dst2,1,'distance','precomputed');
    [sccl2,shoss2] = clust_algn(full(sccl.mFalgn),scclass2,order);
    sccl.D = ifft(fft(sccl2.D(:,scclass)).*fft(sccl.D));
    sccl.Falgn = ifft(fft(sccl.Falgn).*fft(sccl2.D(:,scclass)));
    
    knncl.mFalgn = ifft(fft(full(sccl.D)).*fft(full(knncl.mFalgn)));
    knncl.D = ifft(fft(full(sccl.D(:,scclass))).*fft(full(knncl.D)));
    knncl.Falgn = ifft(fft(full(sccl.D(:,knnclass))).*fft(knncl.Falgn));
end

scclm2 = []; 
for k = 1:nclust
    scclm2(:,k) = mean(x(:,scclass(knnclass0)==k),2); %#ok<*SAGROW>
end

%% Check with phase randomized data
xrand = ifft(fft(randn(size(knncl.mFalgn))).*abs(fft(full(knncl.mFalgn)))); %% Retain PSD
xrand = xrand.*sqrt(mean(abs(hilbert(full(knncl.mFalgn))).^2,2)); %Also try to retain RMS hilbert envelope.
[Dstrand,~,pkcrand] = momdist(xrand,4);
nclust2 = 256;
[scclassrand,Vrand,Erand]  = spectralcluster(Dstrand,nclust2,'distance','precomputed');
[uurand,llrand] = svd(Dstrand);
%% Mean PSDs
PSD = [afsrt.ftnorm];

%% Make plots
knnclass = scclass;
knncm = []; 
mPSD=[];
nclust = size(sccl.mFalgn,2);
for k = 1:nclust
    knncm(:,k) = mean(x(:,knnclass(scclass)==k),2);
    mPSD(:,k) = nanmean(PSD(:,knnclass(scclass)==k),2);
end
knncm(isnan(knncm))=0;
%%

% r = sqrt(knncm(getu==2,:).^2 + knncm(getu==3,:).^2); 
% mang = mod(atan2(knncm(getu==2,:),knncm(getu==3,:)),2*pi);
r = sqrt((uu(:,1)-mean(uu(:,1))).^2 +( uu(:,2)-mean(uu(:,2))).^2)'; 
mang = mod(atan2(uu(:,1)-mean(uu(:,1)),uu(:,2)-mean(uu(:,1))),2*pi)';
rr = r./max(r);
cmap = hsv2rgb([mang'/(2*pi) .9*rr' 0*rr'] + [0 .1 1]);

%   cmap=jt(scclass,:);

nodex = uu(:,1:2)'/sc;
% nodex = knncm./sc;
% nodex = [cos(mang); sin(mang)].*sqrt(rr)*10;

plotxx = sccl.mFalgn;
dxx = 1.1*diff(tref([1 end]));
dyy = diff(minmax(plotxx(:)))*.5;

wlog=linspace(log10(.25),log10(160),sum(w>0&w<160));
plotx = interp1(w,abs(fft(full(sccl.mFalgn))),10.^wlog);
dx = 1.1*diff(wlog([1 end]));
dy = diff(minmax(plotx(:)))*1;

ang=mod(atan2(x(getu==2,:)',x(getu==3,:)'),2*pi)/(2*pi);
rad=sqrt(sum(x(getu==2|getu==3,:).^2));
hsv = [ang,.1+.9*rad'/max(r), ones(size(x,2),1)];
hsv(hsv>1)=1;

cols = hsv2rgb(hsv);
figure, 
mscatter3(x(getu==2|getu==3|getu==4,:)',1,cols,'markerfaceColor','flat')
axis equal vis3d off
set(gcf,'color','black')

figure, 
hold on, 
ploti = 1:nclust;
% ploti = find(abs(ptz2)<.5)';
for k = ploti
    h = plot(tref'+dxx*nodex(1,k),dyy*nodex(2,k) +plotxx(:,k),'color',cmap(k,:));
end
axis off tight
set(gcf,'color','black')

figure, 
hold on, 
for k = ploti
    h = plot(wlog'+dx*nodex(1,k),dy*nodex(2,k) +plotx(:,k),'color',cmap(k,:));
end
axis off tight
set(gcf,'color','black')
plh = plot(wlog'+dx*nodex(1,ploti),dy*nodex(2,ploti) +.5*dy*ismember(wlog',wlog([139  365  594  ])),'color',[1 1 1]*.5);



