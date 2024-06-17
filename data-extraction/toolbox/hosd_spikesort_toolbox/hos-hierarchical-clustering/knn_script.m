order = 4;

% load selectedfeats 
load tempinput
tref = inputdat.tref;
w = (0:length(tref)-1)/length(tref)*500;
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
nclust=400;

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

[knncl0,hoss0] = clust_algn(Fnorm,knnclass0,order);

%% Construct a 2D surface by fitting the PC2-PC3 plane into PC4
[knnclass2,knncm2] = kmeans(x(1:2,:)',nclust);

knncm2 = []; 
for k = 1:nclust
    knncm2(:,k) = mean(x(:,knnclass2==k),2);
end
knncm2(isnan(knncm2))=0;
tri = delaunay(knncm2(1:2,:)');
surf2 = triangulation(tri,knncm2(1:3,:)');
edg = surf2.edges;
edglen= sqrt(sum((knncm2(:,edg(:,1)) - knncm2(:,edg(:,2))).^2));
%%% Cut edges whose length is above a certain threshold
rmtri = surf2.edgeAttachments(edg(isnan(iterz(edglen,5)),:));
tri([rmtri{:}],:)=[];
surf = triangulation(tri,knncm2(1:3,:)');
r2 = sqrt(knncm2(getu==2,:).^2 + knncm2(getu==3,:).^2); 
rr2 = r2./max(r2);
cmap2 = hsv2rgb([mod(atan2(knncm2(getu==2,:),knncm2(getu==3,:))-0,2*pi)'/(2*pi) .9*rr2' 0*rr2'] + [0 .1 1]);

%%%Compute local stats
nn = surf.nearestNeighbor(x(1:3,:)');
vnm = surf.vertexNormal(nn);
ptdev = sum(vnm.*(x(1:3,:)'-surf.Points(nn,:)),2);
for k = 1:size(surf.Points,1)
    vsd(k) = sqrt(mean(ptdev(nn==k).^2));
end
ptz = ptdev./vsd(nn)';

nn2 = surf.nearestNeighbor(knncm2(1:3,:)');
vnm2 = surf.vertexNormal(nn2);
ptdev2 = sum(vnm2.*(knncm2(1:3,:)'-surf.Points(nn2,:)),2);
ptz2 = ptdev2./vsd(nn2)';


%%% In order to obtain a more uniform sammpling of the surface which also respects the z axis, project a
%%% subset of representative points onto the surface and rerun kmeans

kppts = abs(ptz)<Inf; %Retain points representative of the surface

surfprojpts = x'-ptdev.*vnm;

knnclass = zeros(1,size(Fnorm,2));
[knnclass(kppts),knncm] = kmeans(surfprojpts(kppts,:),nclust);


[knncl,hoss] = clust_algn(Fnorm,knnclass,order);


%% Validation clustering and sign alignment for even orders
Dst = momdist(full(knncl.mFalgn),4);
[uu,ll] = svd(Dst);
% if mod(order,2)~=0
%     Dst = (Dst+1)/2;
% end
[scclass,V,E]  = spectralcluster(Dst,20,'distance','precomputed');
jt = jet(8);
% [scclass8,V8,E8]  = spectralcluster(Dst8,12,'distance','precomputed');

[sccl,shoss] = clust_algn(full(knncl.mFalgn),scclass,order);
if mod(order,2)==0
    Dst2 = momdist(full(sccl.mFalgn),order);
    [scclass2,V2,E2]  = spectralcluster(Dst2,1,'distance','precomputed');
    [sccl2,shoss2] = clust_algn(full(sccl.mFalgn),scclass2,order);
    sccl.D = ifft(fft(sccl2.D(:,scclass)).*fft(sccl.D));
    sccl.Falgn = ifft(fft(sccl.Falgn).*fft(sccl2.D(:,scclass)));
    
    knncl.mFalgn = ifft(fft(full(sccl.D)).*fft(full(knncl.mFalgn)));
    knncl.D(:,knnclass>0) = ifft(fft(full(sccl.D(:,knnclass(knnclass>0)))).*fft(full(knncl.D(:,knnclass>0))));
    knncl.Falgn(:,knnclass>0) = ifft(fft(full(sccl.D(:,knnclass(knnclass>0)))).*fft(knncl.Falgn(:,knnclass>0)));
end

scclm2 = []; 
for k = 1:nclust
    scclm2(:,k) = mean(x(:,scclass(knnclass2)==k),2); %#ok<*SAGROW>
end
%% Mean PSDs
PSD = [afsrt.ftnorm];

%% Make plots
knncm = []; 
mPSD=[];
nclust = size(knncl.mFalgn,2);
for k = 1:nclust
    knncm(:,k) = mean(x(:,knnclass==k),2);
    mPSD(:,k) = nanmean(PSD(:,knnclass==k),2);
end
knncm(isnan(knncm))=0;
%%

dd = sqrt(sum((permute(uu(:,1:2),[1 3 2])-permute(uu(:,1:2),[3 1 2])).^2,3));
sc = median(min(dd + 1e9*(dd==0)))*5;
% r = sqrt(knncm(getu==2,:).^2 + knncm(getu==3,:).^2); 
% mang = mod(atan2(knncm(getu==2,:),knncm(getu==3,:)),2*pi);
r = sqrt((uu(:,1)-mean(uu(:,1))).^2 +( uu(:,2)-mean(uu(:,2))).^2)'; 
mang = mod(atan2(uu(:,1)-mean(uu(:,1)),uu(:,2)-mean(uu(:,2))),2*pi)';
rr = r./max(r);
cmap = hsv2rgb([mang'/(2*pi) .9*rr' 0*rr'] + [0 .1 1]);

%   cmap=jt(scclass,:);

nodex = uu(:,1:2)'/sc;
% nodex = knncm./sc;
% nodex = [cos(mang); sin(mang)].*sqrt(rr)*10;

plotxx = knncl.mFalgn;
dxx = 1.1*diff(tref([1 end]));
dyy = diff(minmax(plotxx(:)))*.5;

wlog=linspace(log10(.25),log10(160),sum(w>0&w<160));
plotx = interp1(w,abs(fft(full(knncl.mFalgn))),10.^wlog);
dx = 1.1*diff(wlog([1 end]));
dy = diff(minmax(plotx(:)))*1;

ang=mod(atan2(x(getu==2,:)'-mean(x(getu==2,:)'),x(getu==3,:)'-mean(x(getu==3,:)')),2*pi)/(2*pi);
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



