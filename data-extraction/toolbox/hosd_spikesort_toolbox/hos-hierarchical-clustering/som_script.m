
getu = [2 3];
x = u4(getu,:);

net = selforgmap([20 20],100,3,'hextop','linkdist');

[net,tr] = train(net,x);

y = net(x);

somclassn = sum(y,2);


somclass = (1:size(y,1))*y;

[somcl,hoss] = clust_algn(Falgn,somclass,4);

somcm = []; 
for k = 1:400
    somcm(:,k) = mean(x(:,somclass==k),2);
end
somcm(isnan(somcm))=0;

cmap = hsv2rgb([(atan2(somcm(getu==2,:),somcm(getu==3,:))+pi)'/(2*pi) .9*rr' 0*rr'] + [0 .1 1]);

nodex = net.layers{1}.positions;

plotxx = somcl.mFalgn;
dyy = diff(minmax(plotxx(:)))*.25;

wlog=linspace(log10(.25),log10(160),sum(w>0&w<160));
plotx = interp1(w,abs(fft(full(somcl.mFalgn))),10.^wlog);
dy = diff(minmax(plotx(:)))*.25;

hsv = [atan2(x(getu==2,:)',x(getu==3,:)')/(2*pi)+.5,.1+.9*sqrt(sum(x(getu==2|getu==3,:).^2))'/max(r), ones(size(u4,2),1)];
hsv(hsv>1)=1;

cols = hsv2rgb(hsv);
figure, 
mscatter3(x(getu==2|getu==3|getu==4,:)',1,cols,'markerfaceColor','flat')
axis equal vis3d off
set(gcf,'color','black')

figure, 
hold on, 
for k = 1:400
    h = plot(tref'+dx*nodex(1,k),dyy*nodex(2,k) +plotxx(:,k),'color',cmap(k,:));
end
axis off tight
set(gcf,'color','black')

figure, 
hold on, 
for k = 1:400
    h = plot(wlog'+dx*nodex(1,k),dy*nodex(2,k) +plotx(:,k),'color',cmap(k,:));
end
axis off tight
set(gcf,'color','black')
plh = plot(wlog'+dx*nodex(1,:),dy*nodex(2,:) +.5*dy*ismember(wlog',wlog([139  365  594  ])),'color',[1 1 1]*.5);

