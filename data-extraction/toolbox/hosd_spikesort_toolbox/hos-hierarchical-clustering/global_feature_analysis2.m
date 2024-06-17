inputdir = 'global_feature_anls_clust_PHASERANDOMIZED_out';
% inputdir = xne.subpaths.output.local;
d = dir(fullfile(inputdir,'allfeats*.mat'));

clear ld
for k = 1:length(d)
    ld(k) = load(fullfile(inputdir,d(k).name));
    k
end

allfeats = [ld.allfeats];
%% Compute pairwise bispectral correlation using the efficient technique of McGlaughlin 1968
af1 = allfeats([allfeats.component]==1);

mx = [af1.featfiltmax];
[srt,srti] = sort(mx,'descend');
afsrt = af1(srti);

compno = [afsrt.component];
F = [afsrt.feature];
F(isnan(F)) = 0;
PDF = [afsrt.pdfilt];
PDF(isnan(PDF)) = 0;

ftF = fft(F);
ftPDF = fft(PDF);
% FF = ifft(abs(ftF).^2);
FF = ifft(ftF.*ftPDF);
s3 = sum(FF.^3);
s4 = sum(FF.^4);
K3 = zeros(size(F,2));
K4 = zeros(size(F,2));
ftPDFperm = permute(ftPDF,[1 3 2]);
chunk = 100;
tic
for fk = 1:chunk:size(F,2)
    inds = fk:min(fk+chunk-1,size(F,2));
    ftFsub = ftF(:,inds);
%     krind = kron(ones(1,size(ftFsub,2)),1:size(ftPDF,2));
    
    fF =   ifft(ftFsub.*ftPDFperm);
%     Q = repmat(ftFsub,1,size(ftPDF,2)).*ftPDF(:,krind);
%     fF = ifft(ftF.*conj(ftF(:,fk)));
%     fF = ifft(ftF.*ftPDF(:,fk)); % This implicitly normalizes by SNR.
    K3(inds,:) = squeeze(sum(fF.^3))./sqrt(s3.*s3(inds)');
    K4(inds,:) = squeeze(sum(fF.^4))./sqrt(s4.*s4(inds)');
    fk
    toc
end
K=K3;
% KK = K(1:fk,:)*K(1:fk,:)';
KK = K3;

[u,l] = svd(KK);
v = diag(sqrt(diag(l)).^-1)*u'*K(1:fk,:);

L = linkage(v(1:10,:)','weighted','correlation');
class = cluster(L,'cutoff',.6);

Fft = [afsrt.ftnorm];
Fft(isnan(Fft))=0;
Cfrq = corr(Fft');
Cfrq(isnan(Cfrq))=0;

[ufrq,lfrq] = svd(Cfrq);
vfrq = diag(diag(lfrq).^-.5)*ufrq'*zscore(Fft);

Fnorm = ifft(fft(F).*abs(fft(PDF)));


%% Randomized SVD estimation
KK = K3*K3';
KKKK=KK*KK;
[u,l] = randsvd(KKKK,1000);
unorm = u./sqrt(sum(u.^2,2));
% [u,l] = svd(K3);
%% SOM classificaton
% x = vfrq(1:20,:);
x = u4(1:20,:);

net = selforgmap([20 20],100,3,'hextop','linkdist');

% Train the Network
[net,tr] = train(net,x);
y = net(x);
somclass = (1:size(y,1))*y;
MC = Fft*y'*(y*y'+eps)^-1;
MCx = x*y'*(y*y'+eps)^-1;
[srtc,srtic] = sort(somclass);

%% Data reduction through K means clustering
order = 4;
switch order
    case 4
        u=u4;
    case 3
        u=u3;
end

[kmclass,cm] = kmeans(u(1:50,:)',4096,'distance','correlation','maxiter',1000);

%% Fuzzy c means
[cent,U,obj_fcn] = fcm(u(1:50,:)',400);

idx = dbscan(u(1:20,:)',1,5);

%% Align clusters
hos = hosobject(3,length(tref),500,150,[],[],'window',@rectwin,'glowpass',150);
Falgn = zscore(Fnorm);
fig = [];
D=zeros(size(Fnorm));

class = somclass;
for k = 1:max(class)
    k
%     hos.reset();
    x = Falgn(:,class==k);
    if isempty(x)
        continue
    end
    if mod(order,2)==0
        ftx = fft(x);
        xx = ifft(ftx.*permute(conj(ftx),[1 3 2]));
        mx = max(abs(xx));
        xxsgn = squeeze(sum(xx.*(abs(xx)==mx)));
        [usgn,lsgn] = svd(xxsgn);
        xsgn = sign(usgn(:,1)');
        x = x.*xsgn;
    else
        xsgn = ones(1,size(x,2));
    end
    
    [A,B] = hos.get_block(x);
    fig(mod(k-1,2)+1)=gcf; %#ok<*SAGROW>
%     Falgn(:,class==k) = A;
    [~,FXsh] = hos.apply_filter(x,false,true);
    Xsh = ifft(FXsh);
     Falgn(:,class==k) = Xsh;
%         Xaligned(:,class==k) = Xsh;

    mX(:,k) = mean(Xsh,2);
    D(:,class==k) = (hos.sampt'==-hos.delay).*xsgn;
    if k>1 && ishandle(fig(mod(k,2)+1))
        close(fig(mod(k,2)+1))
    end
end

mFalgn = [];
mF = [];
cN = [];
for k = 1:max(class)
    mFalgn(:,k) = mean(Falgn(:,class==k),2);
    mF(:,k) = mean(F(:,class==k),2);
    cN(k) = sum(class==k);
end

%% 2nd level
[class2,mFalgn2,D,quality2,Filt2,net2] = hierarchical_waveform_clustering(mFalgn,50,.9,4);

%% plot SOM
nodex = net.layers{1}.positions;
dx = 1.1*diff(tref([1 end]));
plotx=zscore(mFalgn);
dy = .5*diff(minmax(plotx(:)));
figure
plot(tref'+dx*nodex(1,:),dy*nodex(2,:)+plotx,'k');

%% Hierarchical 4th order clustering
% oder = 3;
[mclass,mFalgn2,D2,quality,Filt] = hierarchical_waveform_clustering(mFalgn,64,.9,4);
mquality = zeros(size(mclass));
for k =1:size(mclass,2)
    [unq,~,unqi] = unique(mclass(:,k));
    for kk = 1:length(unq) 
        
        mquality(unqi==kk,k) = mean(quality(unqi==kk,k));
        classn(unqi==kk,k) = sum(unqi==kk);
    end
end
[unq,~,unqi] = unique(mclass,'rows');

class2 = zeros(size(class)+[0 size(mclass,2)]);

for k = 1:size(unq,1)
    csubi = unqi(class) == k;
    [~,~,unqi2] = unique(class(csubi));
    class2(csubi,1:end-1) = repmat(unq(k,:),sum(csubi),1);
    class2(csubi,end) = unqi2;
end    

ftD2 = fft(D2);
ftDelay = ftD2(:,class).*fft(D);
Delay = ifft(ftDelay);
ftFalgn2 = fft(zscore(Fnorm)).*ftDelay;
Falgn2 = ifft(ftFalgn2);

% Falgn2 = ifft(fft(Falgn).*ftD(:,class));

%% Weighted Phase 2nd deivative
f = fft(ifftshift(Falgn,1));
ff = f(2:round(end/2),:).*conj(f(1:round(end/2)-1,:));
ff2 = ff(2:round(end/2),:).*conj(ff(1:round(end/2)-1,:));
phslope = angle(sum(ff2)./sum(abs(ff2)))/pi;

% phslope = (w(1:end/2-1)*(w(1:end/2-1)'.*diag(abs(ff))))^-1*(w(1:end/2-1))*ff; %Linear slope of group delay vs frequency, weighted by power
%%
x2 = MCx;
net2 = selforgmap([4 4]);
% Train the Network
[net2,tr2] = train(net2,x2);
y2 = net2(x2);
c2 = (1:size(y2,1))*y2;

%% Atlas classifer

atl = {afsrt.BN246};
[unq,~,unqi] = unique(atl);
atlind = zeros(length(unq),size(atl,2));
atlind(unqi'+(0:size(atl,2)-1)*length(unq))=1;

proto = {afsrt.protocol};
[unqpr,~,unqipr] = unique(proto);
protind = zeros(length(unqpr),size(proto,2));
protind(unqipr'+(0:size(proto,2)-1)*length(unqpr))=1;

blockatl = [atlind;protind];

x = vfrq(1:10,:);
