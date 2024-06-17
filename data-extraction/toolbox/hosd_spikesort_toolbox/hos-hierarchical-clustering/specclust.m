
K = (K4+K4')/2;
%%



% %% Initial reduction with k means
% 

% x = diag(diag(l4(1:50,1:50)).^.25)*u4(1:50,:);
x = u4(1:50,:);
 [kmclass,cm] = kmeans(x',2048,'maxiter',1000,'replicates',5);

for k = 1:2048
    cm(k,:) = mean(x(:,kmclass==k),2);
end

order = 4;


%% Align clusters
hos = hosobject(3,length(tref),500,150,[],[],'window',@rectwin,'glowpass',150);
Falgn = zscore(Fnorm);
fig = [];
D=zeros(size(Fnorm));

% class = somclass;
for k = 1:max(class)
    k
%     hos.reset();
    xf = Falgn(:,class==k);
    if isempty(xf)
        continue
    end
    if mod(order,2)==0
        ftx = fft(xf);
        xx = ifft(ftx.*permute(conj(ftx),[1 3 2]));
        mx = max(abs(xx));
        xxsgn = squeeze(sum(xx.*(abs(xx)==mx)));
        [usgn,lsgn] = svd(xxsgn);
        xsgn = sign(usgn(:,1)');
        xf = xf.*xsgn;
    else
        xsgn = ones(1,size(xf,2));
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
mx= [];
for k = 1:max(class)
    mFalgn(:,k) = mean(Falgn(:,class==k),2);
    mF(:,k) = mean(F(:,class==k),2);
    cN(k) = sum(class==k);
    mx(:,k) = mean(x(:,class==k),2);
end
toc

Dst = momdist(mFalgn,4);
Dst = (Dst+Dst')/2;
[u2,l2] = svd(Dst);

tic
[scclass,V,E]  = spectralcluster(Dst,256,'distance','precomputed');
toc
sccl2 = clust_algn(mFalgn,scclass,4);
sccl3 = clust_algn(Falgn,scclass(kmclass),4);

%%
x2 = u2(:,1:50)';

net = selforgmap([16 16],100,3,'hextop','linkdist');
[net,tr] = train(net,x2);
y = net(x2);
somclass = (1:size(y,1))*y;
somclassn = sum(y,2);
for k = 1:2048
    somcm(k,:) = mean(x2(:,somclass==k),2);
end

[somcl,hoss] = clust_algn(Falgn,somclass,4);
for k = 1:400
    somcm(:,k) = mean(u4(:,somclass==k),2);
end


Dst = momdist(full(somcl.mFalgn),4);
Dst = (Dst+Dst')/2;
[mscclass,mV,mE]  = spectralcluster(Dst,10,'distance','precomputed');
