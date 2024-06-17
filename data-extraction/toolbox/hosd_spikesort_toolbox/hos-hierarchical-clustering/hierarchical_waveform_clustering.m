
function [classout,Xaligned,D,quality,Filt,net] = hierarchical_waveform_clustering(X,nclass,exvar,method)

% Recursive application of kmeans or SOM and skewness-based alignment

step = 5;
cluster_type = 'kmeans';
if nargin < 3 || isempty(exvar)
    exvar = 1;
end

if nargin < 4 || isempty(method)
        distancefun = @(x)momdist(x,3);
        
elseif isa(method,'function_handle')
        distancefun = method;
else
    if isnumeric(method)
        order = method;
    else

        switch method
            case {'skewness','bispect','bispectrum'}
                order = 3;
            case {'kurtosis','trispect','trispectrum'}
                order = 4;
            otherwise
                error('Urecognied method')
        end
    end

    distancefun = @(x)momdist(x,order);
end


    
X(isnan(X))=eps;
[Dst,Dsgn] = distancefun(X);
[u,l] = svd(Dst);


if nclass == 1
   class = ones(1,size(X,2));
else
    pckeep = find(diff([cumsum(diag(l))/sum(l(:))>=exvar;1])==1);
    
    switch cluster_type
        case 'kmeans'
            class = kmeans(u(:,1:pckeep),nclass,'distance','correlation','maxiter',1000);
            net=[];
        case 'som'
            net = selforgmap([floor(sqrt(nclass)) ceil(sqrt(nclass))], 100,3,'hextop','dist');
            [net,tr] = train(net,u(:,1:pckeep));
            y = net(u(:,1:pckeep));
            class = (1:size(y,1))*y;
    end
end

hos = hosobject(3,size(X,1),1,.3,[],[],'window',@rectwin);

fig = [];
    
    
    


D=[];
xsgn = [];
Filt = [];
fig=[nan nan];
for k = 1:max(class)
    fprintf('\nAligning %i of %i',k,max(class))
    x = X(:,class==k);
    if isempty(x)
        continue
    end
    if ~isempty(Dsgn)
       q = Dsgn(class==k,class==k);
       [qu,~] = svd(q);
       x = x*diag(sign(qu(:,1)));
       xsgn = sign(qu(:,1))';
    else
        xsgn=ones(1,size(x,2));
    end
    hos.get_block(x)  ;      
%     quality(class==k,:) = max(hos.xfilt(x))';
    
    [~,FXsh] = hos.apply_filter(x,false,true);
    Xsh = ifft(FXsh);
    fig(mod(k-1,2)+1)=gcf;
%         Xaligned(:,class==k) = Xsh;

    mX(:,k) = mean(Xsh,2);
    D(:,class==k) = (hos.sampt'==-hos.delay).*xsgn;
    Filt(:,class==k) = hos.filterfun*xsgn;
    if k>1 && ishandle(fig(mod(k,2)+1))
        close(fig(mod(k,2)+1))
    end
end

if mod(nclass,step)==0 && nclass>step

    [class2,~,D2,~,Filt2] = hierarchical_waveform_clustering(mX,nclass/step,exvar,method);

    DD2 = D2(:,class(:,1));

%     class(:,end+(1:size(class2,2))) = class2(class(:,1),:);
    FFilt2 =  Filt2(:,class,:).*sum(D);
    Filt = cat(3,FFilt2,Filt);
    classout = [class2(class,:),class];
%     quality = [quality2(class0,:),quality];
    for k = 1:nclass/step
        [~,~,unqi] =unique(class(class2(class)==k,1));
        classout(class2(class)==k,end) = unqi;
    end
    D = ifft(fft(D).*fft(DD2));
else
    classout=class;
end
Xaligned = ifft(fft(D).*fft(X));
% skfun = @(x,y)sum(ifft(abs(fft(x)).*fft(x).*fft(y)).^3)./sqrt(sum(ifft(abs(fft(x)).^2).^3).*sum(ifft((abs(fft(x)).*abs(fft(y))).^2).^3));
%  quality = squeeze(skfun(X,Filt));
skfun = @(x,y)sum(ifft(abs(fft(x)).*fft(x).*fft(y)).^3)./sqrt(sum(ifft(abs(fft(x)).^2).^3).*sum(ifft((abs(fft(x)).*abs(fft(y))).^2).^3));
 quality = squeeze(skfun(X,Filt));

% quality = squeeze(skewness(ifft(fft(X).*fft(Filt))));
% quality = squeeze(max(ifft(fft(Xaligned).*fft(Filt))));
        



function [Dst,Dsgn] = momdist(X,order)

XX = ifft(fft(X).*permute(conj(fft(X)),[1 3 2]));

if mod(order,2)==0
    mx = max(abs(XX));
    Dsgn = squeeze(sign(sum(XX.*(abs(XX)==mx))));
else
    Dsgn =[];
end
sK= sum(ifft(abs(fft(X)).^2).^order);
Dst = squeeze(sum(XX.^order))./sqrt(sK'*sK);
    
        
        
        