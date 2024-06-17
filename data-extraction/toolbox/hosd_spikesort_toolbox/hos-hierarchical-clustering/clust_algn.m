function [out,outhoss] = clust_algn(Fnorm,class,order,do_plot,sampling_rate)

if nargin < 4 || isempty(do_plot)
    do_plot=false;
end
if nargin < 5 || isempty(sampling_rate)
    sampling_rate = 500;
end

hos = hosobject(3,size(Fnorm,1),sampling_rate,200,[],[],'window',@rectwin,'glowpass',150);
Fz = zscore(Fnorm);
fig = [];
Falgn=zeros(size(Fnorm));
D=zeros(size(Fnorm));

mFalgn = sparse(zeros(size(Falgn,1),max(class)));
mF = mFalgn;
cN = [];
fig = [nan nan];
runi = 1:max(class);
for ki = 1:length(runi)
    k =runi(ki)
%     hos.reset();
    x = Fz(:,class==k);
  
    if isempty(x)
        cN(k) = 0;
        
        mFalgn(:,k)=0;
        continue
    elseif   size(x,2) == 1
         D(1,class==k) = 1;
          mFalgn(:,k) = x; %#ok<*SPRIX>
          cN(k) = 1;
          continue
    end
    if mod(order,2)==0
       ftx = fft(x);
        xxsgn = zeros(size(ftx,2));
        for kk =  1:size(ftx,2)
            xx = ifft(ftx.*conj(ftx(:,kk)));
            mx = max(abs(xx));
            xxsgn(kk,:) = squeeze(sum(xx.*(abs(xx)==mx)));
        end
        [usgn,lsgn] = svd(xxsgn);
        xsgn = sign(usgn(:,1)');
          xsgn = xsgn*sign(mean(usgn(:,1))+eps);
        x = x.*xsgn;
    else
        xsgn = ones(1,size(x,2));
    end
    
    [A,B] = hos.get_block(x,25,do_plot);
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

    mFalgn(:,k) = mean(Falgn(:,class==k),2);
%     mF(:,k) = mean(F(:,class==k),2);
    cN(k) = sum(class==k);
    if nargout > 1
        outhoss(k,:) = hosobject(hos);
    end
    quality(k,:) = max(hos.xfilt(hos.feature));
end
out.Falgn = Falgn;
out.D = D;
out.mFalgn = mFalgn;
% out.mF = mF;
out.cN = cN;
out.quality=quality;
