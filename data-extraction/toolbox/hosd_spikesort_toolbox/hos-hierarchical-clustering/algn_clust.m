

if exist('jobindex','var')
    load(fullfile(inputdir,'algndata'))
    runi = jobindex;
    savefname = fullfile(outputdir,sprintf('out%04i.mat',jobindex));
    if exist(savefname,'file')
        return
    end
else
    runi = 1:max(class);
    savefname = '';
end

hos = hosobject(3,length(tref),500,200,[],[],'window',@rectwin,'glowpass',200);
Fz = zscore(Fnorm);
fig = [];
Falgn=sparse(zeros(size(Fnorm)));
D=sparse(zeros(size(Fnorm)));

mFalgn = sparse(zeros(size(Falgn,1),max(class)));
mF = mFalgn;
cN = [];
fig = [nan nan];
% class = somclass;
for ki = 1:length(runi)
    k =runi(ki);
%     hos.reset();
    x = Fz(:,class==k);
    if isempty(x)
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

    mFalgn(:,k) = mean(Falgn(:,class==k),2);
%     mF(:,k) = mean(F(:,class==k),2);
    cN(k) = sum(class==k);
end
out.Falgn = Falgn;
out.D = D;
out.mFalgn = mFalgn;
% out.mF = mF;
out.cN = cN;
out.jobindex = jobindex;
if ~isempty(savefname)
    save(savefname,'-struct','out')
end
