function [lag,t,dbal] = dbtpal(x,y,fs,varargin)

% DBT phase alignment.
%%
wbw = .25;  %2
tbw = .25; %.25


dbxw = dbt(x,fs,wbw,varargin{:});
dbyw = dbt(y,fs,wbw,varargin{:});
dbxt = dbt(x,fs,tbw,varargin{:});
dbyt = dbt(y,fs,tbw,varargin{:});

nfr = length(dbxw.frequency);
nt = length(dbxt.time);

phaw = @(lag)exp(-2*pi*1i*lag*dbxw.frequency);
phat = @(lag)exp(-2*pi*1i*lag*dbxt.frequency);
 
Lw = @(dbx,dby,H,lag,S) -1/2*sum(sum(abs(dby.blrep-dbx.blrep*diag(H).*phaw(lag)).^2*diag(1./S)))-sum(log(S));
Lt = @(dbx,dby,H,lag,S) -1/2*sum(sum(abs(dby.blrep-dbx.blrep*diag(H).*phat(lag)).^2*diag(1./S)))-sum(log(S));
lagt = 0*ones(nt,1);

tol = 1e-6;
k=1;
Hwold =0;
Hstep = Inf;
Lts = [];Lws=[];
while Hstep > tol && k < maxiter
    %%
    flagt = fft(lagt);
    flagt(ceil(length(dbxw.time)/2)+1:length(dbxw.time)) = 0;
    flagt(1) = .5*flagt(1);
    lagtw = 2*real(ifft(flagt));
    %%% Estimate transfer function 
    Hw = mean(dbyw.blrep.*conj(dbxw.blrep.*phaw(lagtw)))./(mean(abs(dbxw.blrep).^2) + eps);
    %%% Estimate noise
    Sw = mean(abs(dbyw.blrep-dbxw.blrep.*phaw(lagtw)*diag(Hw(1:nfr))).^2);
    Sw(isnan(Sw)) = 1e9;
    
    Hw(2*nfr+1) = 0;
    iHw = fft(Hw);
    iHw(end+1:length(dbxt.frequency)*2+1) = 0;
    Ht = ifft(iHw);
    Ht = Ht(1:length(dbxt.frequency));
    
    %%% Sinc interpolation of S
    iSw = fft(Sw);
    iSw(ceil(length(dbxt.frequency)/2+1):end) = 0;
    iSw(ceil(length(dbxt.frequency))+1:end) = [];
    iSw(ceil(length(dbxt.frequency))) = 0;
    iSw(1) = 1/2*iSw(1);
    St = 2*ifft(iSw)*length(iSw)./length(Sw);
    St = real(St(1:length(dbxt.frequency)));
    St(St<=0) = min(Sw);
    
    tol = 1e-6;
    dstep = Inf;
    maxiter=  500;

    i = 1;
    reg = 1e-4;
    Lold = -Inf;
%     optf = @(x) Lt(dbxt,dbyt,Ht,x,St);
%    options = optimoptions('fminunc','GradObj','on'); 
%     lagt = fminunc(@(x)lkfun(dbxt,dbyt,Ht,x,St),lagt,options);
    while dstep > tol && i < maxiter

        phat = exp(-2*pi*1i*lagt*dbxt.frequency);

        %%% First derivative wrt delta
        w = 2*pi*dbxt.frequency;
%        DL = real(nanmean(dbyt.blrep.*conj(dbxt.blrep.*phat)*diag(conj(Ht)./St.*1i.*w) - abs(dbxt.blrep).^2*diag(abs(Ht).^2.*1i.*w),2));
        DL = real(nanmean((dbyt.blrep.*conj(dbxt.blrep.*phat) - abs(dbxt.blrep).^2*diag(Ht))*diag(conj(Ht)./St.*1i.*w),2));

%        D2L = -real(nanmean((dbyt.blrep.*conj(dbxt.blrep.*phat))*diag(conj(Ht)./St.*w.^2),2));
        D2L = 0;

        step = -DL./(D2L+reg);
       
        L = Lt(dbxt,dbyt,Ht,lagt+step,St);
        if L>Lold
            reg = reg/2;
             lagt = lagt+step;
             dstep = sqrt(sum(step.^2)) + reg*1e-6
             Lold = L;
        else
            L = Lold
            reg = reg*2;
        end
        Ls(i,k) = L;
           i = i+1;
    end
    
    Lws(k) = Lw(dbxw,dbyw,Hw(1:nfr),lagtw,Sw);
    
    Hstep = sqrt(sum(abs(Hw-Hwold).^2))
    Hwold = Hw;
  
    k = k+1;
   
end


