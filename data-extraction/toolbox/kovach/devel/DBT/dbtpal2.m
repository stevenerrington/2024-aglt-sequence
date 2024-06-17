function [lag,t,dbal] = dbtpal2(x,y,fs,varargin)

% DBT phase alignment.
%%
wbw = 2;  %2
tbw = wbw;
dbxw = dbt(x,fs,wbw,varargin{:});
dbyw = dbt(y,fs,wbw,varargin{:});
dbxt = dbt(x,fs,tbw,varargin{:});
dbyt = dbt(y,fs,tbw,varargin{:});

nfr = length(dbxw.frequency);
nt = length(dbxt.time);

lagt = 0*ones(nt,1);

w = 2*pi*dbxw.frequency;
phaw = @(lag)exp(-1i*lag*w);
%phat = @(lag)exp(-2*pi*1i*lag*dbxt.frequency);

nt = length(dbxw.time);

%%% Phase adjusted X
phaX = @(lag)dbxw.blrep.*phaw(lag);

%%% Y estimate
Yest = @(lag,H) phaX(lag)*diag(H);

%%% Residual
R = @(lag,H)dbyw.blrep-Yest(lag,H);

%%% Sigma estimate
S =@(lag,H) mean(abs(R(lag,H)).^2)/2;

%%% Log lihelihood
LL = @(lag,H) -.5*sum(sum(abs(R(lag,H)).^2)./S(lag,H) + 2*nt*log(S(lag,H)));

%%% Derivative with respect to H
DLDH = @(lag,H) (sum(R(lag,H).*conj(phaX(lag)))*diag(1./S(lag,H)));

%%% Derivative with respect to lag
DLDT = @(lag,H) -sum( imag(dbyw.blrep.*conj(Yest(lag,H))) * diag(w./S(lag,H)),2);


%Second derivatives 
D2LDH2 = -sum(abs(dbxw.blrep).^2)*diag(1./S(lag,H));
D2LDT2 = @(lag,H) -sum(real(dbyw.blrep.*conj(Yest(lag,H)))* diag(w.^2./S(lag,H)),2);
D2LDTDH = @(lag,H) -imag(dbyw.blrep.*conj(phaX(lag)))* diag(w./S(lag,H));

%Put 1st derivatives into a column vector
DL = @(lag,H) cat(1,conj(DLDH(lag,H)'),DLDT(lag,H));

%Construct the full Hessian matrix 
DL2 = @(lag,H) conj([ diag(D2LDH2)   , D2LDTDH(lag,H)' 
                 D2LDTDH(lag,H) , diag(D2LDT2(lag,H)) ]);

flagt = fft(lagt);
flagt(ceil(length(dbxw.time)/2)+1:length(dbxw.time)) = 0;
flagt(1) = .5*flagt(1);
lagtw = 2*real(ifft(flagt));
Hw =conj( mean(dbyw.blrep.*conj(dbxw.blrep.*phaw(lagtw)))./(mean(abs(dbxw.blrep).^2) + eps))';
%Hw = complexglm(dbxw.blrep,dbyw.blrep);


tol = 1e-6;
k=1;
Hwold =0;
Hstep = Inf;
Lts = [];Lws=[];
LLs= [];
maxiter=1e5;
ll = @(x)LL(x(nfr+(1:nt)),conj(x(1:nfr)'));
%%
reg  = .001;
LLs = [];
while Hstep + reg > tol  && k < maxiter
    %%
    LLs(k) = LL(lagtw,Hw);
    d1 = DL(lagtw,Hw);
    d2 = DL2(lagtw,Hw);
    step = (d2+ eye(size(d2))*reg)\d1;
    % step(real(step.*d1)<0) = -step(real(step.*d1)<0);
    if LL(lagtw+step(nfr+(1:nt)),Hw+step(1:nfr))>=LLs(k)
        Hw = Hw+step(1:nfr);
        lagtw = lagtw+ real(step(nfr+(1:nt)));
        reg = reg/2;
    elseif LL(lagtw-real(step(nfr+(1:nt))),Hw-step(1:nfr))>=LLs(k)
        Hw = Hw-step(1:nfr);
        lagtw = lagtw- real(step(nfr+(1:nt)));
        reg = reg/2;
    else
        reg = reg*2;
    end
    Hstep = sqrt(sum(abs(step).^2));
    k = k+1;
   fprintf('\nreg: %0.4g, step: %0.4g',reg,Hstep)
end


