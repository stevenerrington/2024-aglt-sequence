function [P,stats] = bicpolyreg(Y,X,varargin)


%Stepwise polynomial regression to minimize AIC or BIC
% 
% PX = X;
% Currently Y must be univariate. 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/bicpolyreg.m $
% $Revision: 185 $
% $Date: 2013-04-02 16:12:00 -0500 (Tue, 02 Apr 2013) $
% $Author: ckovach $
% ------------------------------------------------

criterion ='aic';
maxord = 10;
regargs = {};
i = 1;
minord = 0;
while i <= length(varargin)
    
    switch lower(varargin{i})
        
        case 'criterion' 
            criterion = varargin{i+1};
            i = i+1;
        case 'maxord' 
            maxord = varargin{i+1};
            i = i+1;
        case 'minord' 
            minord = varargin{i+1};
            i = i+1;
        case 'order' 
            minord = varargin{i+1};
            maxord = varargin{i+1};
            i = i+1;
        case 'regargs' 
            regargs = cat(2,regargs,varargin(i+1));
            i = i+1;
    
    
     otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
        i = i+1;
end


PXs = [];
AIC = [];
AICc = [];
BIC = [];
if isempty(regargs) || ~iscell(regargs{1})
    regargs = repmat({regargs},1,size(X,2));
end

ords = minord:maxord;
    
for i = 1:length(ords)
    
%    ps = makeregressor(ones(size(X,1),1)); %dc component
      ps = buildpolyreg(X,ords(i),'dc'); %dc component
%       if ords(i) > 1
%         for k = 1:size(X,2)
% 
%             ps(k+1) = buildpolyreg(X(:,k),ords(i)-1,regargs{k}{:},'codeincr',k);
%         end
% %     ps = pool(cat(2,ps,interaction(ps(2:end))));
% %        ps = pool(interaction(ps));
%          ps = pool(ps(1),interaction(ps));
%       end
    PX = ps.value;
%     PX = cat(2,ones(size(PX,1),1),PX);
    P = (PX'*PX)\(PX'*Y);
    
    RSS = (Y-PX*P)'*(Y-PX*P);
    
    k = size(PX,2);
    p = size(RSS,2);
    N = size(Y,1);
    D = N*log(det(RSS)./N) + N*p*log(2*pi+1); %Deviance
    AIC(i) = D + 2*(k*p + p*(p+1)/2);        %AIC
    AICc(i) = D + N*p*log(2*pi./(2*pi+1)) + (N+k)*N*p./(N-k-p-1); % Bias corrected AIC
    
    BIC(i) = D + (p*k + p*(p+1)/2)*log(N);  %BIC
    
    Ps(1:size(P),i) = P;
    Ks(i) = k;
    Ds(i) = D;
end

if isnumeric(criterion)
    mni = criterion+1;
    P = Ps(1:Ks(mni),mni);    
elseif strcmpi(criterion,'aicc')
    [mn,mni] = min(AICc);
    P = Ps(1:Ks(mni),mni);
elseif strcmpi(criterion,'aic')
    [mn,mni] = min(AIC);
    P = Ps(1:Ks(mni),mni);
elseif strcmpi(criterion,'bic')
    [mn,mni] = min(BIC);
    P = Ps(1:Ks(mni),mni);
end

% ps = buildpolyreg(X,mni-1,regargs{:});

  ps = buildpolyreg(X,ords(mni),'dc'); %dc component
%   ps = buildpolyreg(X(:,1),0,'dc'); %dc component
% %   if mni-1 > 0
%     for k = 1:size(X,2)
%      
%             ps(k+1) = buildpolyreg(X(:,k),ords(mni),regargs{k}{:});
%     end
% %   end
% 
%   ps = pool(ps(1),interaction(ps));
% %   ps = pool(cat(2,ps,interaction(ps(2:end))));

PX = ps.value;

RSS = (Y-PX*P)'*(Y-PX*P);
   
FI = N*(PX'*PX)/RSS; %Fisher's information matrix for the means

% E = (PX*inv(FI)*PX');
if size(X,2) > 1
    errfun = @(varargin) sqrt(sum(ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})'.*(inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})'),1))'; %Returns the ML estimation error for the given contrast
%     errfun = @(varargin) sqrt(diag(ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})*inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})')); %Returns the ML estimation error for the given contrast
    estfun = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})*P; %Returns the estimation given contrast
    errmat = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})*inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:})'; %Returns the ML estimation error for the given contrast
else
    errfun = @(varargin) sqrt(sum(ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})'.*(inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})'),1))'; %Returns the ML estimation error for the given contrast
%     errfun = @(varargin) sqrt(diag(ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})*inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})')); %Returns the ML estimation error for the given contrast
    estfun = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})*P; %Returns the estimation given contrast
    errmat = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})*inv(FI)*ps.function(varargin{1},varargin{kron(1:nargin,[1 1])})'; %Returns the ML estimation error for the given contrast
end
stats.Ps = PXs;
stats.D = Ds;
stats.AIC = AIC;
stats.AICc = AICc;
stats.BIC = BIC;
stats.beta = P;
stats.FI = FI;
stats.resvar = RSS./N;
stats.estfun = estfun;
stats.errfun = errfun;
stats.errmatfun = errmat;
% stats.parfun = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])},varargin{:});
stats.parfun = @(varargin) ps.function(varargin{1},varargin{kron(1:nargin,[1 1])});
% stats.parfun =  ps.function;
