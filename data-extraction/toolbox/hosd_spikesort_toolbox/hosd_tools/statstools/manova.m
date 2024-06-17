
function st = manova(X,factor,varargin)

% st = manova(X,factor,varargin)
% Implements multivariate anova based on cell means model,  following pca dimensionality reduction.
% Does not require balanced data. All cells must contain at least 1
% element. Columns of X are observations. Each column of FACTOR is a factor and each unique 
% value within each column corresponds to a level within the factor. i.e.
%
%               FACTOR = [ 1  0  
%                          3  0        represents 2 factors, the first has 3
%                          2  8  ];    levels and the second has 2 levels.
%
% The number of columns of X must equal the number of rows in FACTOR.
%
% Use options 'exvar' to retain a certain amount of variance:
%
%	st = manova(X,factor,'exvar',.9)
%
% Keeps the number of principal components required to explain 90% of the variance.
%
% Use options 'pckeep' to retain a certain number of principle components:
%
%	st = manova(X,factor,pckeep',10)
%
% Keeps the first 10 principal components ordered by variance.
%
% st is a structure with several fields
%
% st.wksp is  p value using Wilk's statistic.
% st.WKS  is a structure with more information on Wilk's test.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



% C. Kovach 2003
%
% christopher-kovach@uiowa.edu
%
%

exvar = .95;
pckeep = [];
discrim = 0;
default_trim = .1;
drop = 0;
trim = [];
ridge_var = [];
ngroups = 2;
Data = {};
filt = [];

pillai = 0;
wilks = 1;
lawhot = 0;
roy = 0;
rootTest = 0; %Keep or drop principle components according to their relative variance explained.
rootXplain = 1; %Threshold in multiples of the average variance.
center = 0;
correl = 0;
i = 1;
compact = 0;
randeff = 0;
interaction_level = min([size(factor,2),2]);
factor_labels = 1:size(factor,2);


varargin{end+1} = 'finis';

while i <= length(varargin)
    switch lower(varargin{i})
       case 'exvar'
            exvar = varargin{i+1};
            i = i+1;
       case 'pckeep'
           pckeep = varargin{i+1};
           i = i+1;
       case 'dropfirst'
           drop =  varargin{i+1};
           i = i+1;
       case 'discrim'
           discrim = 1;
       case 'covar'
           center = 1;
       case 'correl'
            correl = 1;
       case 'trim'
           trim = default_trim;
         if i <length(varargin)
           if ~ischar(varargin{i+1})
               trim = varargin{i+1};
               i = i+1;
           end
         end
       case 'ridge'
           ridge_var = varargin{i+1};
           i = i+1;
       case 'pillai'
           pillai = 1;
       case 'level'
           interaction_level = varargin{i+1};
           i = i+1;
       case 'pillai'
           pillai = 1;
       case '~pillai'
           pillai = 0;
       case 'wilks'
           wilks = 1;
       case '~wilks'
           wilks = 0;
       case 'lawhot'
           lawhot = 1;
       case '~lawhot'
           lawhot = 0;
%        case 'roy'
%            roy = 1;
%        case '~roy'
%            roy = 0;
       case 'exclude'
           exclude = varargin{i+1};
           i = i+1;
           
           excl = find(factor(:,exclude(1)) == exclude(2));
           X(:,excl) = [];
           factor(excl,:) = [];
  
       case 'labels'
           factor_labels = varargin{i+1};
           i = i+1;
       case 'center'
           center = 1;
        case 'filter'
           filt = varargin{i+1};
           i = i+1;
       case 'randeff'
           randeff = 1;
       case 'compact'
           compact = 1;   
       case 'roottest'
           rootTest = 1;
           
           if ~ischar(varargin{i+1}) 
               rootXplain = varargin{i+1};
               i = i+1;
           end
        case 'finis'
       otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end

if interaction_level > size(factor,2)
    interaction_level = min([size(factor,2),2]);
end

if size(factor,1) ~= size(X,2)    
    error('Number of trials and size of factor matrix do not match.')
end


if ~isempty(filt)
    X = bandpassfilter(X,filt(1),filt(2),filt(3));   
end

if ~isempty(trim)
  [X,dsc] = multitrim(X,trim);
  factor(dsc,:) = [];
else
    trim = 0;
end




nxtrials  = size(X,2);


mx = mean(X,2);

totmean = mx;

if center
    X = X - totmean*ones(1,size(X,2));      
end

flip = 0;

if size(X,1) >= size(X,2) 
  S = X'*X;
  flip = 1;
else
  S = X*X';
end




if center & flip
    mp = mean(X,1);
    S = S - mp'*mp;         %Covariance matrix
elseif center
    mp = mean(X,2);
    S = S - mp*mp';
end

if correl
    DS = diag(S).^.5;
    S = diag(DS.^-1)*S*diag(DS.^-1); %correlation matrix
end



[D,L] = svd(S); %Singular value decomposition of S. 

l = diag(L);

if rootTest
    if ~isempty(pckeep)
        warning('specified pckeep overidden by root test')
    end
    propmeanl = l(drop+1:end)/mean(l);    
    [dpmin,pckeep] = min((propmeanl - rootXplain).^2);
end

if isempty(pckeep)
   pckeep = min(find(cumsum(l(drop+1:end)/sum(l))>=exvar))+drop;
    if isempty(pckeep)
        pckeep = rank(S);
    end
else
    exvar = sum(l((drop+1):pckeep))/sum(l);
end


if pckeep > rank(S)-ngroups & size(X,1)>size(X,2)+ngroups
    pckeep = rank(S) - ngroups;
end
    
if flip

    D = X*D*diag(diag(L).^-.5);        %Converts eigenvectors of X'X to nonzero eigenvectors of XX'
end

D = D(:,(1+drop):pckeep);
pckeep = pckeep - drop;

X = D'*X;  %X replaced by pca scores.



unifact = unique(factor,'rows');
ncells = size(unifact,1);
nfactors = size(factor,2);

W = [];
cellmemb = [];

cellN = [];
for i = 1:ncells
    
    cellpos = ismember(factor,unifact(i,:),'rows');
    
    cellN(i) = sum(cellpos); %Number of trials in cell
    
    celltrials = find(cellpos);  %Trials in cell
    
    
    W(celltrials,i) = 1;

    cellmemb(celltrials) = i;    %Cell membership of each trial.
    
    %nM(:,i) = sum(X(:,2),2);    
end


M = X*W*(W'*W)^-1;

E = (X' - W*M')'*(X' - W*M');  %error matrix (within group variance). 
invE = E^-1;

sizem = size(M,1);
%Constructing contrasts

Contrasts = {};
DFs = [];

for i = 1:nfactors;
    
    levels = unique(factor(:,i));
    nlevs(i) = length(levels);

    c = eye(nlevs(i));
    
    c = c - (c*ones(nlevs(i),1))*ones(1,nlevs(i))*1/nlevs(i);
    
    c = c(1:(nlevs(i)-1),:);
    
    redfact = unifact;
    redfact(:,i) = [];
    redunifact = unique(redfact,'rows');                  
    for j = 1:size(redunifact,1)                           
        Contrasts{i}(:,find(ismember(redfact,redunifact(j,:),'rows'))) = c;
    end
    
    %reps = ncells/nlevs(i);    
    %Contrasts{i} = kron(ones(1,reps),c);
    DFs(i) = size(c,1);
end

ve = size(X,2) - prod(nlevs);

index = {};

if wilks
 st.WKS = {};
 st.wksps = [];
end
if pillai
 st.PIL = {};
 st.pilps = [];
end
if roy
  st.ROY = {};
  st.royps = [];
end
if lawhot
  st.LH = {};
  st.lhps = [];
end
if discrim
    discr = {};
    stdiscr = {};
    pcadiscr = {};
end

 J = 1:interaction_level;

if randeff
    J = [interaction_level,J(1:end-1)];
end

for j = J
  
  combinations = factcombs(nfactors,j);  
  for i = 1:nchoosek(nfactors,j)
      
      dfs = DFs(combinations(i,:));
      convert = fliplr(cumprod(fliplr(dfs)));
      convert = [convert(2:end),1];
      
      
      C = [];
      
      for k = 1:prod(dfs)
          
          contracomb = mod(floor(k./convert),dfs) + 1;
          
          Q = ones(1,size(Contrasts{1},2));
          
          for s = 1:length(combinations(i,:))
            Q = Q.*Contrasts{combinations(i,s)}(contracomb(s),:);
          end
          
          C(k,:) = Q;
          
       end
          
          
      
      H = (C*M')'*[C*(W'*W)^-1*C']^-1*(C*M');
      vh = size(C,1);      
      
      if randeff & j == interaction_level  %In random effects model, Error matrix given by SS of interaction.
      
          E = H;
          ve = vh;
          
      else
          
        if wilks
          st.WKS{end+1} = wilksstat(E,H,ve,vh,sizem);
          st.wksps(end+1) = st.WKS{end}.pval;
        end
      
        if pillai
          st.PIL{end+1} = pillaistat(E,H,ve,vh,sizem);
          st.pilps(end+1) = st.PIL{end}.pval;
        end
        
        if lawhot
          st.LH{end+1} = lawhotstat(E,H,ve,vh,sizem);
          st.lhps(end+1) = st.LH{end}.pval;
        end

        if roy
          st.ROY{end+1} = roysstat(E,H,ve,vh,sizem)
        end
      
         index{end+1} = combinations(i,:);
        if discrim
          
          [dscr,d] = svd(invE*H);
          discr{end+1} = D*dscr(:,1:rank(d));
          pcadiscr{end+1} = dscr(:,1:rank(d));
        end
    end
  end
end
st.factormx = factor;     
st.factlist = unifact;
st.factor_labels = factor_labels;
st.index = index;
st.cellN = cellN;
st.cellmemb = cellmemb;
st.pckeep = pckeep;
st.exvar = exvar;
st.scree = l./mean(l);
if discrim
  st.discr = discr;
  st.pcadiscr = pcadiscr;
end

if ~compact
    st.pcascores = X;
    st.D = D;
    st.S = S;
end
if rootTest
    st.rootTest = rootXplain;
end
st.trim = trim;
if trim
    st.discarded = dsc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function wk = wilksstat(E,H,ve,vh,sizem)
     
    condn = 1./(trace(E)./rank(E));; %In case det(E) evaluates to Inf, the numerator and denominator below are divided by this value
    Lambda = det(E/condn)/(det((E+H)/condn));
   % Lambda = det(E)/(det(E+H));
          
    w = ve + vh - 1/2*(sizem + vh + 1);
          
    if sizem*vh == 2
        t = 1;
    else
        t = sqrt((sizem^2*vh^2 - 4)/(sizem^2 + vh^2 - 5));
    end
          
    df1 = sizem*vh;
    df2 = w*t - 1/2*(sizem*vh - 2);
    wk.F = df2/df1*(1 - Lambda^(1/t))/Lambda^(1/t);  
    wk.pval = 1-fcdf(real(wk.F),df1,df2);
    wk.df1 = df1;
    wk.df2 = df2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ry = roysstat(E,H,ve,vh,sizem);
%unfinished
    L = svd(E^-1*H);
    
    theta = L(1)/(L(1)+1);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pl = pillaistat(E,H,ve,vh,sizem)

s = min([vh,sizem]);
d = max([vh,sizem]);
V = trace((E+H)^-1*H);

pl.F =  ((ve - sizem + s)*V)/(d*(s - V));
pl.df1 = s*d;
pl.df2 = s*(ve - sizem + s);
pl.pval = 1-fcdf(real(pl.F),pl.df1,pl.df2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lh = lawhotstat(E,H,ve,vh,sizem)

U = trace(E^-1*H);

a = sizem*vh;
B = (ve+vh-sizem-1)*(ve-1)/((ve-sizem-3)*(ve-sizem));
b = 4 + (a+2)/(B - 1);
c = a*(b-2)/(b*(ve-sizem-1));

lh.F = U/c;
lh.df1 = a;
lh.df2 = b;
lh.pval = 1-fcdf(real(lh.F),lh.df1,lh.df2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,disc] = multitrim(X,trim,exvar);

% [X,Disc] = multitrim(X);
%Removes outliers from X using Multivariate iterative trimming of extreme data points. See Rencher p. 29.

if nargin < 2
    trim = .2;
end

if nargin < 3
  exvar = .95;
end

X = squeeze(X);



Xc = X;


mx = mean(Xc,2);
n = size(Xc,2);
%MX = mx(:,ones(size(Y,2),1));

Xc = Xc-mx*ones(1,size(Xc,2));

if size(X,1) > size(X,2) 
    S = (Xc'*Xc)/(n-1);   %Because S is singular, reduce dimensionality with svd.

    [U,D] = svd(S);

    d = diag(D);

    U = Xc*U*diag((d).^-.5);        %Converts eigenvectors of X'X to nonzero eigenvectors of XX'

    pckeep = min(find(cumsum(d/sum(d))>=exvar));

    r = rank(S);
    if pckeep > r | isempty(pckeep);
        pckeep = r;
    end

    d = d(1:pckeep);
    U = U(:,1:pckeep);

    Xcu = U'*Xc;
else
    Xcu = Xc;
end

zold = inf;
Disc = [];

Xu = Xcu;  
while 1

 Su = Xcu*Xcu';
 R = diag(diag(Su).^-.5)*Su*diag(diag(Su).^-.5);
 rs = R(find(triu(R+1e-15,1)));
 
 znew = .5*log((1+rs)./(1-rs));  % Fisher's z transform
  
 if mean(abs(znew-zold)) < .001
     break
 end
 
 zold = znew;
% T = diag(Xcu'*Su^-1*Xcu); 
  T = diag(Xu'*Su^-1*Xu); 

 [tsort,tind] = sort(T);
 disc = tind(end - ceil(trim*length(tind)) + 1:end);
 %Disc = cat(1,Disc,disc);
 %Xcu(:,disc) = [];
 Xcu = Xu; 
 Xcu(:,disc) = [];
 Xcu = Xcu - mean(Xcu,2)*ones(1,size(Xcu,2));
 if rank(Su) < max(size(Su));
     warning('Matrix is singular.')
     break
 end
end

X(:,disc) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = factcombs(nfacts,level)



P = [];


P(1,:) = 1:level;

index = level;
i = 2;

while i <= nchoosek(nfacts,level) 
               
            
        if P(i-1,index) < nfacts-level + index
            P(i,:) = P(i-1,:);
            P(i,index) = P(i-1,index)+1;            
            P(i,index+1:level) = P(i,index) + [1:(level-index)];
            index = level;
            i = i+1;
        else
            index = index - 1;
        end
     
end    
            


