function bem = BEM(vol,dp)

% computes boundary element forward model
%
% tess - tesselation for each bondary
% X - vertices for each boundary
% s - conductance for each region
% jpX - primary current dipole locations
% jp - primary moments
% M - measured locations

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/BEM.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------



dim = 3;

scale = 1;

if isfield(vol,'bel')
    tess = {vol.bel.tri};
    X = {vol.bel.pnt};
    s = vol.cond;

    ds = diff(s);
    ms = s(1:end-1)+ds/2;
%     nvert = sum(cellfun('length',X));
    allds = [];
    allms = [];
    layer = [];
    for i = 1:length(X)
        if dim == 2
            TR{i} = TriRep(tess{i},X{i}(:,1),X{i}(:,2));
        else
            TR{i} = TriRep(tess{i},X{i}(:,1),X{i}(:,2),X{i}(:,3));
        end
        FN{i} = faceNormals(TR{i});
        FC{i} = incenters(TR{i});
        allds =cat(1,allds,ones(size(FN{i},1),1)*ds(i));    
        allms =cat(1,allms,ones(size(FN{i},1),1)*ms(i));
        layer =cat(1,layer,ones(size(FN{i},1),1)*i);

    end

    allFC = cat(1,FC{:});
    allFN = cat(1,FN{:});
    pnt = [X{:}];
else
     bem = vol;
     allFC = bem.FC;
     allFN = bem.FN;
     allds = bem.ds;
     allms = bem.ms;
     layer = bem.layer;
    pnt = [bem.pnt{:}];
end
bem.scale = scale;
    
nface = size(allFN,1);
nvert = size(pnt,1);

%Edge vectors
E{1} = pnt(bem.tri(:,2),:) - pnt(bem.tri(:,1),:);
E{2} = pnt(bem.tri(:,3),:) - pnt(bem.tri(:,2),:);
E{3} = pnt(bem.tri(:,1),:) - pnt(bem.tri(:,3),:);

normvec = cross(E{1},-E{3});
bem.area = sqrt(sum(normvec.^2,2))./scale;
% normvec = normvec./bem.area(:,[1 1 1]);
% normvecdir = sign(sum(normvec.*bem.FN,2));

% normals to the edges
% for i= 1:3
%     M{i} = cross(E{i},normvec);
%     edgelength(:,i) = sqrt(sum(M{i}.^2,2));
%     M{i} = M{i}./edgelength(:,i*[1 1 1]);
%     
% end



% R= zeros(nface);
% 
% for i = 1:nvert    
%     dd = ( repmat(allFC(i,:),nface,1) - allFC )./scale;    
%     R(i,:) = 1./(4*pi)*bem.area.*allds.*sum(dd.*allFN,2)./sqrt(sum(dd.^2,2)).^dim;
% %     i
% end

%Matrix which accomplishes  RSmat*A = A(tri'(:))
U = speye(nvert,nvert);
tri = bem.tri';
RS1mat = U(:,tri(:))';

%Matrix which accomplishes RS2mat*B = kron(B,ones(3,1);
RS2mat = kron(speye(nface),ones(3,1));


%Matrix to compute teh inner product with normal for each face
% FN = bem.FN'*RS2mat';
% SFN = kron(speye(3*nface),ones(1,3)) * diag(sparse(FN(:)));
FN = (RS2mat*bem.FN)';

%reshapes dd'(:) so that the diagonal is ordered according tri(:)
RS3mat = kron(RS1mat,speye(3));


%sums over adjacent 3 diagonals
RS4mat = kron(speye(3*nface),ones(1,3));

% Matrix operations from from right to left:
%  1. orders d(:) according to tri(:)
%  2. multiplies element wise with FN(:) where FN is face normals
%  3. Takes the sum over three adjacent diagonals, which gives the inner
%      product of the function value  with the normal for each vertex of the given facet.
%  4. Multiply by facet area and difference conductivity aross compartment
%  5. Returns from facet space to vertex space by summing across facets for
%  each vertex 
FNmat = 1./6*RS1mat'*diag(sparse(kron(bem.area.*allds,ones(3,1))))*RS4mat*diag(sparse(FN(:)))*RS3mat./scale.^2;

bem.FNmat = FNmat;

R = zeros(nvert);
for i = 1:nvert 

    dd = (repmat(pnt(i,:)',1,nvert) - pnt')./scale; 
    FD = 1./(4*pi)*dd./repmat(sqrt(sum(dd.^2)).^3,3,1);
    
    R(i,:) = FNmat*FD(:);
end

R(isnan(R(:))) = 0; 

for i = 1:length(R), 
    R(i,i) = allms(i) +R(i,i);
end



[L,U,P] = lu(R,'vector');

bem.Rdata.RLU = L+U;  %combining into one matrix for efficiency
bem.Rdata.dgl = diag(L);
bem.Rdata.dgu = diag(U);
bem.Rdata.permvec = P;

bem.allms = allms;
bem.allds = allds;
bem.allFC = allFC;
bem.allFN  = allFN;
bem.layer = layer;
% bem.RL = L;
% bem.RU = U;
bem.fitforward = @(X,bem) forward(X,bem);
bem.computeV = @(X,bem) computeV(X,bem);

if nargin > 1
 
  bem = forward(dp,bem);
end    

%%% 

function bem = forward(dp,bem)
    
    
DPmat = [];
pnt = [bem.pnt{:}];

for i = 1:length(dp);


    d = (pnt - repmat(dp(i).x,size(pnt,1),1) )./bem.scale;
%     VP = VP +(1./4*pi).*(d*dp(j).p')./sqrt(sum(d.^2,2)).^dim; 

    DPmat = cat(2,DPmat,(1./4*pi)*d./repmat(sqrt(sum(d.^2,2)).^3,1,3));  

end

p = cat(2,dp.p);

L(bem.Rdata.permvec,:) = tril(bem.Rdata.RLU) - diag(bem.Rdata.dgu);
U = triu(bem.Rdata.RLU) - diag(bem.Rdata.dgl);

bem.forward.dp = dp;
bem.forward.LF = U\(L\DPmat); % lead field matrix for vertices
bem.forward.V = bem.forward.LF*p(:); %Vertex voltages



    
%%%%

function bem = computeV(X,bem)
%compute voltage at arbitrary locations

if ~isfield(bem,'forward')
    bem = fitforward(bem);
    
end
[in,on] = insidebem(X,bem);
q = -diff(cat(1,in,[ 0 0]));
inlayer = [1:size(in,1)]*q;
onlayer = [1:size(in,1)]*on;

pnt = [bem.pnt{:}];

nvert = size(pnt,1);
nrec = size(X,1);    
MR = zeros(nrec,nvert);
DPmat = [];

for i = 1: nrec
       
    if onlayer(i)
            cond(i) = unique(bem.ms(bem.layer == onlayer(i)));
    else

            incond = bem.ms + bem.ds./2;
            cond(i) = unique(incond(bem.layer == inlayer(i)));
    end

        ddp = (X(i,:) - repmat(bem.forward.dp(i).x,size(X,1),1) )./bem.scale;
%     VP = VP +(1./4*pi).*(d*dp(j).p')./sqrt(sum(d.^2,2)).^dim; 

    DPmat = cat(2,DPmat,(1./4*pi)*ddp./repmat(sqrt(sum(ddp.^2,2)).^3,1,3));  


    dd = (repmat(X(i,:)',1,nvert) - pnt')./bem.scale; 
    FD = dd./repmat(sqrt(sum(dd.^2)).^3,3,1);

    MR(i,:) = FNmat*FD(:);
end

p = cat(2,dp.p);

bem.forward.MLF = DPmat - MR*bem.forward.LF; % lead field matrix for measured locations
bem.forward.MV = bem.forward.MLF*p(:);

    