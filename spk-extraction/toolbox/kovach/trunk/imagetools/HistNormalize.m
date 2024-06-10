function Zout = HistNormalize(Xin,Yin, window,ywindow,bitres)

% Z = HistNormalize(X,Y, window)
% Equates histogram of image X with Y, returning Z
% by stimulating the flow of the difference between
% the histograms towards a level value.
%
% window is an optional mask

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/imagetools/HistNormalize.m $
% $Revision: 676 $
% $Date: 2016-01-15 14:07:43 -0600 (Fri, 15 Jan 2016) $
% $Author: ckovach $
% ------------------------------------------------

%C. Kovach 2007

% g = normpdf(-20:1:20,0,2);
% G = g'*g;
% G = G./sum(G(:));

if nargin < 5 || isempty(bitres)
    bitres = 255;
    if min(size(Yin)) == 1
        bitres = length(Yin)-1;
    end
end

norm = 0;
if max(Xin(:)) >0 && max(Xin(:)) <=1
    Xin = Xin*bitres;
    norm = 1;
end
if max(Yin(:)) >0 && max(Yin(:)) <=1
    Yin = Yin*bitres;
    norm = 1;
end

for k = 1:size(Yin,3)
    X = Xin(:,:,k);
    Y = Yin(:,:,k);

    % LocalAv = conv2(X,G,'same');
    LocalAv = X;

    if nargin < 3
        xwindow = ones(size(X));
        ywindow = ones(size(Y));
    elseif nargin < 4
        xwindow = window;
        ywindow = window;
    else
        xwindow = window;
    end

    X = round(X);
    Y = round(Y);
    X(X<0) = 0;
    X(X>bitres) = bitres;
    %[h2,X2] = histc(X(:),0:bitres);
    h2 = hist(X(xwindow>0),0:bitres); 
    

    if min(size(Y)) > 1
        X1 = Y;
        %h1 = histc(X1(:),0:bitres)';
        h1 = hist(X1(ywindow>0),0:bitres)';
    else
        h1 = Y;
    end
    h2 = h2(:);
    h1 = h1(:);
    
    h1 = round(h1./sum(h1)*sum(h2));
    
    dh = h1-h2;

    DM = diag(ones(length(dh)+1,1),0) - diag(ones(length(dh),1),1);
    %DM(:,end) = [];

    DM(end,1) = 1;
    % dh(end+1)= 0;
    %psi = (DM'*DM)^-1*(DM'*dh)
    iDM = DM^-1;
    % flow = iDM*dh;

    Z = X;
    %toflow = flow;

    % X2 = Z;
    %[h2,X2] = histc(Z(:),0:bitres);
    h2 = hist(Z(xwindow>0),0:bitres);
    h2 = h2(:);
    dh = h1-h2;
    dh(end+1) = 0;
    toflow = iDM*dh;
    damp = 1;


    fig = gcf;
    hold off
    go = true;
    oldtot = 0;
    while go %sum(abs(toflow))>abs(sum(h1)-sum(h2)) 
    for lev = 0:bitres


        if toflow(lev+1) < 0 
            %match = find(X2 == lev);
            match = find(Z == lev);
            %rp = randperm(length(match));    
            [st,rp] = sort(LocalAv(match));    
            get = match(rp(1:min(round(damp*[length(match),abs(toflow(lev+1))]))));
            Z(get) = Z(get) - 1;
            toflow(lev+1) = toflow(lev+1) + length(get);
        end

        if toflow(lev + 2) > 0
            %match = find(X2 == lev);
            match = find(Z == lev);
            %rp = randperm(length(match));    
            [st,rp] = sort(LocalAv(match),'descend');    
            %rp = rp(end:-1:1);
            get = match(rp(1:min(round(damp*[length(match),abs(toflow(lev+2))]))));
            Z(get) = Z(get) + 1;
            toflow(lev+2) = toflow(lev+2) - length(get);
        end

    end

    h2 = hist(Z(xwindow>0),0:bitres)';
    dh = h1-h2;
    dh(end+1) = 0;
    toflow = iDM*dh;
    sum(abs(toflow))
    if ishandle(fig)
        plot(toflow)
        figure(fig)
        drawnow
    end
 dtot = sum(abs(toflow))-oldtot;
    if dtot == 0
        go = false;
    end
    oldtot = sum(abs(toflow));
    end
   
    Zout(:,:,k) = Z; %#ok<*AGROW>
end
if norm, Zout = Zout/bitres;end
