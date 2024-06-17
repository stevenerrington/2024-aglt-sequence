
function [Dst,Dsgn,pkcorr] = cumdist(X,Y,order)

if isscalar(Y)
    order = Y;
    Y = X;
end

X = cat(1,X,zeros(size(X)));
Y = cat(1,Y,zeros(size(Y)));

FX = fft(X);
FY = fft(Y);

XX = ifft(FX.*permute(conj(FY),[1 3 2]));

if mod(order,2)==0
    mx = max(abs(XX));
    Dsgn = squeeze(sign(sum(XX.*(abs(XX)==mx))));
    if size(X,2)==1
        Dsgn = Dsgn';
    end
else
    Dsgn =[];
end
% sKX= mean(ifft(abs(fft(X)).^2).^order)-mean(X.^order).^2;
% sKY= mean(ifft(abs(fft(Y)).^2).^order)-mean(Y.^order).^2;
% Dst = (squeeze(mean(XX.^order))-mean(X.^order)'*mean(Y.^order))./sqrt(sKX'*sKY);
% sKX= sum(ifft(abs(FX).^2).^order) - mX.^2;
% sKY= sum(ifft(abs(FX).^2).^order)- mY.^2;
% Dst = (squeeze(sum(XX.^order))-mX'*mY)./sqrt(sKX'*sKY);
sKX= cumulant(ifft(abs(FX).^2),order);
sKY= cumulant(ifft(abs(FY).^2),order);

N = sqrt(abs(sKX'*sKY));
if size(X,2)==1 || size(Y,2) ==1
    N = N';
end
Dst = squeeze(cumulant(XX,order))./N;
if size(X,2)==1
    Dst = Dst';
end
%%% Make sure the matrix is symmetric
if isequal(X,Y)
    Dst = (Dst+Dst')/2;
end

if nargout > 2
    mxXX= max(ifft(abs(FX).^2));
    mxYY= max(ifft(abs(FY).^2));
    if mod(order,2) == 0
        mxXY = max(abs(XX));
    else
        mxXY = max(XX);
    end
    N = sqrt(mxXX'*mxYY);
    if size(X,2)==1 || size(Y,2) ==1
        N = N';
    end
    pkcorr = squeeze(mxXY)./N;
    
    if size(X,2)==1
        pkcorr = pkcorr';
    end
    
end