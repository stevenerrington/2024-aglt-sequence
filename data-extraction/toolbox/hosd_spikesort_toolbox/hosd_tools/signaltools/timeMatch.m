
function varargout = timeMatch(samp, ref, P, Q, plotx,useConv2)

%d = timeMatch(samp, ref, P, Q)
% Matches waveform in samp to ref after resampling by P,Q using
% convolution

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

rsampn = 0; %Order of resampling; 0 is nearest neighbor

if nargin < 6
    useConv2 = false;
end

if nargin < 3
    
    P = 2040.4*10;
    Q = 120*10;
    
end

if nargin < 5

%     plotx = 0;
    plotx = 1;
end

samp = resample(samp(:),P,Q,rsampn);


%rzdiv = ifft(fft(rz.^2).*fft([ones(length(sz),1);zeros(ln - length(sz))]) );
if useConv2
    sz = zscore(samp(:));
    %rz = flipud(zscore(ref(:)));
    rz = zscore(ref(:));
    rx = xcorr2(sz,rz);
    %X = conv2(sz,rz)./sqrt(sum(sz.^2).*sum(rz.^2));
else
    
    ln = max([length(samp),length(ref)]);
    sz = zeros(ln,1);
    rz = zeros(ln,1);
    square = zeros(ln,1);
    sz(1:length(samp)) = zscore(samp(:));
    %rz(1:length(ref)) = flipud(zscore(ref(:)));
    rz(1:length(ref)) = zscore(ref(:));
    square(1:length(samp)) = 1;
    rz(rz > 5) = 5;
    sz(sz>5) = 5;
    
    rzdiv = ifft(fft(square).*conj(fft(rz.^2)));
    
    X = ifft(fft(sz).*conj(fft(rz)))./sqrt(sum(sz.^2).*rzdiv);
end


% rzdiv = conv2(rz.^2,ones(size(sz))); 
% X = conv2(rz,flipud(sz))./sqrt(sum(sz.^2).*rzdiv);

[mx,mxi] = max(X);
%d = mxi - length(sz);
%d = mxi - length(sz)./2;
if mxi > length(X)./2  %If the signals are displaced by more than N/2 points, this won't work! 
    d =  length(X) - mxi +1;
else
    d = -mxi +1;
end

if mx < .5
    suspect = true;
else
    suspect = false;
end


if plotx
    %figure,
    subplot(2,1,1);
    plot(flipud(X(:)));
    
    %rz = flipud(rz);
    subplot(2,1,2);
%     rzind = 1:length(samp);
%     szind = [1:length(samp)] + d;
%     szind(rzind < 1 ) = [];
%     rzind(rzind < 1 ) = [];
    %plot([1:length(sz)]+d,sz,'r',[1:length(rz)],rz(end:-1:1),'b')
    plot([1:length(sz)]+d,sz,'r',[1:length(rz)],rz,'b')
    drawnow
end
varargout = {d,suspect};

if nargout >= 3
    varargout{3} = samp;
end

if nargout == 4
    varargout{4} = mx;
end
