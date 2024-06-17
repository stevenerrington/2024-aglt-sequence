

function [R,srtmeas,srti ] = rowsmooth(X,n,sortmeas,downsample,circular)


%Sorts the colums of X according to the magnitude of sortmeas, smooths across rows with a gaussian kernel of width n and downsamples the columns. 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 5
    circular = false;
end

if nargin < 4    
    downsample = 1;
end

if nargin < 3 || isempty(sortmeas)    
    sortmeas = 1:size(X,2);
end

[srtmeas,srti] = sort(sortmeas);

if length(n) ==1 
    g = gausswin(n)';
    g = g./sum(g);
else
    g = n;
end

if ~circular
    R = convn(X(:,srti),g,'same');
else
    
    newg = zeros(1,size(X,2));
    newg([1:ceil(length(g)/2),end-floor(length(g)/2)+1:end]) = fliplr(fftshift(g(:)'));
    fg = fft(newg);
    
    R = (ifft(fft(X(:,srti),[],2).*repmat(fg,size(X,1),1),[],2));
end

if downsample > 1
    R = R(:,1:downsample:end);
    srtmeas = srtmeas(1:downsample:end);
end

