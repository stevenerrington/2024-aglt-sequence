function Y = smfilt(X,window)

%Smooting filter which applies the given window

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if size(X,1) == 1
    X = X';
end

nw = length(window);
nx = length(X);

B = zeros(size(X).*[2 1]);
B(1:ceil(nw./2)) = window(floor(nw/2)+1:end);
B(end-floor(nw./2)+1:end) = window(1:floor(nw/2));


F = fft(B);

Y = real( ifft( fft(cat(1,X,zeros(size(X)))).*conj(F(:,ones(size(X,2),1)) )) );

Y = Y(1:nx,:);

