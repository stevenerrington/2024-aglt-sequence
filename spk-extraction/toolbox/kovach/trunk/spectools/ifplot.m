function varargout  = ifplot(H,x,decim,varargin)

%  plh = ifplot(H,x,decim)
% Plot instantaneous frequency (IF) and amplitude for a matrix containing
% complex analytic signals. IF is plotted as a line whose color represents
% analytic amplitude.
%
% 'decim' is the factor by which to decimate rows before plotting.
% Instantaneous frequency and analytic amplitude are computed before
% decimation.
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/ifplot.m $
% $Revision: 46 $
% $Date: 2011-06-09 12:39:59 -0500 (Thu, 09 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

%
% C. Kovach 2011


if nargin < 2 || isempty(x)
    x = (1:size(H,1))';
else
    x = x(:);
end

if nargin < 3 || isempty(decim)
    decim = 1;
end


PH = H./abs(H); %Unit phasor

ang = atan2(imag(PH),real(PH)); % Phase angle

IF = diff(unwrap(ang))/2/pi/diff(x(1:2));
IF(end+1,:) = 0;


% linewidths = std(IF)*.2;

% %Analytic amplitude
% A = abs(H(1:decim:end,:));
% %Instantaneous frequency
% IF = IF(1:decim:end,:);

%Analytic amplitude
A = resample(abs(H),1,decim);
%Instantaneous frequency
IF = resample(IF,1,decim);

x = x(1:decim:end);

X = cat(1,x,flipud(x));

nf = size(H,2);
% nx = length(x);


% X = repmat(xx,1,nf);

% Y = cat(1,IF,flipud(IF)) + .5*cat(1,repmat(linewidths,nx,1),-repmat(linewidths,nx,1));
Y = cat(1,IF,flipud(IF));

Z = cat(1,A,flipud(A));

for i = 1:nf
    varargout{1}(i) = patch(X,Y(:,i),Z(:,i),...
                    'edgecolor','interp','facecolor','none','linewidth',2,varargin{:});
%    shading interp
end

% colormap hot
% set(gca,'color','k')

if nargout > 1
    varargout{2} = x(1:decim:end);
end





    
    
    
    
    

