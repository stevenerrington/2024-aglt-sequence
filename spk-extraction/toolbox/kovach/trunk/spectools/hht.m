

function [H,x, medf] = hht(x, varargin)

% H = hht(x)
%
% Performs the Hilbert-Huang empirical mode decomposition for the input
% signal. 
%
% [H,Rem] = hht(x,'stopfreq',stopfreq) 
%
% Compute the IMFs down to and including frequency stopfreq. The algorithm
% halts after the first component with median frequency < stopfreq. Rem =
% remainder of the original signal after subtracting sum(real(H),2).
%
%
% [EMD,Rem] = hht(x,'apply_hilbert',false) 
%
% Applies EMD without Hilbert tranform. Real valued IMFs are returned.
% Note that EMD = real(H)
%
% HHT works by decomposing a signal into a finite set of intrinsic mode
% functions (IMFs) and applying the Hilbert transform to each to obtain an
% analytic signal.
%
% IMFs have the following properties:
%
% 1. Number of peaks and number of zero crossings differ by no more than 1.
% 2. The upper and lower envelope of the signal averages to 0.
%
% These properties imply that the phase of the analytic signal is
% monotonically non-decreasing and well behaved.
%
% See N. E. Huang et al., Proceedings: Mathematical, Physical and
%   Engineering Sciences 454, 903 (1998).
%
% This version applies spline interpolation on a circular domain. The
% rationale for doing this is to avoid ambiguity and numerical instability
% at the edges of the domain. This results in edge effects which can
% be addressed through windowing or zero padding, just as when applying
% FFT-based filters.
%
% When adjacent values are constant, peaks are extracted at left-most corners
% except at saddle points for which they are assigned to both corners.
%
% i.e.
%              
%               (+)_____
%                 /     \
%        (+) ____/       \
%     \     /   (-) 
%      \___/ 
%    (-) 



% C. Kovach 2011

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/hht.m $
% $Revision: 763 $
% $Date: 2016-10-25 21:20:19 -0500 (Tue, 25 Oct 2016) $
% $Author: ckovach $
% ------------------------------------------------


stopfreq = 0;
maxiter = Inf;
apply_hilbert = true;
maxiter = 1e3;
norm = false;
fs = 1;

%%% For backward compatibility with earlier version of hht which 
if nargin > 1 && ~ischar(varargin{1})
    stopfreq = varargin{1};
    varargin(1) = [];
    if nargin >2 && ~ischar(varargin{1})
        maxiter = varargin{1};
        varargin(1) = [];
    end
end

i = 1;
while i <= length(varargin)
       
    switch lower(varargin{i})
        case 'stopfreq',   % lowest frequency to extract
            stopfreq= varargin{i+1};
            i = i+1;            
        case 'maxiter'   % Maximum number of iterations for each IMF
            maxiter = varargin{i+1};
            i = i+1;            
        case {'samplingfreq','fs'}   %Sampling frequency
            fs = varargin{i+1};
            i = i+1;            
        case {'apply_hilbert','hilbert','ht'}   %AIf false, Hilbert transform is not applied and only the real-valued IMFs are returned
            apply_hilbert = varargin{i+1};
            i = i+1;            
        case {'normalize','norm'}  %Precondition the signal by removing the mean and setting sd to 1.    
            norm = varargin{i+1};
            i = i+1;            
        otherwise
            error('%s is not a recognized keyword',varargin{i});
    end
    i = i+1;
end



if norm
    sdx = std(x);
    m = mean(x);
    x = (x-m)./sdx;
end

h = x(:);

n = length(x);
xpt = (1:n)';

tol = 1e-6;

imfn = 1;
ds = [];
stop = 0;
imfn = 1;

while sum(x.^2)~=0 && ~stop
    
    
   
    d = 1;
    
    h = x(:);
    iter = 1;
    fprintf('\nIMF %3i, Iter ',imfn)
    fprintf('%4i',iter)
    mm = 0;

    while d>tol
        
%         ds(iter) = d;
        
        [pk,pksign,zeroc] = getpeak(h);

%         nzocs(iter) = sum(zeroc);
%         npks(iter) = sum(pk);
        
%         if isempty(pk)
%             stop = 1;
%             break
%         end
       if sum(pk)<2 || iter > maxiter
           if iter > maxiter
               warning('IMF Failed to converge after %i iterations',iter-1)
           end
          break
       end
        pkipos = find(pk&pksign>0);
        pkineg = find(pk&pksign<0);
        pkhpos = h(pk&pksign>0);
        pkhneg = h(pk&pksign<0);
        
        %%% These functions wrap the first and last peaks to the end and
        %%% beginning, respectively. Splines are computed on a circular
        %%% domain in order to avoid numerical problems associated with
        %%% extrapolation at the edges.     
        npad = 10;
        intpi = @(x) [circnget(x,-npad,n);x;circnget(x,npad,n)];
        intph = @(x) [circnget(x,-npad);x;circnget(x,npad)];
%         
%         intpi = @(x) [x(end-4:end)-n;x;x(1:5) + n];
%         intph = @(x) [x(end-4:end);x;x(1:5)];
%         
        posenv = interp1( intpi(pkipos) ,intph(pkhpos),xpt,'spline');  %  splineinterpolation of upper envelope
        negenv = interp1( intpi(pkineg) ,intph(pkhneg),xpt,'spline');  %#ok<*FNDSB> % spline interpolation of lower envelope
        
        m = (posenv+negenv)/2; %The envelope average
        
%         mm = mm+m;
        
        hnew = h-m; % New candidate signal after subtracting the envelope average    
        
         d = sum((h-hnew).^2)./sum(h.^2); %change
%          ds(end+1) = d;
         h = hnew;
         
         iter = iter+1;
         fprintf('\b\b\b\b%4i',iter)

    end
    
    if apply_hilbert  %%% apply the Hilbert transform to get an analytic signal
        hh = hilbert(h);
    else
        hh = h;
    end
    
    ph =hh./abs(hh); % unit phasor
    iff = diff(unwrap(atan2(imag(ph),real(ph))))./(2*pi)*fs; % insantaneous frequency computed from the derivative of the phase-ange.
    
    medf(imfn) = median(iff);
    
    fprintf('\tMedian freq. %3.2g',medf(imfn))
        
    H(:,imfn) = hh;
    
    x = x-h;
    
    imfn = imfn+1;
    if median(iff) < stopfreq
        break
    end

end


%%%%%%

function xout = circnget(x,n, mlt)

% Return the first n values of a vector, x, if n>0 and the last abs(n) if n<0
% If there are fewer than n values in x, then cycle through x as many times as needed.


if nargin < 3
    mlt = 0;
end

xout = [];

nx = length(x);

if n > 0
    geti = mod(0:n-1,nx)+1;
    xout = x(geti,:) +  mlt*(fix((0:n-1)/nx)+1)';
elseif n < 0
    
    geti = mod(nx+n+1:nx,-nx);
    xout = x(geti+end,:) + mlt*( fix( ( nx+n+1:nx )/nx-1)-1)' ;
end    


   

%%%%%%%%%%
function [pk,pksign,zeroc] = getpeak(x)

% Extract all positive and negative peaks and zero crossings 

tol = eps;

dx = diff([x(end);x;x(1)]);  %Again, compute difference on circular domain

dx = dx.*(abs(dx)>tol); %Apply tolerance threshold

sdx = sign(dx) ;

dsdx = [sdx(end);sdx] - [sdx;sdx(1)];

zeroc = sign([x;x(1)])~=sign([x(end);x]);
% dsdx = diff(sdx);

pk = abs(dsdx)>0;
% pk([1 end]) = true;
pksign = sign(dsdx.*pk);

%%% When one or more adjacent values are equal, assign a peak  to
%%% the first inflection for concave or convex regions, and to both
%%% inflections for saddle points.
if any(dx==0) && any(dx~=0) 
    fpks = pksign(pk);
    pk(pk) = [fpks(end);fpks(1:end-1)] ~= fpks;
    pksign = pk.*pksign;
end
    
    

pk = pk(2:end-1);
pksign = pksign(2:end-1);
zeroc = zeroc(2:end-1);

    



