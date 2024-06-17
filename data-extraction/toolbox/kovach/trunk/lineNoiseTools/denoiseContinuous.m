function [denoised,blocksize] = denoiseContinuous(x,fs, freqrange,windowrange)

% [denoised,blocksize] = denoiseContinuous(x,fs, freqrange,windowrange)
%
% Removes line noise  from a conintuous signal by dividing the signal into a set
% of overlapping windows and searching for peaks near multiples of a given fundamental frequency
% (default 60Hz) within each windowed segment. 
% 
%   x - signal
%  fs - sampling frequency
%  freqrange - which range of fundamental frequencies to search for the
%              periodic signal: [flo fhi] (default is [59.5 60.5])
%  windowrange - range of acceptable window sizes (default 4 - 8 sec)
%  Windows are set to be as near to integer multiples of (Flo + Fhi)/2.
% 
% This function calls denoise60Hz

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/lineNoiseTools/denoiseContinuous.m $
% $Revision: 177 $
% $Date: 2013-03-13 18:00:59 -0500 (Wed, 13 Mar 2013) $
% $Author: ckovach $
% ------------------------------------------------

%C. Kovach 2008

use_parallel = false; %Use parallel computing toolbox, if available
                      % Overhead largely negates the advantage of parallel
x = double(x(:));

if nargin < 3
    freqrange = [59.5 60.5];
end

if nargin < 4
%     windowsize = 360;
    windowrange = [4 8];
end

if length(windowrange) == 2
    wrange = round(windowrange(1)*fs)./fs:1/fs:windowrange(2);
else
    wrange = windowrange;
end
% wrm = rem(wrange,1/mean(freqrange)); % Find nearest point to a complete period within the acceptable range

dvfpart = wrange*mean(freqrange) - fix(wrange*mean(freqrange));
phdiff = abs(1-exp(1i*dvfpart*2*pi));

%%% size of block in sampling points. It should be as close to a multiple of line noise fundamental as possible.
% blocksize = round(windowsize*fs./mean(freqrange)); %size of block in sampling points. It should be as close to a multiple of line noise fundamental as possible.
[mn,mni] = min(phdiff);
blocksize = round(wrange(mni)*fs);            

% edge = 100; % number of points on the edges to ignore due to edge effects
edge = 0; % number of points on the edges to ignore due to edge effects
avgWin = ceil((blocksize - 2*edge - 1)./2);  %Number of points to average with succeding trial to avoid discontinuities.

N = length(x);

stepsize = blocksize - avgWin - 2*edge; 

denoised = zeros(size(x));

dnblock = denoise60Hz(x(1:blocksize,:),fs);

denoised(1:blocksize - edge,:) = dnblock(1:blocksize - edge,:);

steps = [stepsize:stepsize:N-blocksize];

%%% The signal is denoised in overlapping blocks where the center of one block
%%% aligns with the edge of the next.
%%%
%%% The overlapping signals obtained this way are weighted by w = .5+.5*cos(ax+b) where
%%% a and b are chosen such that w = 1 in the center of one block and w = 0
%%% at the edges. The final denoised signal is the weighted average of the
%%% the two, i.e. xdn = w*xdn1 + (1-w)*xdn2. This is done to suppress any
%%% edge effects which might arise from any mismatch between window
%%% size and the periodicity of the noise.



avgFun = .5-.5*cos([0:avgWin-1]'./avgWin.*pi); %Weighting function over half a block

avgFun = avgFun(:,ones(size(x,2),1));

% DN = [];
if ~use_parallel
    for s = steps

        dnblock = denoise60Hz(x(s:s+blocksize,:),fs,freqrange); % Remove fundamental and harmonics from the block. See DENOISE60Hz

        %averaging together the overlapping blciks
        denoised([s:s+avgWin-1]+edge,:) = (1-avgFun).*denoised([s:s+avgWin-1]+edge,:) + avgFun.*dnblock([1:avgWin]+edge,:);

        denoised([s:s + stepsize-1] + edge + avgWin,:) = dnblock([1:stepsize]+edge + avgWin,:);

    %     DN(:,end+1) = dnblock;
    end
    
    %%% Now do the remaining segment
    dnblock = denoise60Hz(x(s+stepsize:end,:),fs); % The last block is as big as the remaining segment of data
    denoised([s:s+avgWin-1]+edge + stepsize,:) = (1-avgFun).*denoised([s:s+avgWin-1]+edge + stepsize,:) + avgFun.*dnblock([1:avgWin]+edge,:);
    denoised(s+edge+avgWin + stepsize:end,:) = dnblock(1+edge + avgWin:end,:);

else
    
%     denoised1 = zeros(size(x));
%     denoised2 = zeros(size(x));
    ST = chopper([0 blocksize],steps,1);
    xx = x(ST);
    denoised1 = zeros(avgWin+edge,length(steps));
    denoised2 = zeros(avgWin+edge,length(steps));
    
    parfor si = 1:length(steps)
        
%         s = steps(si);
        dnblock = denoise60Hz(xx(:,si),fs,freqrange); % Remove fundamental and harmonics from the block. See DENOISE60Hz

        %averaging together the overlapping blciks
        denoised1(:,si) =  avgFun.*dnblock([1:avgWin]+edge,:);

        denoised2(:,si) = (1-avgFun).*dnblock([1:stepsize]+edge + avgWin,:);

    %     DN(:,end+1) = dnblock;
    end
    denoised(ST(1)+edge:ST(stepsize+edge + avgWin,end)) =...
    [denoised1,zeros(size(denoised1,1),1)]+...
    [denoised(ST([1:avgWin]+edge,1)).*(1-avgFun),denoised2(:,1:end)];
    
    s = steps(end);
    
    % Now do the remaining segment
    dnblock = denoise60Hz(x(s+stepsize:end,:),fs); % The last block is as big as the remaining segment of data
    denoised([s:s+avgWin-1]+edge + stepsize,:) = denoised([s:s+avgWin-1]+edge + stepsize,:) + avgFun.*dnblock([1:avgWin]+edge,:);
    denoised(s+edge+avgWin + stepsize:end,:) = dnblock(1+edge + avgWin:end,:);

end
    
