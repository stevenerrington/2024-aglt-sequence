function [xfilt,spike] = spikefilter(x,fs,spike)


%Simple script to remove outliers through iterative thresholding

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/DBT/spikefilter.m $
% $Revision: 491 $
% $Date: 2014-05-06 17:57:17 -0500 (Tue, 06 May 2014) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    spike.threshold = 10; % This is the threshold used in detecting spikes. 
                          % Z-score is computed iteratively until
                          % no points exceed this value. The threshold is set
                          % high by default because the main purpose here is to avoid
                          % distortion of kurtosis used in the
                          % kurtosis-threshold.
    spike.smoothwindow = .2;% Apply hanning window of given duration to smooth the spike filter.

end

spks = false(size(x));
newspks = true;
while any(newspks)
       z = (x-mean(x(~spks)))/std(x(~spks));
       newspks = abs(z)>spike.threshold &~spks; 
       spks = newspks | spks;
end
win = hanning(ceil(spike.smoothwindow.*fs));
spike.filter = exp(convn(log(1-spks+eps),win,'same'));
spike.filter(spike.filter<0)=0;
xfilt = x.*spike.filter;