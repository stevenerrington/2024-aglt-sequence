function [ac,x] = auto_corr(ts, binsize, n_lags)
% INPUT: ts = timestamps (sorted)
%        binsize = binsize for binning timestamps.
%        n_lags = n_lags in the xcorr.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowen 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_it = false;
scaleopt = 'none';

ac = nan(n_lags,1);
x = -100:binsize:100 ; % bin centers.

if isempty(ts)
    return
end

edges = ts(1):binsize:ts(end);
c = histcounts(ts,edges);
[ac] = xcorr(c,n_lags,scaleopt);
ac = ac((n_lags+2):end);