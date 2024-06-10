function xdn = artifactSlayer(x)

% function xdn = artifactSlayer(x)
%
%    x: struct with fields
%         .dat :  data from one channel as a column vector
%         .fs :  sampling frequency
%
%    xdn: same struct as x with denoised data in xdn.dat.

% C. Kovach 2015

%%
%%% Loop through these bandwidths.
filter_bandwidths = [ 1 10 .25 5 ];

do_plot = true;

if do_plot
    fig = figure;
end

xdn = x;

for k = 1:length(filter_bandwidths)
    if do_plot
        figure(fig)
        subplot(length(filter_bandwidths)+1,2,2*k-1)
    end
    
    %%% Apply dbt denoising at successively wider bandwidths
    [xdn.dat,~,dbx] = dbtDenoise(xdn.dat,xdn.fs(1),filter_bandwidths(k),'makeplots',do_plot,'remove spikes',false,'filter above',20);
    
    if do_plot       
        figure(fig)
        subplot(length(filter_bandwidths)+1,2,2*k)
       dbx.specgram;
       title(sprintf('Pass %i, BW %0.2f Hz',k,filter_bandwidths(k)))
       drawnow
    end
end
if do_plot
    figure(fig)
    subplot(length(filter_bandwidths)+1,1,length(filter_bandwidths)+1)
   t = (0:length(x.dat)-1)/x.fs(1);
    plot(t,x.dat,'b');
    hold on
    plot(t,xdn.dat,'r')
    legend({'Before filtering','After filtering'})
   drawnow
end
    