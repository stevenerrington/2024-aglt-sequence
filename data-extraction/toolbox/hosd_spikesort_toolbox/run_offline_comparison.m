outdir = fullfile('hosd_output','online');

for simulation_set=1:3
    for lev=1:4

    %%
        prfigs = [13 14 15 8888];
        fldn = {'hosd','wvlt'};
        tmpdir1 = tempname;
       
        vnames = {'\# Spikes','\# Detected','\# algn','\# Det-algn','TP','FP','FPnoise','FPothclus','det-TP','nrAssgn'};
        pdfnames = {};
        mdtype = [3 1];
        for deti = 1:length(fldn)
            [osort(simulation_set,lev).(fldn{deti}).perf,osort(simulation_set,lev).(fldn{deti}).nrAssigned,osort(simulation_set,lev).(fldn{deti}).assigned,osort(simulation_set,lev).(fldn{deti}).params] = runSimulatedEval(simulation_set, lev, mdtype(deti),false, true);   
     
%        [osort(simulation_set,lev).wvlt.perf,osort(simulation_set,lev).wvlt.nrAssigned,osort(simulation_set,lev).wvlt.assigned,osort(simulation_set,lev).wvlt.params] = runSimulatedEval(simulation_set, lev, 1,false, true);   
            mkdir(fullfile(tmpdir1,fldn{deti}));
            fg=figure(81218)
            cla
            text(.5,.5,sprintf('Performance for %s',fldn{deti}),'Fontsize',20,'horizontalalignment','center')
            axis off
            pdfnames{end+1}=fullfile(tmpdir1,fldn{deti},sprintf('fig%i.eps',0));
            print(fg,'-depsc', pdfnames{end});
            for k = 1:length(prfigs)
             pdfnames{end+1}=fullfile(tmpdir1,fldn{deti},sprintf('fig%i.eps',k));
             print(prfigs(k),'-depsc', pdfnames{end});
             
            end
            
        end
        fg = figure(12321);
        
        nspikes = size(osort(simulation_set,lev).hosd.perf,1);
        x1=[1:nspikes;osort(simulation_set,lev).hosd.perf'];
        x2=[1:nspikes;osort(simulation_set,lev).wvlt.perf'];
           
%           txt1 = sprintf('Sim. set:%i\nNoise lev: %i\n\n\\begin{tabular}{l%s}\\hline\\multicolumn{%i}{c}{\\textbf{HOSD Detection}} \\\\ \\hline &%s \\\\ \\hline %s \\hline \\multicolumn{%i}{c}{\\textbf{Wavelet Detection}} \\\\ \\hline %s \\end{tabular}',...
%               simulation_set,lev,repmat('|l',1,length(vnames)+1),length(vnames)+1,sprintf(' %s & ',vnames{:}),...
%               sprintf([' spike %i & ',repmat(' %i &',1,length(vnames)),' \\\\ '],x1(:)),length(vnames)+1,sprintf([' spike %i & ',repmat(' %i & ',1,length(vnames)),' \\\\ '],x2(:)));
%     
          txt1 = sprintf('Sim. set:%i\nNoise lev: %i\n\n\\begin{tabular}{l%s}\\hline\\multicolumn{%i}{c}{\\textbf{HOSD Detection}} \\\\ \\hline &%s \\\\ \\hline %s \\end{tabular}',...
              simulation_set,lev,repmat('|l',1,length(vnames)+1),length(vnames)+1,sprintf(' %s & ',vnames{:}),...
              sprintf([' spike %i & ',repmat(' %i &',1,length(vnames)),' \\\\ '],x1(:)));
          txt2 = sprintf('\\begin{tabular}{l%s}\\hline\\multicolumn{%i}{c}{\\textbf{Wavelet Detection}} \\\\ \\hline &%s \\\\ \\hline %s \\end{tabular}',...
              repmat('|l',1,length(vnames)+1),length(vnames)+1,sprintf(' %s & ',vnames{:}),...
              sprintf([' spike %i & ',repmat(' %i &',1,length(vnames)),' \\\\ '],x2(:)));
       

          th2(1) = text(-.1,.66,txt1,'interpreter','latex')  ; 
          th2(2) = text(-.1,.33,txt2,'interpreter','latex')  ; 
          axis off
          pdfnames{end+1}=fullfile(tmpdir1,'perf.eps');
          print(fg,'-depsc', pdfnames{end});
          outfn = fullfile(outdir,sprintf('summary_set%i_lev%i.pdf',simulation_set,lev));
         com1 = sprintf('/usr/local/bin/gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s %s',...
                           outfn,sprintf(' %s',pdfnames{[end,1:end-1]}));

           com2=sprintf('/Library/TeX/texbin/pdfcrop %s %s',outfn,outfn);
%  [err2,out2]=system(com2);

        fid = fopen('temp.sh','w');
        fprintf(fid,'LD_LIBRARY_PATH=\n%s\n\n%s',com1,com2);
        fclose(fid);
        system('chmod 755 temp.sh');
%         [err,out] = system('/opt/X11/bin/xterm -e "./temp.sh;read;"');
        [err,out] = system('./temp.sh');
         close all
         
    end
end