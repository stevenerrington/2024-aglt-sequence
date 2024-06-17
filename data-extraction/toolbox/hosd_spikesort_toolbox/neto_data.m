

params.tolerance = 1e-3;

params = HOSD_default_params(params);

Fs = 3e4;

ddir = 'neto_data/2015_09_03_Pair_9_0';
adcfile = dir(fullfile(ddir,'adc*bin'));
ampfile = dir(fullfile(ddir,'amp*bin'));


fid = fopen(fullfile(adcfile.folder,adcfile.name),'r');
datatype = 'uint16';

adc = fread(fid,[datatype,'=>double']);
adc = reshape(adc,8,length(adc)/8)';
fclose(fid);


% [~,adcchan] = max(kurtosis(adc));
[~,adcchan] = max(var(adc));

z= zscore(adc(:,adcchan));
[adc_detect,hosadc,zadc] = HOSD_spike_detection(struct('dat',z,'fs',Fs,'chan',3),struct('ncomp',3,'lowpass',6000,'windur',.01,'highpass',100));


adc_cluster = sort_spikes(adc_detect);

% zresid = zadc-sum(hosadc.xrec(zadc),2);


keepcl = (max(adc_cluster.avg_waves)-min(adc_cluster.avg_waves))./nanstd(zadc)>5;
target = adc_cluster.spike_sep(keepcl);

nchan = 128;
chan = 100;

switch datatype
    case 'single'
        nbytes = 4;
    case 'double'
        nbytes = 8;
    otherwise
        nbytes = round(log2(double(intmax(datatype))-double(intmin(datatype)))/8);
end

fid = fopen(fullfile(ampfile.folder,ampfile.name),'r');
clear amp;
for k = 1:length(chan)
    fseek(fid,nbytes*(chan(k)),-1);
    amp(:,k) = fread(fid,[datatype,'=>double'],nbytes*nchan-2,'n');
  if k ==1 && length(chan)>1
        amp(end,length(chan))=0;
    end
end
    fclose(fid);
  
ampfilt = iterz(zscore(hpfilt(amp,[Fs 300])),20);


[Ttru,tt] = chopper([-1 1]*.001,[target.spike_times],Fs,length(ampfilt));%-mean(ttav),Fs);
MA = [];
for k = 1:size(amp,2)
    MA(:,k) = nanmean(ampfilt(Ttru + (k-1)*length(ampfilt)),2);
end
msq = nanmean(ampfilt.^2);
mamsq = nanmean(MA.^2)- msq/size(Ttru,2);
mamsq(mamsq<0)=0;

%%% Find the channel with the greatest SNR
% SNRs = 10*log10(mamsq./msq);
SNRamp = 20*log10(((max(MA)-min(MA))/2)./sqrt(msq));
SNRsd = 10*log10(std(MA)./sqrt(msq));
[~,keepch] = max(SNRamp);
zamp = ampfilt(:,keepch);


zamp(isnan(zamp))=0;
[amp_detect,amp_hos,ampz] = HOSD_spike_detection(struct('dat',zamp,'fs',Fs,'chan',chan),params);
amp_cluster = sort_spikes(amp_detect);


 [srt,srti] = sort(amp_cluster.cl);
 
 [Tdet,tt] = chopper([-1 1]*.015 ,(amp_detect.spike_indices-1)/Fs,Fs,length(zamp));

view = length(unique(Tdet(:)))/length(zamp);
%  Tdet = mod(Tdet-1,length(zamp))+1;
qdet=zeros(size(zamp));
qdet(amp_detect.spike_indices) = amp_cluster.cl;
detindx = zeros(size(zamp));
detindx(amp_detect.spike_indices)=(1:length(amp_detect.spike_indices));
truindx = zeros(size(zamp));
truindx(round([target.spike_times]*Fs+1))=(1:length([target.spike_times]));

% qtru = zeros(size(zamp));
% qtru(round([target.spike_times]*Fs)) =1;

%  dbq = dbt(detrend([qdet',qtru'],0),Fs,10,'centerDC',false,'upsampleFx',4);
%  xc = ifft(mean(dbq.blrep(:,:,1).*conj(dbq.blrep(:,:,2))));
 
% [Ttru,tt] = chopper([-1 1]*params.tolerance,[target.spike_times],Fs,length(zamp));%-mean(ttav),Fs);
   
 for k = 1:length(target)
    [Ttru,tt] = chopper([-1 1]*params.tolerance,[target(k).spike_times],Fs,length(zamp));%-mean(ttav),Fs);
    
    evtt = tt.*(qdet(Ttru)>0) + 0./qdet(Ttru);
    [mnt,mni] = min(abs(evtt));
    
    target(k).clust = [0 1]'*sum(((1:length(tt))'==mni).*qdet(Ttru)); 
    target(k).clust(1,:) = k;
    target(k).detindx = sum(((1:length(tt))'==mni).*detindx(Ttru));
    target(k).truindx = sum((abs(tt)==min(abs(tt))).*truindx(Ttru));
  
 end
 cl = [target.clust];
 qtru = zeros(size(zamp));
 qtru(round([target.spike_times]*Fs)+1) = 1;
%  Q = (qtru(Tdet)>0).*hann(length(tt)).^4;
% ttav = (tt'*Q)./sum(Q);
% qt = sum(qtru(Tdet).*Q)./sum(Q);
 for k = 1:length(amp_cluster.spike_sep)
    Tdet = chopper([-1 1]*params.tolerance,[ amp_cluster.spike_sep(k).spike_times],Fs,length(zamp));
    
    evtt = tt.*(qtru(Tdet)>0);
    [mnt,mni] = min(abs(evtt) + 0./evtt);
   
   
    amp_cluster.spike_sep(k).clust = sum(((1:length(tt))'==mni).*qtru(Tdet)); 
    amp_cluster.spike_sep(k).clust(2,:) = k;
    amp_cluster.spike_sep(k).truindx = sum(((1:length(tt))'==mni).*truindx(Tdet));
    amp_cluster.spike_sep(k).detindx = detindx(Tdet(abs(tt)==min(abs(tt)),:))';
   
 end
% cl = [ amp_cluster.spike_sep.clust];
cl = [target.clust,  amp_cluster.spike_sep.clust];
inds = [[target.detindx;target.truindx],[amp_cluster.spike_sep.detindx;amp_cluster.spike_sep.truindx]];

%%% Make sure there is no double counting.
%%%First, make sure not to discard any adc events.
[unq,unqi] = unique(inds([2 1],:)','rows');
cl = cl(:,unqi);
inds = inds(:,unqi);


[unq1,~,unqi1] = unique(inds(1,:)');
mult1 = find(crosstab(unqi1)>1 & unq1~=0);
multi1=find(ismember(unqi1,mult1));

[unq2,~,unqi2] = unique(inds(2,:)');
mult2 = find(crosstab(unqi2)>1 & unq2~=0);
multi=find(ismember(unqi1,mult1) | ismember(unqi2,mult2));
minds = inds(:,multi);

dscminds = find(diff([0,minds(2,:)])==0|minds(2,:)==0);
inds(:,multi(dscminds)) = [];
cl(:,multi(dscminds)) = [];




[a,b] = ndgrid(0:length(target),0:amp_cluster.Nclust);
clpad = [cl,[a(:)';b(:)']]; %Padding to make sure no cells are omitted
ctab = crosstab(clpad(1,:),clpad(2,:))-1; %Effect of padding is removed by subtracting 1
 
 out.ctab = ctab;
 out.allhits = ctab(2:end,:)*[0,ones(1,amp_cluster.Nclust)]';
 out.targetN = [0 ones(1,length(target))]*ctab*[1,ones(1,amp_cluster.Nclust)]';
 out.hitrate = sum(out.allhits)./out.targetN;
 out.detectN = [1 ones(1,length(target))]*ctab(:,2:end);
 out.precision =  sum(sum(ctab(2:end,2:end)))./sum(sum(ctab(:,2:end)));
 out.classhitrate = sum(ctab(2:end,2:end),1)./sum(sum(ctab(2:end,:)));
 out.classprecision = sum(ctab(2:end,2:end),1)./sum(ctab(:,2:end),1);
 out.view = view;
 out.SNRamp = SNRamp;
 out.SNRsd = SNRsd;
 out.F1 = 2*out.hitrate.*out.precision./(out.hitrate+out.precision);
 out.FMindx = sqrt(out.hitrate*out.precision);
 out.classF1 = 2*out.classhitrate.*out.classprecision./(out.classhitrate+out.classprecision);
 out.classFM = sqrt(out.classhitrate.*out.classprecision);
 
 