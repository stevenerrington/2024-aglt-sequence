

Fs = 20e3;

fid = fopen('yger_data/patch_2_MEA.juxta.raw');
jux = fread(fid,'single=>double');
fclose(fid);

[jux_detect,jux_hos,zjux] = HOSD_spike_detection(struct('dat',jux,'fs',Fs),struct('ncomp',3,'windur',.01));

jux_cluster = sort_spikes(jux_detect);


fid = fopen('yger_data/patch_2_MEA.raw');
amp = fread(fid,'uint16');
recast = @(amp) mod(double(amp)+double(intmax('uint16')/2),double(intmax('uint16')))-double(intmax('uint16')/2);
% amp = mod(double(amp)+double(intmax('uint16')),intmax('uint16'))-double(intmax('uint16'));
fclose(fid);
amp = reshape(amp(1:6e6*256),256,6e6)';

zamp = zscore(double(recast(amp(:,166))));

[amp_detect,hosamp,zamp] = HOSD_spike_detection(struct('dat',zamp,'fs',Fs,'chan',8),struct('ncomp',5,'lowpass',4000,'windur',.015,'highpass',100));
amp_cluster = sort_spikes(amp_detect);

qjux = zeros(size(zamp));
qjux(jux_detect.spike_indices) = jux_cluster.cl;

qamp = zeros(size(zamp));
qamp(amp_detect.spike_indices) = amp_cluster.cl;

[Tjux] = chopper([-1 1]*.015 ,jux_detect.spike_indices/Fs,Fs,length(zamp));
[Tamp,tt] = chopper([-1 1]*.015 ,amp_detect.spike_indices/Fs,Fs,length(zamp));


