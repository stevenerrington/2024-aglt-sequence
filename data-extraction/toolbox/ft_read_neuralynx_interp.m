function [data_all] = ft_read_neuralynx_interp(fname)

% FT_READ_NEURALYNX_INTERP reads a cell-array of NCS files and performs interpolation
% and gap-filling of the timestamp axis to correct for potentially different offsets
% across channels, and potential gaps within the recordings. All samples are being
% read out.
%
% Input
%  input FNAME is the list of CSC filenames
%  All of these channels should have the sample sampling frequency.
% Output
%  data_all is a raw data structure containing interpolated data and NaNs at the gaps,
%  based on all the available samples in a recording.
%
% Copyright (c) Martin Vinck, SILS, Center for Neuroscience, University of Amsterdam, 2013.


cfg.showcallinfo = 'no';

% first check if these are indeed ncs files
ftype = zeros(length(fname), 1);
for i=1:length(fname)
if     ft_filetype(fname{i}, 'neuralynx_ncs')
  ftype(i) = 1;
end
end
if ~all(ftype==1), error('some files do not correspond to ncs files'); end

% number of channels
nchans = length(fname);

% get the original headers
for i=1:nchans
orig(i) = ft_read_header(fname{i});
end

% check if they all have the same sampling frequency, otherwise return error
for i=1:length(orig), SamplingFrequency(i) = orig(i).Fs; end
if any(SamplingFrequency~=SamplingFrequency(1))
error('not all channels have the same sampling rate');
end
%%
% get the minimum and maximum timestamp across all channels
for i = 1:nchans
ts    = ft_read_data(fname{i}, 'timestamp', true);
mn(i) = ts(1);
mx(i) = ts(end);
ts1   = ts(1);
md(i) = mode(diff(double(ts-ts1)));

% get the minimum and maximum across all channels
if i>1 && mn(i)<min_all
  min_all = mn(i);
else
  min_all = mn(i);
end

if i>1 && mx(i)>max_all
  max_all = mx(i);
else
  max_all = mx(i);
end
end
%%
% issue some warning if channels don't start or end at the same time
startflag = 0; endflag = 0;
if any(mn~=mn(1)),
warning('not all continuous channels start at the same time');
startflag = 1;
end
if any(mx~=mx(1)),
endflag = 1;
warning('not all continuous channels end at the same time');
end
if any(md~=md(1)), warning('not all channels have same mode'); end

%% take the mode of the modes and construct the interpolation axis
% the interpolation axis should be casted in doubles
mode_dts  = mode(md);
rng       = double(max_all-min_all); % this is small num, can be double
offset    = double(mn-min_all); % store the offset per channel
offsetmx  = double(max_all-mx);
tsinterp  = [0:mode_dts:rng]; % the timestamp interpolation axis

% loop over the channels if the channels have different timestamps
for i = 1:nchans
cfg         = [];
cfg.dataset = fname{i};
data        = ft_preprocessing(cfg);
ts          = ft_read_data(cfg.dataset, 'timestamp', true);

% original timestamaps in doubles, with the minimum ts subtracted
ts          = double(ts-ts(1)) + offset(i);
%%
% check if there are gaps to correct
gaps     = find(diff(ts)>2*mode_dts); % skips at least a sample
if isempty(gaps) && startflag==0 && endflag==0
  fprintf('there are no gaps and all channels start and end at same time, no interpolation performed\n');
else
%   % interpolate the data
%   datinterp   = interp1(ts, data.trial{1}, tsinterp);
% 
%   % you can use NaN to replace the data in the gaps
%   gaps     = find(diff(ts)>2*mode_dts); % skips at least a sample
%   for igap = 1:length(gaps)
%     sel = tsinterp < ts(gaps(igap)+1) & tsinterp > ts(gaps(igap));
%     datinterp(sel) = NaN;
%   end
% 
%   % set data at the end and beginning to nans
%   if startflag==1
%     n = floor(offset(i)/mode_dts);
%     if n>0
%       datinterp(1:n) = NaN;
%     end
%   end
% 
%   % set data at the end and beginning to nans
%   if endflag==1
%     n = floor(offsetmx(i)/mode_dts);
%     if n>0
%       datinterp(end-n+1:end) = NaN;
%     end
%   end
% 
%   % update the FieldTrip data structure
%   data.trial{1} = datinterp; clear datinterp
%   data.time{1}  = [0:length(tsinterp)-1].*(1./data.hdr.Fs);
end

% append all the data
if i==1
  data_all = data;
  clear data
else
  data_all      = ft_appenddata([], data_all, data);
  clear data
end
end

% correct the header information and the sampling information
data_all.hdr.FirstTimeStamp     = min_all;
data_all.hdr.LastTimeStamp      = uint64(tsinterp(end)) + min_all;
data_all.hdr.TimeStampPerSample = mode_dts;
len = length(tsinterp);
data_all.hdr.nSamples           = len;
data_all.sampleinfo             = [1 len];