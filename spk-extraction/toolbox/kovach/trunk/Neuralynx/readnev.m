function evnt = readnev(fnames,pth)

% Basic script to read  neuralynx event data.
%
% C Kovach 2017
params.header_bytes = 16384;
params.data_precision = 'uchar';
params.samphead_bytes=28;
%params.samptime_precision='uint64';
params.samphead  =   struct('nstx','int16',...
                            'npkt_id','int16',...
                            'npkt_data_size','int16',...
                            'qwTimeStamp','uint64',...
                            'nevent_id','int16',...
                            'nttl','int16',...
                            'ncrc','int16',...
                            'ndummy1','int16',...
                            'ndummy2','int16');%,...
%                             'dnExtra','int32');


persistent path

if nargin < 2
    pth = '';
end


% if nargin > 0 
%     [pth2,fnames] = fileparts(fn);
%     pth = fullfile(pth,pth2);
% end

if ~isempty(pth)
    path = pth;
end


if nargin < 1 || isempty(fnames)
    
    [fnames,pth] = uigetfile({'*.nev','NEV files';'*.*','All files'},'Select NLX Event data file',path,'multiselect','on');
    if isnumeric(fnames)
        return
    else
        path = pth;
    end
    
    if ~iscell(fnames)
        fnames = {fnames};
    end
    
  
else

    if ~iscell(fnames)
        fnames = {fnames};
    end
end
%Valid fieldname characters
fnchars ='abcdefghijklmnopqrstuvwxyz1234567890_';
fnchars = unique([fnchars,upper(fnchars)]);

evnt = struct('evnt',[],'time',[],'name','','header',[]);

%parse sample header
sampheadfn = fieldnames(params.samphead);

for k = 1:length(sampheadfn)
    sampheadbyteN(k) =  cellfun(@str2double,regexp(params.samphead.(sampheadfn{k}),'(\d+)','tokens','once'))/8;
end
params.samphead_bytes = sum(sampheadbyteN);

for k = 1:length(fnames)
  
    fn = fullfile(path,fnames{k});
    fid = fopen(fn);
    txt = fread(fid,16384,'uchar=>char')';
    txt = regexprep(txt,'µ','u');
    flds = regexp(txt,'-([^\s]*)\s*([^\n]*)','tokens');
    flds = cat(1,flds{:})';
    flds(1,:) = cellfun(@(x)x(ismember(x,fnchars)),flds(1,:),'uniformoutput',false);
    flds(2,:) = cellfun(@(x)deblank(x),flds(2,:),'uniformoutput',false);
        
    data.header = struct(flds{:});

    
    params.ndatabytes = 1;
    
    params.nsamp  = (str2num(data.header.RecordSize)-params.samphead_bytes)./params.ndatabytes;
    
    % get the sampleheader_values
    skip = 0;
    databytes = params.ndatabytes*params.nsamp;
    for kk = 1:length(sampheadfn)
        fseek(fid,params.header_bytes+skip,-1);
        data.(sampheadfn{kk}) = fread(fid,params.samphead.(sampheadfn{kk}),databytes+params.samphead_bytes-sampheadbyteN(kk));
        skip = skip+sampheadbyteN(kk);
    end
      % Get the events
    fseek(fid,params.header_bytes+skip,-1);    
     evntdata= fread(fid,sprintf('%i*%s=>char',params.nsamp,params.data_precision),params.samphead_bytes);
     evntdata = reshape(evntdata,databytes,length(evntdata)/databytes)';
     evnt(k).dnExtra = evntdata(:,1:32);
     evnt(k).log = evntdata(:,32+1:end);
     evnt(k).evnt = data.nttl;
     evnt(k).time = data.qwTimeStamp/1e6; 
     evnt(k).time_fmt = 'unix_sec';
     evnt(k).header = data.header;
     [~,evnt(k).name] = fileparts(fn);
     evnt(k).file = fnames{k};
end
fclose(fid);
