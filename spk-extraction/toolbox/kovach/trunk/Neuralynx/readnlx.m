function [data,evnt] = readnlx(fnames,pth,load_directory,chan_types,varargin)

% Basic script to read neuralynx continuous data (*.ncs) and event (*.nev) files.
%
% Usage:
%
%  data = readnlx
%
%      With no input arguments, files are selected through a GUI.
%
%  data = readnlx(filename,[path])
%
%       Filename can be either a string or cell array of filenames.
%
%  data = readnlx(filename,[path],load_directory)
%
%       Load all files in the directory if true
%
%  data = readnlx(filename,[path],true,chan_types)
%
%       Load all files in the directory of channel_type, where channel_type
%       is a cell array {'TYPEA','TYPEB',...}. The ncs files are assumed to
%       be prepended with 'TYPEX', etc. 
%
%  [cdata,evnt] = readnlx({'file1.ncs','event1.nev'},...)
%
%       If both continuous and event data are selected, continuous data are
%       returned in the first output argument and events in the second.  
%
% C Kovach 2017

persistent path

if nargin < 2
    pth = '';
end

if nargin < 3 || isempty(load_directory)
    load_directory=false;
end

if nargin < 4 || isempty(chan_types)
    chan_types = {''};
    load_directory=false;
else
    if ischar(chan_types)
        chan_types = {chan_types};
    end
end

if ~isempty(pth) 
    path = pth;
end

if isnumeric(path)
    path='';
end

file_size_filter = 16384; % Ignore files that appear to be empty. 


if nargin < 1 || isempty(fnames)
    [fnames,sz] =  get_nlx_files_sorted(load_directory,pth,chan_types);
    fnames = fnames(sz>file_size_filter | isnan(sz));

else

    if ~iscell(fnames)
        fnames = {fnames};
    end
end
kd = 0;
ke = 0;

%% Remove PDes97 onwards
pattern = 'PDes(97|98|99|\d{3,})'; % Matches 'PDes97' or any number greater than or equal to 97
indicesToRemove = cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), fnames);
fnames(indicesToRemove) = [];
%%

for k = 1:length(fnames)
    [pth,fn,ext] = fileparts(fnames{k});  
    if ~isempty(pth)
        path = pth;
    end
    fn = [fn,ext];
    switch lower(ext(2:end))
    
        case 'ncs'
            kd = kd+1;
            data(kd) = readncs(fn,path,varargin{:});
        case 'nev'
            ke = ke+1;
            evnt(ke) = readnev(fn,path);
     end
                
end

if ke==0
    evnt=[];
end

if kd == 0 && nargout == 1
    data = evnt;
end
