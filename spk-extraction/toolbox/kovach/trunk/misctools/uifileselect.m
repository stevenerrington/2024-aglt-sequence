function fn = uifileselect(varargin)

% filename = uifileselect
%
% File selection dialog that does not throw a warning if a file
% exists like uiputfile.
% 
% filename = uiselectfile(filetypes,starting_path,title)
% 
% Opens in starting_path and filters by files with extension 'filetype'. 
%
% The filter syntax is the same as for uigetfile. 
%
% The title for the dialog box is specified with 'title'.
%
% See also UIGETFILE and UIPUTFILE
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/misctools/uifileselect.m $
% $Revision: 151 $
% $Date: 2013-01-18 02:34:05 -0600 (Fri, 18 Jan 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2013

persistent fpath

if nargin < 1 || isempty(varargin{1})
%     ffilt  = {'*.kov','Kovach Data Format (*.kov)'; '*.*','All Files'};
      ffilt  = {};
else
    ffilt = varargin{1};
    if ischar(ffilt)
        ffilt = {ffilt};
    end
end

if nargin < 2 || isempty(varargin{2})
    if isempty(fpath)
        fpath = cd;
    end
else
    fpath = varargin{2};
end

if nargin < 3 || isempty(varargin{3})
    title = 'Which file?';
else
    title = varargin{3};
end

parframe = com.mathworks.hg.peer.utils.DialogUtilities.createParentWindow;
assert(~isempty(parframe));

obj = javahandle_withcallbacks.com.mathworks.mwswing.MJFileChooser;

if ~isempty(ffilt)
    javaFileExtensionFilters = getPeer(matlab.ui.internal.dialog.FileExtensionFilter(ffilt));
    for i = 1:length(javaFileExtensionFilters)
        obj.addChoosableFileFilter(javaFileExtensionFilters{i});
    end
end

obj.setCurrentDirectory(java.io.File([fpath,filesep]));

obj.setDialogTitle(title);

obj.showOpenDialog(parframe);

fn = char(obj.getSelectedFile);

fpath = char(obj.getCurrentDirectory);

delete(obj)
