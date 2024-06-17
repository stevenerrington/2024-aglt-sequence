
function   initialize_stimulator(obj,stimstruct,getinput)



% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/TDT/@stimulator/initialize_stimulator.m $
% $Revision: 84 $
% $Date: 2012-03-01 18:42:19 -0600 (Thu, 01 Mar 2012) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    getinput = false;
end
if nargin <2 || isempty(stimstruct)
    stimstruct=obj;
end
if nargin <2 || getinput

    quest = {  'Pos. polarity channels: ', num2str(stimstruct.pos_chan)
               'Neg. polarity channel: ',  num2str(stimstruct.neg_chan)
               'Duration (s): ',  num2str(stimstruct.duration)
               'Level (uA/V):',  num2str(stimstruct.level)
               'Mode (const. curr. = 0; const. volt. = 1): ', num2str(stimstruct.stim_mode)
               'Frequency (Hz): ',  num2str(stimstruct.frequency) 
               'Do stimulation (1=yes,0=no): ', num2str(stimstruct.do_stimulate)};

    fields = {'pos_chan', 'neg_chan', 'duration','level','stim_mode','frequency','do_stimulate'}';   

    inp =inputdlg(quest(:,1),'Set and verify stimulation parameters',1,quest(:,2));

    ninp = cellfun(@(x)str2num(x),inp,'uniformoutput',0);
    ninp([5 7])=inp([5 7]);
    

    for i = 1:length(fields)
        obj.(fields{i}) = ninp{i};
    end

else
    fields = fieldnames(stimstruct);
    
    for i = 1:length(fields)
        obj.(fields{i}) = stimstruct.(fields{i});
    end
end

