function [FIX,RAW] = readEDF(filename,varargin)

% [FIX,RAW] = readEDF(filename,varargin)
%
%Uses edfmex to read edf and puts the data into a matlab data structure. 
%
%
% FIX is a data structure with the following fields:
%
% FIX.recdata   : recording parameter data. See eyelink manual for details.
% FIX.eyetracker: A label for the eye tracker used in the experiment.
% FIX.header    : EDF file header, which contains such as recording date, eye
%           tracker model, camera type, etc.
% FIX.filename  : name of the file from which the data were extracted.
% FIX.seg       : A structure array for each recording segment with the fields:
%        .xdat      : a structure array containing message events with fields:
%             .startT   : Time message was received
%             .id       : a numeric identifier for the message.
%             .code: the string which was received.
%        .fix       : structure array containing data for each fixation:
%            .meanPos   : average gaze position during fixation.
%            .startT    : fixation onset time.
%            .dur       : fixation duration.
%            .endT      : fixation end time.
%            .updconv   : conversion to degrees from screen units, given by the eyetracker.
%            .eye       : which was recorded - 1 left, 2 - right.
%            .xdhist    : numeric id for the last two message events [earlier, later]
%            .xdat      : most recent message event.
%            .xdindex   : ids for two most recent message events (indexes into .xdatCodes)
%            .lastcode  : string for the most recently received message.
%            .shiftvec  : difference between meanPos for the previous fixation and the current one.
%            .dt: time between the previous fixation and the present one.
%            .sac       : index into the .sac array for the preceding saccade.
%         .sac      : structure array containing data for each saccade:
%            .startT    : start of saccade.
%            .dur       : duration of saccade.
%            .startPos  : gaze position at start of saccade.
%            .shiftVec  : difference between findal and starting position.
%            .eye       : eye for the saccade event.
%            .lastfix   : index into .fix array for the last fixation.
%            .nextfix   : index into .fix array for the next fixation.
%            .xdhist    : codes for two most recent message events (indexes into .xdatCodes)
%            .xdindex   : index into .xdat for the most recent message event.
%            .lastcode  : string for the most recently received code.
%         .xdatCodes: a cell array of all unique messages received during
%                    the recording segment. The .xdhist fields in .fix and .sac 
%                    and the id field in .xdat index into this array.
% FIX.units     : label for position units. 
% FIX.fs        : sampling frequency.
%
%
% RAW is a data structure with the following fields:
% 
% RAW.recdata   : same as in FIX.
% RAW.eyetracker:   "
% RAW.header    :   "
% RAW.filename  :   "
% RAW.seg       : array structure for recording segments with the following fields:
%        .horz      : horizontal gaze position.
%        .vert      : vertical gaze position.
%        .pupil     : pupil diameter.
%        .input     : external input, such as parallel port.
%        .xdat      : the id for the most recently received message at each time
%               point.
%        .degConversion  : horz,vert units-to-degrees conversion as reported
%                       by eyelink.
%        .sample    : time in ms.
%        .xdatCodes : same as in FIX.seg.
% RAW.fs        : sampling frequency.
%
%
% See also EDFMEX



% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C Kovach 2010

nhist = 2;
i = 1;
subtractFirstXdat = false;
eyeinp = 'x';
while i <= length(varargin)
   switch lower(varargin{i})
       case 'subtract_first_xdat_time' %Subtract the time of the first message event. 
           subtractFirstXdat = true;
       case 'eye' %'Use EYE: 'L' for left,'R' for Right, A for average
           eyeinp = varargin{i+1};
           i = i+1;
       case 'nhist' %'Use EYE: 'L' for left,'R' for Right
           eyeinp = varargin{i+1};
           i = i+1;
        otherwise
           error([varargin{i},' is not a valid option.']);
   end         
    i = i+1;
end

% if nargin < 2 || isempty(nhist)
%     nhist = 2;
% end

if ~exist(filename,'file')
    error('Unable to find: %s',filename);
end

try 
    fprintf('%s ',filename)
    fprintf('\n')
    EDF = edfmex(filename);
    fprintf('\n')
catch
    le = lasterror;
    error(sprintf(['\nIt appears the pre-compiled file edfmex.mex is broken or incompatible with your system.',...
            '\nTry recompiling edfmex.cpp in the source directory with mex ''edfmex.cpp edfapi.lib''']))
end

RAW.recdata = EDF.RECORDINGS;
RAW.eyetracker = 'Eyelink';
RAW.header = EDF.HEADER;
RAW.filename = EDF.FILENAME;
FIX.recdata = EDF.RECORDINGS;
FIX.eyetracker = 'Eyelink';
FIX.header = EDF.HEADER;
FIX.filename = EDF.FILENAME;

eye = EDF.RECORDINGS(1).eye;
etfs = EDF.RECORDINGS(1).sample_rate;
if eye == 3
%     inp = 'x';
    while ~ismember(eyeinp,{'r','l','b','a'})
        eyeinp = lower(inputdlg('Which eye to load: (L)eft, (R)ight, (B)oth or (A)Average?'));
    end
    
    eye = 3*isequal(eyeinp{1},'b')+ 2*isequal(eyeinp{1},'l')+isequal(eyeinp{1},'r');
    
    if eye == 0  
       eye  = [1 2];
       avgeye=true;
    elseif eye == 3
        eye = [1 2];
        avgeye=false;
    else
        avgeye=false;
        
    end
else
    avgeye=false;

end

    
code = {EDF.FEVENT.codestring};
msg = strcmp(code,'MESSAGEEVENT');
msgevnt = EDF.FEVENT(msg);
msgcodes = {msgevnt.message};
msgt = [msgevnt.sttime] - EDF.FSAMPLE.time(1)+1;

[unqmsg,q,msgmap] = unique(msgcodes);

msgsig  = zeros(size(EDF.FSAMPLE.time));

msgsig(ceil(single(msgt)*etfs./1e3)) = sum(repmat(msgt,length(msgt),1) == repmat(msgt',1,length(msgt))) ; 
xdatraw = msgmap(cumsum(msgsig));

if subtractFirstXdat
    startTime = msgevnt(1).sttime;
else
    startTime = 0;
end

if nargout > 1
    RAW.seg.horz = double(EDF.FSAMPLE.gx(eye,:));
    RAW.seg.vert = double(EDF.FSAMPLE.gy(eye,:));
    RAW.seg.pupil = double(EDF.FSAMPLE.pa(eye,:));
    if eye == 0
        RAW.seg.horz =mean(RAW.seg.horz );
        RAW.seg.vert = mean(RAW.seg.vert);
        RAW.seg.pupil = mean(RAW.seg.pupil);        
    end
    RAW.seg.input = double(EDF.FSAMPLE.input); % Parallel port input
    RAW.seg.xdat = xdatraw;
    RAW.seg.degConversion(:,1) = double(EDF.FSAMPLE.rx);
    RAW.seg.degConversion(:,2) = double(EDF.FSAMPLE.ry);
    RAW.seg.time = double(EDF.FSAMPLE.time) - double(startTime) ;
    RAW.seg.xdatCodes = unqmsg; %1:length(unqmsg);
%     RAW.seg.xdatLabels = unqmsg;
end

%fixation data
endfix = strcmp(code,'ENDFIX');
endsacc = strcmp(code,'ENDSACC');



fxs = EDF.FEVENT(endfix);

if length(msgevnt) ~= length(msgcodes)
    warning('Uh Oh. Some of the messages fields are empty.\nEverything will be out of alignment if I continue.')        
    inp = input(sprintf('1. Stop\n2. Continue\n: '));
    if inp == 1
         error('Stopping...');
    end

end


FIX.seg.xdat= struct('startT',{},'id',{},'code',{});
for i = 1:length(msgmap)
    FIX.seg.xdat(i).startT = double(msgevnt(i).sttime) - double(startTime);
    FIX.seg.xdat(i).id = msgmap(i);
    FIX.seg.xdat(i).code= msgcodes{i};    
end
msgts = [FIX.seg.xdat.startT];

inpt = diff(EDF.FSAMPLE.input);
inpti = find(inpt);
inptT = EDF.FSAMPLE.time(inpti);
FIX.seg.input = struct('startT',{},'id',{});
for i = 1:length(inptT)
     FIX.seg.input(i).startT = double(inptT(i)) - double(startTime);    
     FIX.seg.input(i).id = double(inpt(inpti(i)));     
end

saccs = EDF.FEVENT(endsacc);
if ~isempty(saccs)
    saceye = double([saccs.eye]);
    saccs = saccs(ismember(saceye,(eye - 1)));
    sact = [saccs.sttime]-startTime;
end

    
fxeye=double([fxs.eye]);

if sum(eye) == 3
    
    if avgeye
        tthresh = 100; %% combine fixations only if they occur within this interval
        dthresh = 5; %%% threshold normalized intrafixation variance
        tthresh_keep = 200; %%% Keep monocular fixations if they are outside this time range

        fxt=double([fxs.sttime]);

        nfx = length(fxs);

        for k = 1:nfx
            [mn,mni] = min(abs( fxt(k) - fxt + 1e10*(fxt(k) - fxt==0 )));
            assoc(k) = mni;
            kp(k)=(mn <= tthresh) & (fxeye(k)~=fxeye(mni)); 
            kpmonoc(k)= (mn > tthresh_keep); 
        end   

        spd = spalloc(length(assoc),length(assoc),2*numel(assoc));
        spd( sub2ind(size(spd),find(assoc),assoc(assoc~=0)) )= 1;

        sym = spd.*spd'; %% ignore fixations that aren't reciprocally associated 
        kp = kp & any(sym);

        fxx=double([fxs.gavx]);
        fxy=double([fxs.gavy]);

        stfx= double([fxs.gstx]);
        stfy= double([fxs.gsty]);
        enfx= double([fxs.genx]);
        enfy= double([fxs.geny]);

        %%% within fixation varance 
        fvar = mean( (stfx/2+enfx/2-fxx).^2 + (stfy/2+enfy/2-fxy).^2);

        %%%% inter eye distance
        dfx=0*kp;
        dfx(kp) = sqrt(detrend( (fxx(kp)-fxx(assoc(kp))).*(fxeye(kp)-.5)*2).^2 +...
              detrend(( fxy((kp))-fxy(assoc((kp)))).*(fxeye(kp)-.5)*2).^2);

        %%% normalized inter-fixation distance
        ndfx = dfx./sqrt(fvar);
        kp = kp & ndfx<=dthresh;
        assoc(kpmonoc) = find(kpmonoc);
    else
        kpmonoc = true(size(fxs));
        kp = false(size(fxs));
        assoc= 1:length(fxs);
    end
    kpall =kp | kpmonoc;
    fxs=fxs(kpall);

    eyes = false(2,sum(kpall));
    eyes(:,kp(kpall)) = true;
    eyes((find(kpmonoc(kpall))-1)*2 +fxeye(kpmonoc(kpall))+1) = true;
    eyes = eyes./repmat(sum(eyes),2,1);
   
    indxkp(kpall) =1:sum(kpall);
    assoc = indxkp(assoc(kpall));
else
    assoc= 1:length(fxs);
    eyes = ones(size(assoc));
    
end

fixstr = struct('meanPos',[],'startT',[],'dur',[],'endT',[],'updconv',...
                     [],'eye',[],'xdhist',[],'xdat',[],'xdindex',[],'lastcode',[],...
                     'shiftvec',[],'dt',[],'sac',[]);
sacstr = struct('startT',[],'dur',[],'startPos',[],'shiftVec',[],'eye',[],...
                'lastfix',[],'nextfix',[],'xdhist',[],'xdindex',[],'lastcode',[]);

FIX.seg.fix = fixstr;
FIX.seg.sac = sacstr;

i = 0;
for k = 1:length(fxs)
    if ~ismember(double(fxs(k).eye), eye - 1)
        continue
    end
    i = i+1;
    
    if avgeye
        inds = [k assoc(k)];
    else 
        inds = k;
    end
    
    
    FIX.seg.fix(i).meanPos = eyes(:,k)'*[fxs(inds).gavx; fxs(inds).gavy]';
    FIX.seg.fix(i).startT = double([fxs(inds).sttime]-startTime)*eyes(:,k);
    FIX.seg.fix(i).dur   = double([fxs(inds).entime] -[fxs(inds).sttime])*eyes(:,k); 
    FIX.seg.fix(i).endT = double([fxs(inds).entime]-startTime)*eyes(:,k);
    FIX.seg.fix(i).updconv = eyes(:,k)'*double([fxs(inds).supd_x; fxs(inds).supd_y]' + [fxs(inds).eupd_x; fxs(inds).eupd_y]')./2;
    FIX.seg.fix(i).eye = [fxs(inds).eye];

    recentmsg = msgts <= double([fxs(inds).sttime]-startTime)*eyes(:,k);
    frecentmsg= find(recentmsg );
    FIX.seg.fix(i).xdhist  = zeros(1,nhist);
    FIX.seg.fix(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    FIX.seg.fix(i).xdat = FIX.seg.fix(i).xdhist(end);
    cr = cumsum(recentmsg);
    FIX.seg.fix(i).xdindex = cr(end);    
    FIX.seg.fix(i).lastcode = msgcodes{frecentmsg(end)};
    
    if i>1
        FIX.seg.fix(i).shiftvec = diff(cat(1,FIX.seg.fix(i + (-1:0)).meanPos));
        FIX.seg.fix(i).dt = diff(cat(1,FIX.seg.fix(i + (-1:0)).startT)); % difference between fixation onsets of 
    
                                                                                         %the current and previous fixation
    else
        FIX.seg.fix(i).shiftvec =[nan nan];
        FIX.seg.fix(i).dt =nan;
    end
    
    FIX.seg.fix(i).sac = find( sact <= FIX.seg.fix(i).startT,1,'last');
    if isempty(FIX.seg.fix(i).sac)
        FIX.seg.fix(i).sac = nan;
    end
     FIX.units = 'xy Screen';        
end


i = 0;
fxts = [FIX.seg.fix.startT];
for k = 1:length(saccs)
    if ~ismember(double(saccs(k).eye), eye-1)
        continue
    end
    i = i+1;
    
    stt = double(saccs(k).sttime)- double(startTime);
    dfst = fxts - stt;
    lastfix = find(dfst < 0,1,'last'); 
    nextfix = find(dfst > 0,1,'first'); 
    
    FIX.seg.sac(i).startT= stt;
    FIX.seg.sac(i).dur= saccs(k).entime-saccs(k).sttime;
    FIX.seg.sac(i).startPos= [saccs(k).gstx saccs(k).gsty];
%     FIX.seg.sac(i).shiftVec= [saccs(k).gstx saccs(k).gsty] - [saccs(k).genx saccs(k).geny];    
    FIX.seg.sac(i).shiftVec= [saccs(k).genx saccs(k).geny] - [saccs(k).gstx saccs(k).gsty];    
    FIX.seg.sac(i).eye= saccs(k).eye;
    if ~isempty(lastfix)
        FIX.seg.sac(i).lastfix= lastfix;
    else
        FIX.seg.sac(i).lastfix= nan;        
    end
    if ~isempty(nextfix)
        FIX.seg.sac(i).nextfix= nextfix;
    else
        FIX.seg.sac(i).nextfix= nan;        
    end
    recentmsg = msgts <= saccs(k).sttime;
    frecentmsg = find(msgts <= saccs(k).sttime);
    FIX.seg.sac(i).xdhist  = zeros(1,nhist);
    FIX.seg.sac(i).xdhist(end:-1:end-min([length(frecentmsg) nhist]-1)) = msgmap(frecentmsg(end:-1:end-min([length(frecentmsg) nhist]-1)));
    cr = cumsum(recentmsg);
    FIX.seg.sac(i).xdindex = cr(end);
    FIX.seg.sac(i).lastcode = msgcodes{frecentmsg(end)};
    
end
FIX.seg.xdatCodes = unqmsg;
RAW.fs = etfs;
FIX.fs = etfs;
