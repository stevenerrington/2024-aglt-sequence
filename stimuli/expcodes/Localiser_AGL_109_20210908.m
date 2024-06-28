function Localiser_AGL_109(subInit,runNo)
% FixBasic(subInit,runNo,hand,eyeTrackingOn)
% this is for opto + faces/voices
% - SubInit - subject initials in string form 'example'
% Following parameters not required.
% - runNo - which run of the session is this?
% - hand - (1) right, (2) left
% - eyeTrackingOn - (1) to record eye tracking, (0) to not record it

%Last update: 2021-June-2 (YK)

%% Initial checks
clearvars -except subInit runNo eyeTrackingOn par.eyeTrackAbort;

%% juice code for convenience here
% ARG.labjack_ver = 6; % Set for U3 (3) or U6 (6)
% intialiseLabjack; % run subscript to initialise LabJack
% Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,5,0); % juice ON
% Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0); % juice OFF

%% EXPERIMENT DETAILS
% Update these for each experiment.
AGL_type = 2; %1: Exposure, 2: Test
%%% Experiment name %%%
par.Expt = 'AGL'; % for your information%% Session parameters.

%%% Number of Reps of all trials %%%
% Check input arguments
if nargin == 0
    %     help EYE1_LJ_nonMulti_preAllocate
    return
elseif nargin == 1
    runNo = 0;
    par.eyeTrackAbort = 1;
elseif nargin == 2
    par.eyeTrackAbort = 1;
    subInit='test';
end

clc;
disp('NOTE: If errors try running ''PsychJavaTrouble'' without any arguments');

rand('state',sum(100*clock)); % seed random number generator

%% Set directories

parentdir = '\\campus\rdw\ion03\03\lcn\EPHYS\RAWDATA\NHP\Neuralynx\AGL\Troy\';
scriptsDir = [parentdir 'scripts\'];
stimDir = [parentdir 'stim\'];
ddir = [parentdir 'data\'];
helpersDir = 'C:\Users\nlcn109\MATLAB\Optogenetics_109\helpers';

if AGL_type == 1
    sdir =  ['C:\Users\nlcn109\MATLAB\stimuli\AGL\stimuli\Exposure'];
    par.numReps = 20; %Change to 20
elseif AGL_type == 2
    sdir =  ['C:\Users\nlcn109\MATLAB\stimuli\AGL\stimuli\SingleMulti'];
    par.numReps = 10; 
elseif AGL_type == 2
    sdir =  ['C:\Users\nlcn109\MATLAB\stimuli\AGL\stimuli\SingleMulti'];
    par.numReps = 10;     
    
    
end

addpath('C:\Users\nlcn109\MATLAB\INSTALL\ViewPoint\2.9.5.117\2.9.5.117\Burn to disk\ViewPoint\Interfaces\Windows\MATLAB\ViewPoint_EyeTracker_Toolbox');
addpath(helpersDir);
addpath('C:\Users\nlcn109\MATLAB\INSTALL\LabJackMatlab\MATLAB_LJUD\LJUD_Functions');

if ~exist(ddir,'dir'),
    mkdir(ddir);
end

%% Experiment Parameters - DO NOT CHANGE DURING EXPERIMENTS
[~, machine] = system('hostname'); % get machine name
machine = strtrim(machine); % remove redundant blank characters

ARG.lab = '109';
ARG.getKeyResponses = 0;
ARG.MRI = 0;
ARG.input = 0; % 0 = using mouse (testing); 1 = Using touch lever (111)

%%% Monitor display %%%
ARG.whichScreenFig = 2; % 1 = Primary, 2 = secondary - which screen to display eye trace fig
ARG.screenNumberExperiment = 1; % 1 = Primary, 2 = secondary - which to show experiment

%%% SET LABJACK VERSION %%%
ARG.labjack_ver = 6; % Set for U3 (3) or U6 (6)


%% %% Timing details, all in milliseconds %%%
par.fixOn       = 10;       % when fixation starts
par.fixRequired = 500;      % when fixation required %% IF NEEEDED CAN CODE TO WAIT FOR HIM TO FIXATE TO START TRIAL - CIP
par.stimOn      = 1000;     % when stimOn, stimuli less that 1 second; !! Best not to do anything else during stimulus presentation !!
par.fixDur      = par.stimOn + 2667;     %% when to turn off fixation [413 976 1540 2103 2667] (NEED TO CHANGE EVERY TRIAL DEPENDING ON SEQ LENGTH -YK)
par.rewardDelay = par.fixDur + 500;     % when to reward
par.rewardDur   = 300;      % how long of a reward (from rewardDelay time)
par.trialEnd    = par.rewardDelay + par.rewardDur + 500;     % when does trial end

par.abortTime = 1000;       % timeout after abort

%% EVENT TIMING %%%
% Event                            % Time of onset                                             % Flag
events(1,:) = [{'trialStart'},     0,                                                            0];
events(2,:) = [{'fixSpotOn'},      par.fixOn,                                                    0];
events(3,:) = [{'fixRequiredOn'},  par.fixRequired,                                              0]; % 700
events(4,:) = [{'stimPlay'},       par.stimOn,                                                   0]; %another way to grab time: events{strcmp(events(:,1),'fixSpotOn'),2}+par.fixDur
%events(4,:) = [{'optoOn'},      600,                                                             0]; %another way to grab time: events{strcmp(events(:,1),'fixSpotOn'),2}+par.fixDur
events(5,:) = [{'fixRequiredOff'}, par.fixDur                                                    0];
events(6,:) = [{'fixSpotOff'},     par.fixDur                                                    0];
events(7,:) = [{'rewardOn'},       par.rewardDelay,                                                 0];
events(8,:) = [{'rewardOff'},      par.rewardDelay + par.rewardDur,                                                0];
events(9,:) = [{'trialEnd'},       par.trialEnd,             0]; % Trial end is hardcoded!!!!!
%% Display a figure showing the event timings, to check that everything is correct.
trialFigure(events)

%%
%%% Eye Tracking %%%
eyeTrackingOn = 1;
par.eyeTrackAbort = 1;      % use eyetracking?
par.eyeSampleRate = 0.01;   % How often to sample the eye tracking position (secs)
par.blinkAllowance = 4;    % blink allowance in samples (for now)
par.fixSize = 0.1;          % Size of fix window
par.centre = [0.5,0.5];     % Centre location
par.minVal = par.centre(1) - par.fixSize; % used to calculate quality of fixation
par.maxVal = par.centre(1) + par.fixSize; % used to calculate quality of fixation
par.spotSize = 40;          % Initial size of fixaiton spot
par.xAdjust = 0;            % adjustment to manually calibrate eye traces (move manual slip correction)
par.yAdjust = 0;
par.xOffset = 0;            % adjust position of fixation window based on the position of the fixation spot
par.yOffset = 0;
par.reviewEyeTraces = 1;    % subplots at end of experiment: 0 - no eye trace plots, 1 - all trials plotted, % 2 = every 10th trial, 3 = every 20th trial
par.zoom = 1;


% disp('are these trial timings correct?')
% pause

%% INITALISE REQUIRED HARDWARE AND LIBRARIES
% Initalise Eyetacking toolbox
% if ARG.MRI == 0
vpx_Initialize; % initialise the eye tracking toobox
% end

% Run labjack setup
intialiseLabjack; % run subscript to initialise LabJack

% Stop immediate juice delivery
pause(0.05);
ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0); % juice OFF

%% Set up datafile
datetime = datestr(now);
datetime = strrep(datetime,':','_'); % Replace colon with underscore
datetime = strrep(datetime,'-','_'); % Replace minus sign with underscore
datetime = strrep(datetime,' ','_'); % Replace space with underscore

if runNo == 0
    datafile = [subInit,'_',datetime]; % int2str for subInit
else
    datafile = [subInit,'_Run#',int2str(runNo),'_',datetime]; % int2str for subInit
end

%% set up dat structure
dat = [];
dat.Expt = par.Expt;
dat.sub  = subInit;
dat.events = events(:,1:2);
dat.run  = runNo;
dat.par = par;
dat.ARG = ARG;

%% Pre-allocating
dat.eye = [];
dat.eye.eye = 0;
dat.eye.sampling = par.eyeSampleRate; % how often to sample eye tracker


par.num_videostims = 8;

%conditions 1:8 audio only, 7:12 video, 13:18 both

par.condsStr = '1 = Audio; 2 = Video; 3 = AV';


audiostims = [];
% Humans:
% 2 Human: Chris and Christoph
% Monkeys (3 monkeys):
% 2 coos Stim_Emil; Stim_Heini;
% 2 grunts: Heini, Prince

if AGL_type ==1
    par.num_audiostims = 8;
    par.stimNums = [1:8];
    par.conds = [1,1,1,1,1,1,1,1]; % 1 = Audio only; 2 = Video only; 3 = AV
    audiostims{1} = ['Grammatical_2.wav'];
    audiostims{2} = ['Grammatical_3.wav'];
    audiostims{3} = ['Grammatical_4.wav'];
    audiostims{4} = ['Grammatical_8.wav'];
    audiostims{5} = ['Grammatical_6.wav'];
    audiostims{6} = ['Grammatical_7.wav'];
    audiostims{7} = ['Grammatical_9.wav'];
    audiostims{8} = ['Grammatical_11.wav'];
   
elseif AGL_type ==2
    par.num_audiostims = 16;
    par.conds = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; % 1 = Audio only; 2 = Video only; 3 = AV
    par.stimNums = [1:16];
    audiostims{1} = ['Grammatical_3.wav'];
    audiostims{2} = ['Grammatical_8.wav'];
    audiostims{3} = ['Grammatical5_novel.wav'];
    audiostims{4} = ['Grammatical10_novel.wav'];
    
    audiostims{5} = ['Grammatical_3.wav'];
    audiostims{6} = ['Grammatical_8.wav'];
    audiostims{7} = ['Grammatical5_novel.wav'];
    audiostims{8} = ['Grammatical10_novel.wav'];
    
    audiostims{9} = ['Multiviolation_1.wav'];
    audiostims{10} = ['Multiviolation_2.wav'];
    audiostims{11} = ['Multiviolation_3.wav'];
    audiostims{12} = ['Multiviolation_4.wav'];
    
    audiostims{13} = ['Single_Rule_break_1.wav'];
    audiostims{14} = ['Single_Rule_break_2.wav'];
    audiostims{15} = ['Single_Rule_break_3.wav'];
    audiostims{16} = ['Single_Rule_break_4.wav'];
end



% removed dummy sound - cip
movstims = [];
movstims{1} = ['Stim_Chris.avi']; %Human_1 face
movstims{2} = ['Stim_Christoph.avi'];%Human_2 face
movstims{3} = ['Stim_Heini.avi'];%Monkey_1 face
movstims{4} = ['Stim_Emil.avi']; %Monkey_2 face
movstims{5} = ['HeiniGrunt.avi']; %Monkey_1 face
movstims{6} = ['PrinceGrunt.avi']; %Monkey_3 face

par.useAudDelay = 0; % use the audio delay?
aud_delay{1} = 219; %Chris, audio delay relative to video onset, in ms
aud_delay{2} = 205; %Christoph
aud_delay{3} = 109; %Heini
aud_delay{4} = 129; %Emil
aud_delay{5} = 89; %HeiniGrunt
aud_delay{6} = 183; %PrinceGrunt

%Creating a list of random conditions and what values will be loaded for
%video and/or audio
%tmp = 3 * par.num_audiostims;
tmp = 1 * par.num_audiostims;
par.numTrials = par.numReps * tmp;
stimOrder = [];
for i = 1:par.numReps
    stimOrder = [stimOrder, randperm(tmp)];
end

% save all of the par data into the dat structure
par.condOrder = par.conds(stimOrder);

%par.sndur = cell2mat(aud_dur(stimOrder));
par.stimOrder = stimOrder;
par.audiostims = audiostims;
par.movstims = movstims;
par.aud_delay = aud_delay;
dat.par = par;



%% Reads sound data and checks sampling rate against expected
clear sounddata

disp('Loading sounds');
for n = 1:length(audiostims) % there are always 16 stimuli repeated twice
    [sounddata{n},fs_check] = audioread([sdir '\' audiostims{n}]);
    dat.sound_samplefrequency = fs_check;
end
disp('Sounds loaded');






%% Back to the original code here %%%
%% START EXPERIMENT
% intialisePTB % run subscript to initialise psychtoolbox screen
ListenChar(2);

% check for opengl compatability
AssertOpenGL;

% we do not need perfect timing.
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'ConserveVRAM', 64);

% get screen
screens = Screen('Screens');
% HideCursor;

pause(0.1); % these are here to possibly help with loading all the PTB stuff. I don't know if they work.

% set window
pixdepth = 32;
buffermode = 2; % double buffer
[w, wRect] = Screen('OpenWindow', ARG.screenNumberExperiment, 0, [], pixdepth, buffermode);
[width, height] = Screen('WindowSize', w);
Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

pause(0.1); % these are here to possibly help with loading all the PTB stuff. I don't know if they work.

% clear screen
Screen('FillRect',w, [0 0 0]);
Screen('Flip', w);

% set text parameters
yoffset = 0;
textsize = 32;

drawText('Loading and initialising ...', w, yoffset + 60, [255 255 255], textsize);
Screen('Flip', w);

pause(0.2);

%% loop through trials
for trial = 1:par.numTrials
    
    
    %% Inital trial setup
    % Set the blinkCounter to 0 to give a chance to fixate
    if eyeTrackingOn == 1
        blinkCounter = 0;
        blinkCounterABS = 0;
        timedAbortFlag = 0;
        dat.correct(trial) = 0; % placeholder
        dat.abort(trial) = 0; % placeholder
    end
    
    % flags
    % Initially set all of the event flags to 0
    for xx = 1:size(events,1)
        events{xx,3} = 0;
    end
    % Sounds buffer loading
    par.condOrder(trial) == 1, %A trial
    p = audioplayer(sounddata{par.stimOrder(trial)},dat.sound_samplefrequency);
    aud_delay_tmp = 0;
    %     elseif par.condOrder(trial) == 3, %AV trial
    %         p = audioplayer(sounddata{par.stimOrder(trial)-12},dat.sound_samplefrequency);
    %         aud_delay_tmp = aud_delay{par.stimOrder(trial)-12};
    %     end
    % set up the monitoring figure
    %     drawEyeMonitorFigure
    drawEyeMonitorFigure_local
    
    %% Start trial
    disp(['Started trial: ' int2str(trial), ' ; CondNo:', int2str(par.condOrder(trial)),' ; StimNo: ',int2str(par.stimOrder(trial))]);
    LJ_byte_out(1,ljHandle,LJ_ioPUT_DIGITAL_BIT);
    % Counters
    abortCounter = 1; % Counter for use if trial aborts
    responseCounter = 1; % Allows collecting and storage of all lever responses
    eyeCounter = 1; % Allows collecting and storage of all eye data (gaze & pupil)
    
    %trial start
    trialStart = GetSecs;
    dat.T.start(trial) = GetSecs - trialStart;
    dat.trialNum(trial) = trial;
    
    % Get initial eye position
    dat.eye.dat{trial} = [];
    dat.eye.dat{trial}.timepoint(1) = GetSecs-trialStart; % grab first timepoint
    
    if ARG.MRI == 1
        [Error dat.eye.dat{trial}.gx(eyeCounter)] = ljud_eGet(ljHandle,LJ_ioGET_AIN,0,0,0); % get X eye voltage
        dat.eye.dat{trial}.gx(eyeCounter) = 0 + [((dat.eye.dat{trial}.gx(eyeCounter)-(-5.25))*(1-0))/(5.25-(-5.25))]; % Normalise -5v to 5v voltage between 0 and 1
        [Error dat.eye.dat{trial}.gy(eyeCounter)] = ljud_eGet(ljHandle,LJ_ioGET_AIN,1,0,0); % get Y eye voltage
        dat.eye.dat{trial}.gy(eyeCounter) = 0 + [((dat.eye.dat{trial}.gy(eyeCounter)-(-5.25))*(1-0))/(5.25-(-5.25))]; % Normalise -5v to 5v voltage between 0 and 1
        dat.eye.dat{trial}.fixWindow(eyeCounter) = par.maxVal-centre(1);
        
        newEyePos = scatter(dat.eye.dat{trial}.gx(eyeCounter),dat.eye.dat{trial}.gy(eyeCounter));
        
    elseif ARG.MRI == 0
        [dat.eye.dat{trial}.gx(eyeCounter),dat.eye.dat{trial}.gy(eyeCounter)] = vpx_GetGazePointSmoothed(dat.eye.eye);
        dat.eye.dat{trial}.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter);
        dat.eye.dat{trial}.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter)* -1 + 1;
        
        % Calibrated eye position with manual correction
        dat.eye.dat{trial}.raw.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter);
        dat.eye.dat{trial}.raw.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter);
        dat.eye.dat{trial}.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter) + par.xAdjust;
        dat.eye.dat{trial}.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter) + par.yAdjust;
        
        % other eye parameters
        [dat.eye.dat{trial}.px(eyeCounter),dat.eye.dat{trial}.py(eyeCounter)]=vpx_GetPupilSize(dat.eye.eye);
        dat.eye.dat{trial}.fixWindow(eyeCounter) = par.maxVal-par.centre(1);
        dat.eye.dat{trial}.xAdjust(eyeCounter) = par.xAdjust;
        dat.eye.dat{trial}.yAdjust(eyeCounter) = par.yAdjust;
        
        newEyePos = scatter(dat.eye.dat{trial}.gx(eyeCounter),dat.eye.dat{trial}.gy(eyeCounter));
    end
    
    %% collect eye tracking data while running the trial
    while (GetSecs - trialStart < (events{end,2}/1000))
        
        % Check keys while waiting
        [KeyIsDown, endrt, KeyCode] = KbCheck;
        
        %% Check whether it is time for an event, then run that event
        for e = 1:size(events,1)
            if GetSecs - trialStart > (events{e,2}/1000) && events{e,3} == 0
                
                %%%%%%%%%% Display fixation spot %%%%%%%%%%
                if strcmp(events{e,1},'fixSpotOn')
                    Screen('FillOval', w, [255 255 0], CenterRect([0 0 par.spotSize par.spotSize],wRect) + [0 0 0 0]); %hard circles
                    %Screen('FillOval', w, [255 255 0], CenterRect([0 0 25 25],wRect) + [0 0 0 0]); %hard circles
                    Screen('Flip', w);
                    %%%TIMESTAMP%%%
                    dat.T.fixOn(trial) = GetSecs - trialStart;
                    % set the flag to show the event has occurred.
                    events{e,3} = 1;
                    LJ_byte_out(5+stimOrder(trial),ljHandle,LJ_ioPUT_DIGITAL_BIT)
                    %%%%%%%%%% Remove fixation spot %%%%%%%%%%
                elseif strcmp(events{e,1},'fixSpotOff')
                    Screen('Flip', w);
                    %%%TIMESTAMP%%%
                    dat.T.fixOff(trial) = GetSecs - trialStart;
                    % set the flag to show the event has occurred.
                    events{e,3} = 1;
                    
                    %%%%%%%%%% Start requirement to fixate %%%%%%%%%%
                elseif strcmp(events{e,1},'fixRequiredOn')
                    timedAbortFlag = 1;
                    %%%TIMESTAMP%%%
                    dat.T.fixRequiredOn(trial) = GetSecs - trialStart;
                    % set the flag to show the event has occurred.
                    events{e,3} = 1;
                    
               %%%%%%%%%% Stimulus %%%%%%%%%%
                elseif strcmp(events{e,1},'stimPlay')
                    dat.stimCond(trial) = par.condOrder(trial);
                    dat.stim(trial) = par.stimOrder(trial);
                    if par.condOrder(trial) == 1, % if Audio only trial, just play sound
                        play(p);
                        LJ_byte_out(2,ljHandle,LJ_ioPUT_DIGITAL_BIT);
                        dat.T.stimOn(trial) = GetSecs - trialStart ;
                    end
                    % PsychPortAudio('Start', pahandle,[],[],[]); % last [] to 1 for PTB timing info
                    Error = ljud_ePut(ljHandle, LJ_ioPUT_DIGITAL_PORT, 14, 5, 0);
                    %%%TIMESTAMP%%%
                    dat.stimCond(trial) = par.condOrder(trial);
                    dat.stim(trial) = par.stimOrder(trial);
                    %dat.T.stimOn(trial) = GetSecs - trialStart ;
                    
                    events{e,3} = 1;     
                    %%%%%%%%%% remove requirement to fixate %%%%%%%%%%
                elseif strcmp(events{e,1},'fixRequiredOff')
                    stop(p);
                    timedAbortFlag = 0;
                    %%%TIMESTAMP%%%
                    dat.T.fixRequiredOff(trial) = GetSecs - trialStart;
                    % set the flag to show the event has occurred.
                    events{e,3} = 1;
                    
                    
                    %%%%%%%%%% Reward start %%%%%%%%%%
                elseif strcmp(events{e,1},'rewardOn') % && respType == 1
                    %                     if ARG.MRI == 0
                    Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,5,0);
                    LJ_byte_out(4,ljHandle,LJ_ioPUT_DIGITAL_BIT)
                    %                     elseif ARG.MRI == 1
                    %                         Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0);
                    %                     end
                    %%%TIMESTAMP%%%
                    dat.T.rewardOn(trial) = GetSecs - trialStart ;
                    dat.T.reward(trial) = 1;
                    events{e,3} = 1;
                    
                    %%%%%%%%%% Reward stop %%%%%%%%%%
                elseif strcmp(events{e,1},'rewardOff') % && respType == 1
                    if ARG.MRI == 0
                        Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0);
                    elseif ARG.MRI == 1
                        Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0);
                    end
                    LJ_byte_out(5,ljHandle,LJ_ioPUT_DIGITAL_BIT)
                    %%%TIMESTAMP%%%
                    dat.T.rewardOff(trial) = GetSecs - trialStart ;
                    events{e,3} = 1;
                    
                    %%%%%%%%%% end the trial %%%%%%%%%%
                elseif strcmp(events{e,1},'trialEnd')
                    dat.T.trialEnd(trial) = GetSecs - trialStart;
                    events{e,3} = 1;
                    Screen('FillRect', w, [0 0 0]);
                    Screen('Flip', w);
                    ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0);
                    
                    %                     pause(0.1)
                    %                     LJ_byte_out(50+stimOrder(trial),ljHandle,LJ_ioPUT_DIGITAL_BIT)
                    
                    pause(0.1)
                    
                    LJ_byte_out(7,ljHandle,LJ_ioPUT_DIGITAL_BIT)
                end
            end
        end
        
        %% EYE TRACKING
        if eyeTrackingOn == 1
            
            % Save eye position at prespecifed sampling rate
            if (GetSecs - trialStart > dat.eye.dat{trial}.timepoint(eyeCounter) + dat.eye.sampling) && eyeTrackingOn == 1
                eyeCounter = eyeCounter + 1;
                
                if ARG.MRI == 1
                    [Error dat.eye.dat{trial}.gx(eyeCounter)] = ljud_eGet(ljHandle,LJ_ioGET_AIN,0,0,0); % get X eye voltage
                    dat.eye.dat{trial}.gx(eyeCounter) = 0 + [((dat.eye.dat{trial}.gx(eyeCounter)-(-5))*(1-0))/(5-(-5))]; % Normalise -5v to 5v voltage between 0 and 1
                    [Error dat.eye.dat{trial}.gy(eyeCounter)] = ljud_eGet(ljHandle,LJ_ioGET_AIN,1,0,0); %get Y eye voltage
                    dat.eye.dat{trial}.gy(eyeCounter) = 0 + [((dat.eye.dat{trial}.gy(eyeCounter)-(-5))*(1-0))/(5-(-5))]; % Normalise -5v to 5v voltage between 0 and 1
                    
                elseif ARG.MRI == 0
                    [dat.eye.dat{trial}.gx(eyeCounter),dat.eye.dat{trial}.gy(eyeCounter)]=vpx_GetGazePointSmoothed(dat.eye.eye);
                    
                    % Correct the flipped eye positions
                    dat.eye.dat{trial}.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter); % the eye tracker was mirrored, I think only in the y directions. This flips it.
                    dat.eye.dat{trial}.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter)* -1 + 1;
                    
                    % Calibrated eye position with manual correction
                    dat.eye.dat{trial}.raw.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter);
                    dat.eye.dat{trial}.raw.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter);
                    dat.eye.dat{trial}.gx(eyeCounter) = dat.eye.dat{trial}.gx(eyeCounter) + par.xAdjust;
                    dat.eye.dat{trial}.gy(eyeCounter) = dat.eye.dat{trial}.gy(eyeCounter) + par.yAdjust;
                    
                    % other eye parameters
                    [dat.eye.dat{trial}.px(eyeCounter),dat.eye.dat{trial}.py(eyeCounter)]=vpx_GetPupilSize(dat.eye.eye);
                    dat.eye.dat{trial}.fixWindow(eyeCounter) = par.maxVal-par.centre(1);
                    dat.eye.dat{trial}.xAdjust(eyeCounter) = par.xAdjust;
                    dat.eye.dat{trial}.yAdjust(eyeCounter) = par.yAdjust;
                    
                end
                dat.eye.dat{trial}.fixWindow(eyeCounter) = par.maxVal - par.centre(1); % Fix window size
                dat.eye.dat{trial}.timepoint(eyeCounter) = GetSecs - trialStart;
                
                % change window size/position
                if ( KeyIsDown==1 && KeyCode(187)==1 ) % 187 is keycode for + key, ASCII is 61
                    par.maxVal = par.maxVal + 0.001;
                    par.minVal = par.minVal - 0.001;
                end
                if ( KeyIsDown==1 && KeyCode(189)==1 ) % 187 is keycode for - key, ASCII is 45
                    par.maxVal = par.maxVal - 0.001;
                    par.minVal = par.minVal + 0.001;
                end
                
                %change Centre position (crude, but useful for now
                if ( KeyIsDown==1 && KeyCode(37)==1 ) % 37 is keycode for left arrow
                    par.xAdjust = par.xAdjust - 0.001;
                end
                if ( KeyIsDown==1 && KeyCode(38)==1 ) % 38 is keycode for up arrow
                    par.yAdjust = par.yAdjust + 0.001;
                end
                if ( KeyIsDown==1 && KeyCode(39)==1 ) % 39 is keycode for right arrow
                    par.xAdjust = par.xAdjust + 0.001;
                end
                if ( KeyIsDown==1 && KeyCode(40)==1 ) % 40 is keycode for down arrow
                    par.yAdjust = par.yAdjust - 0.001;
                end
                
                % Counter for eye position increases outside of specified window
                if par.eyeTrackAbort == 1 && timedAbortFlag == 1
                    if dat.eye.dat{trial}.gx(eyeCounter) > par.maxVal || dat.eye.dat{trial}.gy(eyeCounter) > par.maxVal || dat.eye.dat{trial}.gx(eyeCounter) < par.minVal || dat.eye.dat{trial}.gy(eyeCounter) < par.minVal
                        blinkCounter = blinkCounter + 1;
                        blinkCounterABS = blinkCounter + 1;
                    else
                        blinkCounter = 0;
                    end
                end
                
                % If poor fixation, abort, send MRI trig so it matches scanner trial and write some values
                if blinkCounter >= par.blinkAllowance && timedAbortFlag == 1 %&& events(stimOnPosIdx,3) == 1
                    Error = ljud_ePut(ljHandle,LJ_ioPUT_DAC,0,0,0); % juice OFF
                    stop(p); %sound off
                    Screen('FillRect', w, [0 0 0]); %Screen turns black when making error not red
                    Screen('Flip', w);
                    dat.abort(trial) = 1;
                    dat.correct(trial) = 0;
                    dat.T.abort_time(trial) = GetSecs - trialStart;
                    disp(['Aborted trial: ' int2str(trial)]);
                    pause(par.abortTime/1000); % give a timeout
                    Screen('FillRect', w, [0 0 0]);
                    Screen('Flip', w);
                    disp('-------------------------');
                    break
                    
                else
                    dat.eye.blinkCounter(trial) = blinkCounter;
                    dat.eye.blinkCounterABS(trial) = blinkCounterABS;
                    dat.abort(trial) = 0;
                end
            end
            
            %% plot the eye position every 5 samples.
            plotEvery = 5;
            if rem(eyeCounter,plotEvery) == 0
                % update eye position figure
                axes(ax2);
                oldEyePos = newEyePos;
                newEyePos = scatter(dat.eye.dat{trial}.gx(eyeCounter),dat.eye.dat{trial}.gy(eyeCounter),'k','filled');
                if ( KeyIsDown==1 && KeyCode(187)==1 ) || ( KeyIsDown==1 && KeyCode(189)==1 )
                    delete(lh1);delete(lh2);delete(lh3);delete(lh4);
                    lh1=line([par.minVal+par.xOffset,par.minVal+par.xOffset],[par.minVal+par.yOffset,par.maxVal+par.yOffset],'Color','b','LineWidth',0.5);
                    lh2=line([par.maxVal+par.xOffset,par.maxVal+par.xOffset],[par.minVal+par.yOffset,par.maxVal+par.yOffset],'Color','b','LineWidth',0.5);
                    lh3=line([par.minVal+par.xOffset,par.maxVal+par.xOffset],[par.minVal+par.yOffset,par.minVal+par.yOffset],'Color','b','LineWidth',0.5);
                    lh4=line([par.minVal+par.xOffset,par.maxVal+par.xOffset],[par.maxVal+par.yOffset,par.maxVal+par.yOffset],'Color','b','LineWidth',0.5);
                end
                
                delete(oldEyePos)
                viscircles([par.centre(1),par.centre(2)],0.001);
                viscircles([par.centre(1) + par.xOffset,par.centre(2) + par.yOffset],0.001);
                drawnow;
                
                % update trial figure (SLOW - maybe avoid axis switching in loop)
                axes(ax1)
                delete(lineTime);
                lineTime = line([dat.eye.dat{trial}.timepoint(eyeCounter)*1000,dat.eye.dat{trial}.timepoint(eyeCounter)*1000],[0,0.5],'Color','r','LineWidth',1);
                drawnow;
            end
        end
        
        %% record correct trials.
        if dat.abort(trial) == 0
            dat.correct(trial) = 1;
        else
            dat.correct(trial) = 0;
        end
        
        %% Check keys
        if ( KeyIsDown==1 && KeyCode(27)==1 ) %was "esc" pressed?
            Screen('FillRect', w, [0 0 0]);
            Screen('Flip', w);
            pause(0.500);
            Screen('CloseAll');
            ShowCursor;
            ListenChar;
            Priority(0);
            disp('Experimenter exit ...');
            vpx_Unload;
            error('exit');
        end
    end
    
    %% save the data
    tic
    %save(['C:\TESTDATA\',par.Expt,'_',datafile],'-v6','dat');
    save([ddir,par.Expt,'_',datafile],'dat');
    if (trial) == 0
        disp('-------------------------');
    end
    
end % trial loop

%% save the data
save(['C:\TESTDATA\',par.Expt,'_',datafile],'-v6','dat');
save([ddir,par.Expt,'_',datafile],'dat');

%% cleanup at end of experiment
PsychPortAudio('Close');
disp('Normal exit')
Screen('CloseAll');
ShowCursor;
fclose('all');
ListenChar;
Priority(0);
vpx_Unload;
close(12);

%% Print final results to screen, to store in lab book spreadsheet
clc

summary = cell(2,4);
summary{1,1} = 'Correct';
summary{1,2} = 'Attempt';
summary{1,3} = 'Percent Correct';
summary{1,4} = 'Data';

summary{2,1} = length(dat.trialNum) - sum(dat.abort);
summary{2,2} = length(dat.trialNum);
summary{2,3} = summary{2,1}/summary{2,2}*100;
summary{2,4} = [par.Expt,'_',datafile];

fprintf(['Variables to save in lab book spreadsheet: \n \n']);
fprintf([summary{1,1},' = ',num2str(summary{2,1}),'\n']);
fprintf([summary{1,2},' = ',num2str(summary{2,2}),'\n']);
fprintf([summary{1,3},' = ',num2str(summary{2,3}),'\n']);
fprintf([summary{1,4},' = ',summary{2,4},'\n','\n']);