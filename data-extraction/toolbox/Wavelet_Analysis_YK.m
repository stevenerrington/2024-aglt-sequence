%Wavelet_Analysis_YK.m


%From newtimef.m
% Outputs:
%            ersp   = (nfreqs,timesout) matrix of log spectral diffs from baseline
%                     (in dB log scale or absolute scale). Use the 'plot' output format
%                     above to output the ERSP as shown on the plot.
%            itc    = (nfreqs,timesout) matrix of complex inter-trial coherencies.
%                     itc is complex -- ITC magnitude is abs(itc); ITC phase in radians
%                     is angle(itc), or in deg phase(itc)*180/pi.
%          powbase  = baseline power spectrum. Note that even, when selecting the
%                     the 'trialbase' option, the average power spectrum is
%                     returned (not trial based). To obtain the baseline of
%                     each trial, recompute it manually using the tfdata
%                     output described below.
%            times  = vector of output times (spectral time window centers) (in ms).
%            freqs  = vector of frequency bin centers (in Hz).
%         erspboot  = (nfreqs,2) matrix of [lower upper] ERSP significance.
%          itcboot  = (nfreqs) matrix of [upper] abs(itc) threshold.
%           tfdata  = optional (nfreqs,timesout,trials) time/frequency decomposition
%                      of the single data trials. Values are complex.

clear all
clc;

task_flag = 18;
%1:tonotopy(tonos); 2:coo; 3 tonocoo; 4:hrm;   5:snd40;
%7: syntax_single;   8: syntax_multi;
% 9. FRA;     10.FRAlow
%11: pitchNRN    12: pitchNRN512
%13:pitchNRNitr  14:pitchRNR512
%15:WM
%16: pitchNRN all stimuli; 17: pitch NRN512 all stimuli
%18 SoundON (syntax task)

Linuxformat =0; %0: Windows; 1: Linux

freq_range=[4 100];%[1.2 100]


% Put the contvars event name to align [t1(1);t2(1) ...] format
if  task_flag == 1 %tonotopy task
    event2analyse = {'AllsndonCorrect';'snd00On_Correct';...
        'snd01On_Correct';'snd02On_Correct';'snd03On_Correct';'snd04On_Correct';'snd05On_Correct';...
        'snd06On_Correct';'snd07On_Correct';'snd08On_Correct';'snd09On_Correct';'snd10On_Correct';...
        'snd11On_Correct';'snd12On_Correct';'snd13On_Correct';'snd14On_Correct';'snd15On_Correct';...
        'snd16On_Correct';'snd17On_Correct';'snd18On_Correct';'snd19On_Correct';'snd20On_Correct';...
        'snd21On_Correct';'snd22On_Correct';'snd23On_Correct';'snd24On_Correct';'snd25On_Correct';...
        'snd26On_Correct';'snd27On_Correct';'snd28On_Correct';'snd29On_Correct';'snd30On_Correct';...
        'snd31On_Correct';'snd32On_Correct';'snd33On_Correct';'snd34On_Correct';'snd35On_Correct';...
        'snd36On_Correct';'snd37On_Correct';'snd38On_Correct';'snd39On_Correct';'snd40On_Correct';...
        'snd41On_Correct';'snd42On_Correct';'snd43On_Correct';'snd44On_Correct';'snd45On_Correct';...
        'snd46On_Correct';'snd47On_Correct';'snd48On_Correct';'snd49On_Correct';'snd50On_Correct';...
        'snd51On_Correct';'snd52On_Correct';'snd53On_Correct';'snd54On_Correct';'snd55On_Correct';...
        'snd56On_Correct';'snd57On_Correct';'snd58On_Correct';'snd59On_Correct';'snd60On_Correct';...
        'snd61On_Correct';'snd62On_Correct';'snd63On_Correct';'snd64On_Correct';'snd65On_Correct';...
        'snd66On_Correct';'snd67On_Correct';'snd68On_Correct';'snd69On_Correct';'snd70On_Correct';...
        'snd71On_Correct';'snd72On_Correct';'snd73On_Correct';'snd74On_Correct';'snd75On_Correct';...
        'snd76On_Correct';'snd77On_Correct';'snd78On_Correct';'snd79On_Correct';'snd80On_Correct';...
        'snd81On_Correct';'snd82On_Correct';'snd83On_Correct';'snd84On_Correct'; 'snd00_Correct_BARUP'
        };
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -0.5; %start trial relative to sound onset (s)
    epochmax =1; %end trial relative to sound onset (s)
    plot_marktimes = [0];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag == 2 %coo task
    event2analyse = {'AllsndonCorrect';...
        'snd01On_Correct';'snd02On_Correct';'snd03On_Correct';'snd04On_Correct';'snd05On_Correct';...
        'snd06On_Correct';'snd07On_Correct';'snd08On_Correct';'snd09On_Correct';'snd10On_Correct';...
        'snd11On_Correct';'snd12On_Correct';'snd13On_Correct';'snd14On_Correct';'snd15On_Correct';...
        'snd16On_Correct';'snd17On_Correct';'snd18On_Correct';'snd19On_Correct';'snd20On_Correct';...
        'snd21On_Correct';'snd22On_Correct';'snd23On_Correct';'snd24On_Correct';'snd25On_Correct';...
        'snd26On_Correct';'snd27On_Correct';'snd28On_Correct';'snd29On_Correct';'snd30On_Correct';...
        'snd31On_Correct';'snd32On_Correct';'snd33On_Correct';'snd34On_Correct';'snd35On_Correct';...
        'snd36On_Correct';'snd37On_Correct';'snd38On_Correct';'snd39On_Correct';'snd40On_Correct';...
        'snd41On_Correct';'snd42On_Correct';'snd43On_Correct';'snd44On_Correct';'snd45On_Correct';...
        'snd46On_Correct';'snd47On_Correct';'snd48On_Correct';'snd49On_Correct';'snd50On_Correct';...
        'snd51On_Correct';'snd52On_Correct';'snd53On_Correct';'snd54On_Correct';'snd55On_Correct';...
        'snd56On_Correct'; 'snd00_Correct_BARUP'};
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -0.5; %start trial relative to sound onset (s)
    epochmax =1; %end trial relative to sound onset (s)
    plot_marktimes = [0];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag == 3 %tonocoo task
    event2analyse = {'AllsndonCorrect'; 'snd00On_Correct';...
        'snd01On_Correct';'snd02On_Correct';'snd03On_Correct';'snd04On_Correct';'snd05On_Correct';...
        'snd06On_Correct';'snd07On_Correct';'snd08On_Correct';'snd09On_Correct';'snd10On_Correct';...
        'snd11On_Correct';'snd12On_Correct';'snd13On_Correct';'snd14On_Correct';'snd15On_Correct';...
        'snd16On_Correct';'snd17On_Correct';'snd18On_Correct';'snd19On_Correct';'snd20On_Correct';...
        'PT_On_Correct'; 'Coo_On_Correct'; ; 'snd00_Correct_BARUP'};
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -0.5; %start trial relative to sound onset (s)
    epochmax =1; %end trial relative to sound onset (s)
    plot_marktimes = [0];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag == 4 %hrm task
    event2analyse = {'AllsndonCorrect';...
        'snd00On_Correct';'snd01On_Correct';'snd02On_Correct';'snd03On_Correct';'snd04On_Correct';'snd05On_Correct';...
        'snd06On_Correct';'snd07On_Correct';'snd08On_Correct';'snd09On_Correct';'snd10On_Correct';...
        'snd11On_Correct';'snd12On_Correct';'Coo_Correct_On';...
        'PT_On_Correct'; 'BPN_On_Correct'; 'BPcoo_On_Correct';'BScoo_On_Correct';...
        'F0_On_Correct'; 'F1_On_Correct'; 'F0F1_On_Correct'};
    
    
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -0.5; %start trial relative to sound onset (s)
    epochmax =1; %end trial relative to sound onset (s)
    plot_marktimes = [0];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
    
elseif task_flag == 5 %snd40
    event2analyse = {'AllsndonCorrect';...
        'snd01On_Correct';'snd02On_Correct';'snd03On_Correct';'snd04On_Correct';'snd05On_Correct';...
        'snd06On_Correct';'snd07On_Correct';'snd08On_Correct';'snd09On_Correct';'snd10On_Correct';...
        'snd11On_Correct';'snd12On_Correct';'snd13On_Correct';'snd14On_Correct';'snd15On_Correct';...
        'snd16On_Correct';'snd17On_Correct';'snd18On_Correct';'snd19On_Correct';'snd20On_Correct';...
        'snd21On_Correct';'snd22On_Correct';'snd23On_Correct';'snd24On_Correct';'snd25On_Correct';...
        'snd26On_Correct';'snd27On_Correct';'snd28On_Correct';'snd29On_Correct';'snd30On_Correct';...
        'snd31On_Correct';'snd32On_Correct';'snd33On_Correct';'snd34On_Correct';'snd35On_Correct';...
        'snd36On_Correct';'snd37On_Correct';'snd38On_Correct';'snd39On_Correct';'snd40On_Correct';...
        'MC_On_Correct'; 'OA_On_Correct'};
    
    basemin = -1000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =2; %end trial relative to sound onset (s)
    plot_marktimes = [0];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag ==7 %syntax_single
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    freq_range=[2.5 100];%[1.2 100]
    event2analyse ={'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';'snd05_on_correct';...
        'snd06_on_correct';'snd07_on_correct';'snd08_on_correct'};
    %Just for test
    %event2analyse ={'snd01_on_correct'};
    
    %The figure only shows 500ms
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    %Original setting
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =4; %end trial relative to sound onset (s)
    %epochmin = -4; %start trial relative to sound onset (s)
    %epochmax = 6; %end trial relative to sound onset (s)
    plot_marktimes = [0 413 563 976 1126 1539 1689 2102 2253 2666];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
    
elseif task_flag == 8 %syntax_multi
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    freq_range=[2.5 100];%[1.2 100]
    event2analyse ={'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';'snd05_on_correct';...
        'snd06_on_correct';'snd07_on_correct';'snd08_on_correct'};
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =4; %end trial relative to sound onset (s)
    plot_marktimes = [0 413 563 976 1126 1539 1689 2102 2253 2666];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag ==9 %FRA
    event2analyse = {'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';...
        'snd05_on_correct';'snd06_on_correct';'snd07_on_correct';'snd08_on_correct';...
        'snd09_on_correct';'snd10_on_correct';'snd11_on_correct';'snd12_on_correct';...
        'snd13_on_correct';'snd14_on_correct';'snd15_on_correct';'snd16_on_correct';...
        'snd17_on_correct';'snd18_on_correct';'snd19_on_correct';'snd20_on_correct';...
        'snd21_on_correct';'snd22_on_correct';'snd23_on_correct';'snd24_on_correct';...
        'snd25_on_correct';'snd26_on_correct';'snd27_on_correct';'snd28_on_correct';...
        'snd29_on_correct';'snd30_on_correct';'snd31_on_correct';'snd32_on_correct';...
        'snd33_on_correct';'snd34_on_correct';'snd35_on_correct'};
    basemin = -1000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =2; %0.1 s shorter than poststim %end trial relative to sound onset (s)
    plot_marktimes = [0 300];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
    
elseif task_flag ==10 %FRAlow
    event2analyse = {'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';...
        'snd05_on_correct';'snd06_on_correct';'snd07_on_correct';'snd08_on_correct';...
        'snd09_on_correct';'snd10_on_correct'};
    basemin = -1000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =2; %0.1 s shorter than poststim %end trial relative to sound onset (s)
    plot_marktimes = [0 300];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag ==11   %pitchNRN
    
    freq_range=[2.5 150];
    %freq_range=[2.5 100]; %Based on the cycle [3 0.7] I currently set, the lowest fq I can extract is 4 Hz.
    event2analyse = {'HC8_correct_on';'HC16_correct_on';'HC32_correct_on';'HC64_correct_on';...
        'HC128_correct_on'; 'HC256_correct_on';'RIN8_correct_on';'RIN16_correct_on';'RIN32_correct_on';...
        'RIN64_correct_on'; 'RIN128_correct_on'; 'RIN256_correct_on'};    %add 'HC512' and 'RIN512_correct_on' when necesarry
    
    %event2analyse = {'RIN16_correct_on'};    %add 'HC512' and 'RIN512_correct_on' when necesarry
    
    % basemin = 150; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    % basemax=450;
    
    basemin = 0; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=500;
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap significance prob. level.
    
    
elseif task_flag ==12    %pitchNRN512
    
    freq_range=[2.5 150];
    %freq_range=[2.5 100];
    event2analyse = {'HC8_correct_on';'HC16_correct_on';'HC32_correct_on';'HC64_correct_on';...
        'HC128_correct_on'; 'HC256_correct_on';'HC512_correct_on';'RIN8_correct_on';'RIN16_correct_on';'RIN32_correct_on';...
        'RIN64_correct_on'; 'RIN128_correct_on'; 'RIN256_correct_on'; 'RIN512_correct_on';};    %add 'HC512' and 'RIN512_correct_on' when necesarry
    
    
    
    % basemin = -3000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    %epochmin = -3.0; %start trial relative to sound onset (s)
    %epochmax =3; %end trial relative to sound onset (s)
    
    %Previous settings
%     basemin = 150; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
%     basemax=450;
    
    basemin = 0;%0 %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=500;%500
    
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
    
    
elseif task_flag ==13    %pitchNRNitr
    alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
    freq_range=[2.5 200];
    event2analyse = {'RINitr1';'RINitr4';'RINitr8';'RINitr16';'RINitr32';'RINitr64';'HC'};
    % basemin = -3000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    %epochmin = -3.0; %start trial relative to sound onset (s)
    %epochmax =3; %end trial relative to sound onset (s)
    
    %Previous settings
    basemin = 150; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=450;
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
    eventnames = {'RINitr1';'RINitr4';'RINitr8';'RINitr16';'RINitr32';'RINitr64';'HC'};
    
elseif task_flag ==14    %pitchRNR512
    freq_range=[2.5 200];
    event2analyse = {'HC8_correct_on';'HC16_correct_on';'HC32_correct_on';'HC64_correct_on';...
        'HC128_correct_on'; 'HC256_correct_on';'HC512_correct_on';'RIN8_correct_on';'RIN16_correct_on';'RIN32_correct_on';...
        'RIN64_correct_on'; 'RIN128_correct_on'; 'RIN256_correct_on'; 'RIN512_correct_on';};    %add 'HC512' and 'RIN512_correct_on' when necesarry
    % basemin = -3000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    %epochmin = -3.0; %start trial relative to sound onset (s)
    %epochmax =3; %end trial relative to sound onset (s)
    
    %Previous settings
    basemin = 650; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=950;
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag ==15    %WM
    event2analyse = {'Go_correct_on';'Go_incorrect_on';'NoGo_correct_on';'NoGo_incorrect_on'};
    %Previous settings
    preevetime = 0.5; %time before the event to align
    postevetime =5; %time after the event to align
    
    
    % basemin = -3000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    %epochmin = -3.0; %start trial relative to sound onset (s)
    %epochmax =3; %end trial relative to sound onset (s)
    
    %Previous settings
    basemin = -1000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    epochmin = -0.5; %start trial relative to sound onset (s)
    epochmax =5; %end trial relative to sound onset (s)
    plot_marktimes = [0 100 1100 1200];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
elseif task_flag==16 %pitchNRNallstimuli
    %pitchNRN
    event2analyse = {'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';'snd05_on_correct';...
        'snd06_on_correct';'snd07_on_correct';'snd08_on_correct';'snd09_on_correct';'snd10_on_correct';...
        'snd11_on_correct';'snd12_on_correct';'snd13_on_correct';'snd14_on_correct';'snd15_on_correct';...
        'snd16_on_correct';'snd17_on_correct';'snd18_on_correct';'snd19_on_correct';'snd20_on_correct';...
        'snd21_on_correct';'snd22_on_correct';'snd23_on_correct';'snd24_on_correct';'snd25_on_correct';...
        'snd26_on_correct';'snd27_on_correct';'snd28_on_correct';'snd29_on_correct';'snd30_on_correct';...
        'snd31_on_correct';'snd32_on_correct';'snd33_on_correct';'snd34_on_correct';'snd35_on_correct';...
        'snd36_on_correct'};
    basemin = 150; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=450;
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap significance prob. level.
    
elseif task_flag==17 %pitchNRN512allstimuli
    event2analyse = {'snd01_on_correct';'snd02_on_correct';'snd03_on_correct';'snd04_on_correct';'snd05_on_correct';...
        'snd06_on_correct';'snd07_on_correct';'snd08_on_correct';'snd09_on_correct';'snd10_on_correct';...
        'snd11_on_correct';'snd12_on_correct';'snd13_on_correct';'snd14_on_correct';'snd15_on_correct';...
        'snd16_on_correct';'snd17_on_correct';'snd18_on_correct';'snd19_on_correct';'snd20_on_correct';...
        'snd21_on_correct';'snd22_on_correct';'snd23_on_correct';'snd24_on_correct';'snd25_on_correct';...
        'snd26_on_correct';'snd27_on_correct';'snd28_on_correct';'snd29_on_correct';'snd30_on_correct';...
        'snd31_on_correct';'snd32_on_correct';'snd33_on_correct';'snd34_on_correct';'snd35_on_correct';...
        'snd36_on_correct';'snd37_on_correct';'snd38_on_correct';'snd39_on_correct';'snd40_on_correct';...
        'snd41_on_correct';'snd42_on_correct'};
    % basemin = -3000; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    %epochmin = -3.0; %start trial relative to sound onset (s)
    %epochmax =3; %end trial relative to sound onset (s)
    
    %Previous settings
    basemin = 150; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=450;
    epochmin = -0.9; %start trial relative to sound onset (s)
    epochmax =2.75; %end trial relative to sound onset (s)
    plot_marktimes = [0 500 1400 2300];
    alpha_val = 0.05; %compute two-tailed bootstrap s, ignificance prob. level.

elseif task_flag ==18 %Sound On for syntax tasks
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    freq_range=[2.5 100];%[1.2 100]
    event2analyse ={'snd99'};
    %Just for test
    %event2analyse ={'snd01_on_correct'};
    
    %The figure only shows 500ms
    basemin = -500; %duration prior to sound onset (ms): needs to follow the matrix size (basemin = epochmin)
    basemax=0;
    %Original setting
    epochmin = -1; %start trial relative to sound onset (s)
    epochmax =4; %end trial relative to sound onset (s)
    %epochmin = -4; %start trial relative to sound onset (s)
    %epochmax = 6; %end trial relative to sound onset (s)
    plot_marktimes = [0 413 563 976 1126 1539 1689 2102 2253 2666];
    alpha_val = 0.01; %compute two-tailed bootstrap s, ignificance prob. level.
    
end

maxfreq = max(freq_range);

%% INPUT---------------------------------------------
% PathName=uigetdir('Choose a directory');
% cd (PathName);
dirnames =uipickfiles( 'Prompt','Pick directorie(s) that includes Level mat files');
ndir = size(dirnames,2);

for jj = 1:ndir
    %Prep for a directory for results
    loaddir=char(dirnames(jj));
    cd(loaddir);
    
    if Linuxformat == 1
        savdir =[loaddir '/alpha' num2str(alpha_val) 'maxfreq' num2str(maxfreq) '_baseline' num2str(basemin) '-' num2str(basemax) '_noiseremoved'];
    elseif Linuxformat == 0
        savdir =[loaddir '\alpha' num2str(alpha_val) 'maxfreq' num2str(maxfreq) '_baseline' num2str(basemin)  '-' num2str(basemax)];
    end
    mkdir(savdir);
    
    %Parameter settings
    num_lev =length(event2analyse); %This needs to be matched the event extracted as Level
    num_channels = 1;
    srate = 1000;
    padratio = 2; %Original 2
    
    Elecs=1;
    maxersp = 6; %db in power
    outtimes=1500; %200 for data points in time for cohrence and power analysis %Original 1500
    
    itcplot = 'on'; %intertrial coherence plot
    
    %     for level=1:num_lev
    %         mkdir([savdir '\Level_ref2_' num2str(level)]);
    %      end
    
    
    for level = 1:num_lev
        
        load(['Level_' num2str(level) '_' char(event2analyse(level)) '.mat'])
        display([loaddir 'Level_' num2str(level) '_' char(event2analyse(level)) '.mat']);
        %load([loaddir '\Level_' num2str(level)]);
        if strcmp('eventnames', event2analyse(level))
            error('Error. eventnames does not match to the event2analyse');
        end
        
        data =X(Elecs,1:(abs(epochmin)+epochmax)*1000,:);
        %         data=diff(data,1,2);
        elecs = size(data,1);
        elec=1;
        b=size(data,3); %b = number of trials
        in(:,1)=[1:b]';
        in(:,2)=abs(epochmin)*ones(b,1); 
        
        eeglab
        EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'Level', 'srate',srate, 'pnts',0, 'xmin',0, 'nbchan',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_importepoch( EEG, 'in', { 'Epoch', 'stim'}, 'latencyfields',{ 'stim'}, 'timeunit',1, 'headerlines',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'stim'  }, [epochmin         epochmax], 'newname', 'Level epochs', 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, [basemin    0]);
        EEG = eeg_checkset( EEG );
        clear  X_C delay_data_comb delay_data_ind data in temp
        close all
        %set(gcf,'Position',[1 1 1600 1000]);
        
        
        %Figure settings-------------------------------------
        title([loaddir  event2analyse(level)], 'Fontsize', 6);
        fprintf('\n\n');
        fprintf('Level number (of 4) %d\n', level)
        fprintf('\n\n');
        fprintf('Electrode number (of 16) %d\n', elec)
        

        
        switch task_flag
            case {7,8} %For syntax MS
                %Usual trial-average analysis
                [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
                    1, elec, [EEG.xmin EEG.xmax]*srate,[2 0.5], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
                    'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 200, 'baseboot',1,'rmerp','off', ...
                    'erspmax', maxersp, 'plotersp','on', 'plotitc',itcplot,'baseline',[basemin basemax],'marktimes',plot_marktimes);
%             case{11,12,13,14,16,17} %For pitch MS (Note wavelet cycles are different)
%                 [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
%                     1, elec, [EEG.xmin EEG.xmax]*srate,[3 0.5], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
%                     'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 300, 'baseboot',1,'rmerp','off', ...
%                     'erspmax', maxersp, 'plotersp','on', 'plotitc',itcplot,'baseline',[basemin basemax],'marktimes',plot_marktimes);
                case{9,10,11,12,13,14,16,17,18} %For pitch MS (Note wavelet cycles are different)
                [ersp,itc,powbase,times,freqs,erspboot,itcboot,alltfX] = pop_newtimef(EEG, ...
                    1, elec, [EEG.xmin EEG.xmax]*srate,[3 0.7], 'maxfreq',maxfreq, 'freqs',freq_range,'padratio', padratio, ...
                    'plotphase', 'off', 'timesout', outtimes, 'alpha', alpha_val, 'naccu', 200, 'baseboot',1,'rmerp','off', ...
                    'erspmax', maxersp, 'plotersp','on', 'plotitc',itcplot,'baseline',[basemin basemax],'marktimes',plot_marktimes);
        end
        
        
        
        
        
        if Linuxformat == 1
            save([savdir  '/Level_' num2str(level) '_' char(event2analyse(level))],'ersp','itc','powbase','freqs','times','erspboot','itcboot','alltfX',...
                'X', 'event2analyse','preevetime','postevetime','xPos', 'yPos');
        elseif Linuxformat == 0
            save([savdir  '\Level_' num2str(level) '_' char(event2analyse(level))],'ersp','itc','powbase','freqs','times','erspboot','itcboot','alltfX',...
                'X', 'event2analyse','preevetime','postevetime','xPos', 'yPos');
        end
        
        
        clear alltfX times freqs erspboot itcboot ersp itc powbase X preevetime postevetime xPos yPos
        
        if Linuxformat == 1
            saveas(gcf, [ savdir '/Level' num2str(level)], 'fig')
        elseif Linuxformat == 0
            saveas(gcf, [ savdir '\Level' num2str(level)], 'fig')
        end
    end
end

