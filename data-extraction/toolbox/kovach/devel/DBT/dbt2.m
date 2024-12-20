classdef dbt2

% The band limited data class for efficient representation of
% band-limited data using the Nyquist sampling theorem. 
% 
% Use:
%
%  B = dbt(X,Fs,BW)
%
%  X  - Signal as a column vector
%  Fs - Sampling frequency
%  BW - Bandwidth of decomposition
% 
%
%  B = dbt(X,Fs,Bw, ['option'], [value])
%  Options:
%
%  	offset   -  offset of the first band from 0hz (default = 0)
%       padding  - 'time': pad signal in the time domain, changing duration (default)
%                  'frequency': pad in the frequency domain, changing sampling rate
%       shoulder - (0 - 1) degree of overlap between neighboring bands (default = 0)           
%
%  The overlapping portions of the bands are windowed with a  taper:
%	     ____   __________   ________   _________
%	         \ /	      \ /        \ /   
% .  .  .	  X            X          X           . . . 
% 	         / \	      / \        / \
%            |-----------|           |-|
%    	               BW              shoulder
%      
% The taper is defined so that squares sum to 1


%     C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/DBT/dbt2.m $
% $Revision: 440 $
% $Date: 2013-12-13 19:45:36 -0600 (Fri, 13 Dec 2013) $
% $Author: ckovach $
% ------------------------------------------------

    properties 
        
        blrep = [];  %%% Band-limited analytic representaiton of the signal (time x frequency band)
        bandwidth = [];  %%% bandwidth parameter
        sampling_rate = []; %%% Sampling rate after downsampling
        offset = 0;	%%% Offset if the first band doesn't begin at 0 (a highpass cutoff)
        time = [];     %%% Time points for the rows in blrep
        frequency = [];  %%% Center frequencies for the columns in blrep
        bands =[];     %%% Bands for the columns of blrep 
        fullN = 0;    %%% Length of the signal after padding
        fullFS = 0;   %%% Sampling frequency of the reconstructed signal
        Norig = 0;    %%% Original signal length before padding
        nyqval = 0;   %%% fft value at the nyquist frequency
        shoulder = 0; %%% Degree of frequency overlap between neighboring bands (0 - 1)
        lowpass = []; %%% Lowpass cutoff
%        taperfun =[];
        taper = [];
%         taper = 'quadratic'; 
    end
    
    
    methods
    
        function me = dbt2(varargin)
            
            padding = 'time';
            i = 4;       
%             me.taperfun = @(x)x;
%            me.taperfun = @(x)(1-cos(x*pi))/2; % Function that defines tapering in the overlapping regions of the window
            me.taper = taper; %#ok<CPROP>
            
           while i < length(varargin)
              switch lower(varargin{i})
                  
                  case {'offset','highpass'}
                      me.offset = varargin{i+1};
                      i = i+1;
                  case 'band'
                      band = varargin{i+1};
                      me.offset = band(1);
                      me.lowpass = band(end);                      
                      i = i+1;
                  case 'padding'
                      padding = varargin{i+1};
                      i = i+1;
                  case 'shoulder'
                      me.shoulder = varargin{i+1};
                      i = i+1;
                  case 'lowpass'
                      me.lowpass = varargin{i+1};
                      i = i+1;
                 case 'taper'
                      me.taper = varargin{i+1};
                      i = i+1;
                  otherwise
                     error('Unrecognized keyword %s',varargin{i})
              end
              i = i+1;
           end
               
           if isempty(varargin)
               return
           end
           
           fs  = varargin{2};
           bw =  varargin{3};
%            me.decim = decim;
           
           bandwindows = 0:bw:fs/2-bw;
           
           fullsig  = varargin{1};
           nyq = bandwindows(end)+bw;
           
           n = length(fullsig);
           me.Norig = n;
           %Resample signal so that everything factors
           %Keeping signal duration fixed and allowing bandwidth and sampling
           %rate to vary.
           
           T = n./fs; % signal duration
           %%% K and M need to be integers
           
           K = ceil(1/2*(n./T - me.offset)/bw); 
           M = ceil(bw*T/2)*2;
           newnf = 2*ceil(K*M  );
           newn = 2*ceil(newnf/2 + me.offset*newnf/(fs-2*me.offset) );
           noffset = (newn-newnf)/2;
            
       
           %%% Adjusting bandwidth based on new signal length 
           newbw = 1/2*(newnf)/(T*K);           
           %%% Likewise for new sampling frequency
           newfs = newn./T;
           me.fullFS = newfs;
           nsh = floor(me.shoulder*newbw./newfs*newn/2)*2;
           me.shoulder = nsh*newfs./newn./newbw;
           me.fullN = newn+2*nsh;
         
           if isempty(me.lowpass)
              me.lowpass = newfs/2; 
           end
           
           nwin = K;
           winN = M+nsh;
         
%            noffset = floor(me.offset*newn/newfs);
           me.offset = noffset*newfs/newn;
           
            switch padding
               case 'time'
                  T = newn./fs;
                 fullsig(newn+2*nsh) = 0; 
              
                 
               case 'frequency'
                   
                   
           
           end
           
           F = fft(fullsig./sqrt(length(fullsig)));
          
           F([1 newnf/2+1])= F([1 newnf/2+1])/2;
           rsmat = noffset + repmat((0:nwin-1)*(winN-nsh),winN,1) + repmat([floor((winN)/2)+1:(winN),1:floor((winN)/2)]',1,nwin);
           tp = me.taper.make((1:1:nsh)/nsh); 
           invtaper = me.taper.make(1-(1:1:nsh)/nsh);

           Frs = double(F(rsmat));
%            Frs(floor(winN/2),1)=Frs(floor(winN/2),1 )/2;
%             Frs(winN/2-1,end)=Frs(winN/2-1,end )/2;
           %%% Now add the taper
%                Frs(end+(1-nsh:0),1:nwin-1) = diag(sparse(taper))*Frs(end+(1-nsh:0),1:nwin-1);
           if nsh>0
               Frs(winN/2+(1-nsh:0),1:nwin-1) = diag(sparse(tp))*Frs(winN/2+(1-nsh:0),1:nwin-1);
               Frs(winN/2+(1:nsh),2:nwin) = diag(sparse(invtaper))*Frs(winN/2+(1:nsh),2:nwin); % - Frs(end+(1-nsh:0),1:nwin-1);
           end
           winN = size(Frs,1);
%            end
%            me.nyqval = newF(newn/2+1);
           
%            Frs(winN*2,:) = 0;
           me.blrep = ifft(2*Frs)*sqrt(winN);
           
           me.sampling_rate = 2*winN/T;
           
                              
           me.bandwidth = newbw;
           me.bands = [me.offset:newbw:me.lowpass-newbw;(me.offset+newbw:newbw:me.lowpass)]';
           me.time = ((1:size(me.blrep,1))-.5)*T./size(me.blrep,1);
           w = ((0:K-1)+.5)*newbw + me.offset;          
           me.blrep(:,w>me.lowpass) = [];
           w(w>me.lowpass) = [];
           me.frequency = w;
           
        end

        %%%%

        function [data,fs] = signal(me,columnfilter,hilbert)
            
            %%% Reconstruct the signal from its band-limited analytic representation
            if nargin < 2 || isempty(columnfilter)
                mult = 1;
            elseif islogical(columnfilter)
                mult = diag(sparse(columnfilter));
            elseif min(size(columnfilter)) == 1
                
               mult = diag(sparse( ismember(1:size(me.blrep,2),columnfilter)));
               
            else
                mult = columnfilter;
            end
            
            if nargin < 3 || isempty(hilbert)
               hilbert = false; 
            end
            
            n = me.fullN;
            noffset = me.offset./me.fullFS*n;            
            F = fft(me.blrep )*mult/sqrt(size(me.blrep,1));          
            nsh = round(me.shoulder*me.bandwidth./me.fullFS*me.fullN);
            nnyq = size(F,1)/2;
            
            % Taper is normally defined so that h(k).^2 + h(k+bw).^2 = 1
            tp = me.taper.make((1:1:nsh)/nsh); 
            invtaper = me.taper.make(1-(1:1:nsh)/nsh);
            
            F = F([nnyq+1:end,1:nnyq],:);
            
            sh = diag(sparse(tp))*F(end-nsh+1:end,1:end-1);
            F(1:nsh,2:end) = diag(sparse(invtaper))*F(1:nsh,2:end)+sh;
            fend = F(end-nsh+1:end,end);
            F(end-nsh+1:end,:) = [];
             F(1:nsh,end+1) = fend;
            
            Ffull(noffset+(1:numel(F))) = F(:)*sqrt(me.fullN);
            Ffull(me.fullN) = 0;
            
            
%             Ffull(me.fullN/2+1) = imag(Ffull(1))/2;
%             Ffull(1) = real(Ffull(1))/2;
            
            if hilbert
                data = ifft(Ffull);
            else
                data = real(ifft(Ffull));
            end
            data = data(1:me.Norig);
            fs = me.fullFS;
            
        end
        
       
    end
end
    
