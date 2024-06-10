function pk2contact = simanneal(D,T,TXYZs,XYZ,pk2contact)

% pk2contact = simanneal(D,T,TXYZs,XYZ)
% 
% Stimulate annealing to match distances in matrix D with template T.
% D - N x N matrix of inter-peak distances.
% T - M x M matrix of inter-electrode distances based on contacts. Wherever
%           the distance is unknown (eg grid to strip), put a value of -1.
% TXYZs - X-Y-Z coordinates of electrodes in template (needed only for plotting)
% XYX - X-Y-Z coordinates of peaks (needed only for plotting)
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/simanneal.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

costscale = 5;   %Parameters of the cost function
costreg = .5;     

wgtscale = 30; %differences are weighted by a function of the distances in the template
               %Errors in large differences will be treated as less important
               %than errors in small differences, given that the grids can
               %deform

               
plotprog = true; %Generates a plot showing the fit if true

allow_multiple_peaks = true;

add_null_contact = true; % Add a pseudo contact for unmapped peaks

if add_null_contact
     T(end+1,:) = -1; 
     T(:,end+1) = -1;
     TXYZs{end+1} = [0 0 0 0];
     
end

% WT = exp(-T/wgtscale).*(T>=0); %Weight matrix


%%% Determines how the cost function is computed for simulated
%%% annealing. It returns a function of the error between the
%%% template and the observed inter-peak distances.

%%% This is one possible cost function (maybe others would work better?)
costfun = @(t,d) -log(abs(t-d)./costscale+costreg).*exp(-t/wgtscale).*(t>=0);
% proxscale = 10;
% costfun = @(t,d) -log(abs(t-d)./costscale+costreg).*exp(-t/wgtscale).*(t>=0) - exp(-d/proxscale).*(t<0); %penalize nearby peask that are not mapped to template

% % costfun = @(t,d) -((t-d).^2/costscale^2 - 1)/costscale^2.*exp(-((t)./(wgtscale)).^2).*(t>=0); %2nd derivitave of the guassian
% % costfun = @(t,d) -log(abs(t-d)./costscale+costreg).*exp(-t/wgtscale).*(t>=0)./sum(exp(-t(:)/wgtscale).*(t(:)>=0));

% %correlation
% costfun = @(t,d) (sum(t.*d.*(t>=0),2) - sum(t.*(t>=0),2).*sum(d.*(t>=0),2)./sum(t>=0,2))./...
%             sqrt((sum(t.^2.*(t>=0),2)-sum(t.*(t>=0),2).^2./sum(t>=0,2)).*(sum(d.^2.*(t>=0),2)-sum(d.*(t>=0),2).^2./sum(t>=0,2)));
        


ncontacts = size(T,1); %number of contacts in the template
npeaks = size(D,1);    % number of extracted peaks


%%% pk2contact contains the map from peaks to contacts -- initially this is
%%% random. Note that we allow multipe peaks to map to a single contact
%%% (but not vice versa)

% %First match by histogram
% [Dhist,b] = hist(D,0:1:200);
% Thist = hist(T,[-1,0:1:200]);
% Thist = Thist(2:end,:);
% 
% DThR = corr(Dhist,Thist); %Correlation between the histograms of starting points

if nargin < 5
    if allow_multiple_peaks
        pk2contact = ceil(rand(1,npeaks)*ncontacts); 
    %     rthresh = .1;
    %     [mx, pk2contact] = max(DThR,[],2);
    %     pk2contact(mx< rthresh) = ncontacts;

    else
         pk2contact = cat(2,1:min(ncontacts,npeaks),ones(1,max(npeaks-ncontacts,0))*ncontacts); 

         shuffle = @(a) a(randperm(length(a)));
         pk2contact = shuffle(pk2contact);
    end
    locked = zeros(size(pk2contact));
else
    
    locked = pk2contact ~=0;
    
    pk2contact(pk2contact==0) = ncontacts;
    if add_null_contact
        pk2contact(ncontacts) = ncontacts;
    end
end

% maxtemp = 5000;
maxtemp = 2; %Maximum inverse temperature in the simulated annealing-- the best value depends on the cost function
starttemp = 0; %Initial inverse temperature (0 - complete randomness)
nstep = 2e2;  %Number of iterations between starttemp and maxtemp


invtemps = linspace(starttemp,maxtemp,nstep); %The inverse temperatures
invtemps = (invtemps./maxtemp).^2*maxtemp; %quadratic change


% nullp = linspace(.5,0,nstep); % The probability of removing a link 
%                               % This is added to reduce the probability of
%                               % a poor connection permanently displacing a
%                               % better connection.
nullp = zeros(1,nstep);

Ttemp = T(pk2contact,pk2contact); %Reorder the template to match the mapping in pk2contact

C =  costfun(Ttemp,D);  %Initial matrix of costs
E = sum(C(:));


% Es0 = [];

% tempicon = {'[X^b','[X^<','[X^O','[.^O','[:^O','[:^o','[:^(','[:^|','[:^)','[:^D','[;^D'};

if plotprog   %Initialize plot
    Es = [];
    fig = figure;
    TXY = cat(1,TXYZs{:});
%     null = size(TXY,1)+1;
%     TXY(null,:) = 0;
    subplot(1,2,2)
    plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ko')
    hold on
    getpl = ~add_null_contact  | pk2contact~=ncontacts;
    plind = pk2contact(getpl);
    pl = plot3([XYZ(getpl,1),TXY(plind,1)]',[XYZ(getpl,2),TXY(plind,2)]',[XYZ(getpl,3),TXY(plind,3)]','r');
    plot3(TXY(1:end-1,1),TXY(1:end-1,2),TXY(1:end-1,3),'+')
    grid on
    axis equal
    drawnow
    
    ui(1) = uicontrol(fig,'style','pushbutton','units','normalized','position',[.5 .9 .1 .05 ],'string','Halt','Callback',@haltcallback);
    ui(2) = uicontrol(fig,'style','pushbutton','units','normalized','position',[.6 .9 .1 .05 ],'string','Restart','Callback',@restartcallback);
    ui(3) = uicontrol(fig,'style','pushbutton','units','normalized','position',[.7 .9 .1 .05 ],'string','Interpolate','Callback',@interpolate);
    ui(4) = uicontrol(fig,'style','checkbox','units','normalized','position',[.7 .9 .1 .05 ],'string','Show links','value',1);


end

setappdata(fig,'halt',false);
setappdata(fig,'restart',false);

rotor = {'-','\\','|','/'};

fprintf('\n-')
% for t = 1:length(invtemps)  %Step through each temperature
t= 0;
halt =  false;
while t <length(invtemps) && ~halt  %Step through each temperature
    
    t = t+1;
    
    fprintf('%4i',t)
%     Es0(end+1) = E;
    
    ks = setdiff(randperm(npeaks),find(locked));
%     for k = 1:npeaks
    for ki = 1:length(ks)
        %%% For each peak, evaluate the incremental change (DE) in the cost function
        %%% from switching to each of the possible contacts in the template, or the null contact. 
        %%% Then choose one of the contacts using softmax( DE/temperature)
        
%         setind =  setdiff(1:npeaks,k);
%         R = costfun(T(:,pk2contact(setind)),repmat(D(k,setind),ncontacts,1));
        fprintf(['\b',rotor{mod(ki,length(rotor))+1}])
        
        k = ks(ki);
            
        sever = rand< nullp(t);
        
        if sever
           pk2contact(k) = ncontacts;
            E = E - 2*sum(C(k,:)) + C(k,k); 
            C(k,:) = 0;
            C(:,k) = 0;
        else
            R = costfun(T(:,pk2contact),repmat(D(k,:),ncontacts,1)); %Cost of switching peak k to each of the contacts        
            DE =  2*sum(R,2) - 2*sum(C(k,:)) + C(k,k) - R(:,k); %Incremental cost

            if ~allow_multiple_peaks
    %             R = costfun(T(pk2contact,pk2contact),repmat(D(k,:),npeaks,1)); %Cost of switching peak k to each of the contacts        
                R2 =  diag(pk2contact < ncontacts)*costfun(repmat(T(pk2contact(k),pk2contact),npeaks,1),D);

                DEditch = DE;
                DEditch(pk2contact) = DE(pk2contact) - 2*sum(C,2) +diag(C);
                DEswap = DEditch;
                DEswap(pk2contact) =  DEditch(pk2contact) + 2*sum(R2,2) - R2(:,k);


                DE = [DEswap; DEditch];
            end

            DEc = DE-max(DE);           %Make exp more numerically stable by subtracting the max 


            W = exp(DEc*invtemps(t));  %Relative weight of each alternative

            rp = randpick(W);           %Randomply pick one of the alternative with probability proportional to W.

            if ~allow_multiple_peaks
                if rp < ncontacts
                        swk = pk2contact==rp;
                        pk2contact(swk) = pk2contact(k);            
                 elseif rp > ncontacts 
                        rp = rp-ncontacts;
                        swk = pk2contact==rp;
                        pk2contact(swk) = ncontacts;                    
                end
            end
    %         if sum(swk)>1&& rp < ncontacts, keyboard, end
             pk2contact(k) = rp;     %Updating the mapping with the chosen value    


    %         C(k,setind) = R(rp,:);
    %         C(setind,k) = R(rp,:);

            C(k,:) = R(rp,:);
            C(:,k) = R(rp,:);
            if ~allow_multiple_peaks && rp < ncontacts
                C(swk,:) = R2(swk,:);
                C(:,swk) = R2(swk,:);
            end        
              E = E + DE(rp);   %Update the cost function
        end
    end
    fprintf('\b')
        
    
    if plotprog 
        Es(end+1) = E;
        
%         plot(Es0);
        subplot(2,2,1)  %Plot the next value of the cost function
        cla
        plot(Es,'r');
%          tt = title(sprintf('Inverse Temperature: %s',tempicon{ceil(t./length(invtemps)*length(tempicon))}));
%          tt = title(sprintf('%s',tempicon{ceil(t./length(invtemps)*length(tempicon))}),'interpreter','none');
%         set(tt,'rotation',-90,'position',get(tt,'position')+[0 150 0])
        subplot(1,2,2)  %Plot the current mapping
        delete(pl)
%        pl = plot3([XYZ(:,1),TXY(pk2contact,1)]',[XYZ(:,2),TXY(pk2contact,2)]',[XYZ(:,3),zeros(size(XYZ,1),1)]','r');
      getpl = (~add_null_contact |  pk2contact~=ncontacts) & ~locked;
      plind = pk2contact(getpl);
      
      if get(ui(4),'value')
          pla = plot3([XYZ(getpl,1),TXY(plind,1)]',[XYZ(getpl,2),TXY(plind,2)]',[XYZ(getpl,3),TXY(plind,3)]','r-');
      end
      pl = [pla;plot3([XYZ(getpl,1),TXY(plind,1)]',[XYZ(getpl,2),TXY(plind,2)]',[XYZ(getpl,3),TXY(plind,3)]','g.')];
      
      getpl = locked;
      plind = pk2contact(getpl);
      pl = [pl;plot3([XYZ(getpl,1),TXY(plind,1)]',[XYZ(getpl,2),TXY(plind,2)]',[XYZ(getpl,3),TXY(plind,3)]','k-')];
      pl = [pl;plot3([XYZ(getpl,1),TXY(plind,1)]',[XYZ(getpl,2),TXY(plind,2)]',[XYZ(getpl,3),TXY(plind,3)]','g.')];
      
      
      figure(fig)
       
      halt = getappdata(fig,'halt');
      restart = getappdata(fig,'restart');
      if restart
          fprintf('\nRestarting...')
          t = 0;
         setappdata(fig,'restart',false)
      end
    end

    
    fprintf('\b\b\b\b')
end

pk2contact(add_null_contact & (pk2contact == ncontacts) ) = 0;



%%%%%%%%%%

function haltcallback(h,e,handle)
fig  = get(h,'parent');
setappdata(fig,'halt',true);

%%%%%%%%%%

function restartcallback(h,e,handle)

fig  = get(h,'parent');
setappdata(fig,'restart',true);



function interpolate(h,e,handles)

