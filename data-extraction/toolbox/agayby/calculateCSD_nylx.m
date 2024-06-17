function [varargout] = calculateCSD_nylx(data,varargin)

if ~isempty(varargin)
    brokenElektrode = varargin{1};
    electrodeOrder  = varargin{2};
else
     brokenElektrode = logical(ones(1,length(data.label)));
     electrodeOrder  =  1:length(data.label);
end

%filter data for cleaner CSD
for t = 1:length(data.trial)
    data.trial{t}    = ft_preproc_dftfilter(data.trial{t},data.fsample, 50);
%     data.trial{t}    = ft_preproc_bandpassfilter(data.trial{t}, data.fsample, [2 50]);
    data.trial{t}    = ft_preproc_lowpassfilter(data.trial{t}, data.fsample, 50);
end

%redifinetrials
% cfg                 = [];
% cfg.toilim          = [-0.1 0.1];
% data                = ft_redefinetrial(cfg,data);

%% iSplineCSD
cfg                 = [];
cfg.vartrllength    = 1;
cfg.channel         = 1:length(data.label);
CSDdata             = ft_timelockanalysis(cfg, data);
% 
cfg = [];
cfg.baseline = [-0.50 0];
CSDdata = ft_timelockbaseline(cfg,CSDdata);
cfg = [];
cfg.latency = [-0.5 1.2];
CSDdata = ft_selectdata(cfg,CSDdata);


for i = 1:length(electrodeOrder)
    reOrderedCSD(i,:) = CSDdata.avg(electrodeOrder(i),:);
end

CSDdata.avg         = reOrderedCSD(brokenElektrode,:);

dt                  = 0.5;
diam                = 0.5;
cond                = 0.4;
cond_top            = 0.4;
el_d                 =  0.2;
el_pos              = (0.7:el_d:0.7+(length(data.label)*el_d))*1e-3;
% el_pos              = 0.1:0.1:2;
el_pos              = el_pos(brokenElektrode);

gauss_sigma         = 0.1*1e-3;
filter_range        = 5*gauss_sigma;

Fcs                 = F_cubic_spline(el_pos,diam,cond,cond_top);
[zs,CSD_cs]         = make_cubic_splines(el_pos,CSDdata.avg,Fcs);
[zs,CSD_cs]         = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
iCSD                = CSD_cs/10^6;

%% classic CSD
dt                  = 0.5;
b0                  = 0.54;
b1                  = 0.22;
cond                = 0.4;
pot                 = CSDdata.avg;

% el_pos              = (0.3:0.1:0.3+(length(data.label)*0.1))*1e-3;
% el_pos              = el_pos(brokenElektrode);

N                   = length(el_pos);
h                   = mean(diff(el_pos));
[m1,m2]             = size(pot);

el_pos_plot         = el_pos;
pot(1,:)            = pot(1,:);
pot(2:m1+1,:)       = pot;
pot(m1+2,:)         = pot(m1,:);

CSD                  = -cond*D1(length(pot(:,1)),h)*pot;

[n1,n2]             = size(CSD);
CSD_add(1,:)        = zeros(1,n2);   %add top and buttom row with zeros
CSD_add(n1+2,:)     = zeros(1,n2);
CSD_add(2:n1+1,:)   = CSD;           %CSD_add has n1+2 rows
CSD                 = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows

%%
varargout{1}        = CSD;
varargout{2}        = iCSD;
varargout{3}        = CSDdata.time;
end