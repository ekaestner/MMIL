%% SL Analysis Script
clear; clc;

sbj_num = 1;

% Setting up Variables
sbj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
sbj  = sbj{sbj_num};

fprintf('Starting on sbject %s \n\n\n',sbj)

indir       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/indir/' sbj]);  
cln_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cln_dir/' sbj]); 
end_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/end_dir/' sbj]); 
tsk         = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/tsk/' sbj]);     
rmv_chn     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/rmv_chn/' sbj]),'uni',0)';
lne_nse     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/bse_frq/' sbj]),'uni',0)';
bse_frq     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/bse_frq/' sbj]),'uni',0)';
trl_fun = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/trialfun/' sbj]);

chn_nse_grp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' sbj]);
tmp = strfind(chn_nse_grp,';');
if any(~isempty(tmp{:}))
    for i = 1:numel(chn_nse_grp)
        chn_nse_grp_dat_nme(i) = str2num(chn_nse_grp{i}(tmp{1}+1:end))';
        chn_nse_grp{i} = str2num(chn_nse_grp{i}(1:tmp{1}-1))'; 
    end
else
    chn_nse_grp = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' sbj]),'uni',0)';
end

% In Script Switches 
plt_spc = [1 1];

%% Data Paths
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/';

for iIN = 1:numel(indir)
    inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
end

inpath = [];
for iIN = 1:numel(indir)
    inpath = [inpath inpath_holder.(tsk{iIN})];
end
inpath = sort(inpath);

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sll_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0);

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
if any(~cellfun(@isempty,strfind(end_dir,'set')) | ~cellfun(@isempty,strfind(end_dir,'edf'))); cfg.continuous = 'yes'; else cfg.continuous = 'no'; end
cfg.dataset    = inpath;
sll_dat        = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Remove Unimportant Channels - EJK
rmv_chn_uni = find(~cellfun(@isempty,strfind(sll_dat.(sll_dat.data_name{1}).label,'DC')))';
rmv_chn_uni = [rmv_chn_uni find(~cellfun(@isempty,strfind(sll_dat.(sll_dat.data_name{1}).label,'EKG')))'];

if ~isempty(rmv_chn_uni)
    cfg = [];
    cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{1}).label(rmv_chn_uni))'];
    sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Remove Line Noise
if ~isempty(lne_nse)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sll_dat      = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Examine Data for Noise
if plt_spc(1) == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'freq' '/' sbj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Look through bad channels (and good channels) of data


%% Remove identified Bad Channels
if ~isempty(rmv_chn)
    if numel(rmv_chn) == 1
    cfg = []; 
    cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{1}).label(rmv_chn))'];
    sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
    else
        for iRM = 1:numel(rmv_chn)
            cfg = [];
            cfg.data_name = iRM;
            cfg.channel = ['all',strcat('-',sll_dat.(sll_dat.data_name{iRM}).label(rmv_chn{iRM}))'];
            sll_dat = ft_func(@ft_preprocessing,cfg,sll_dat);
        end
    end
end

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sll_dat      = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Remove noise through removing common noise
% if ~isvar('chn_nse_grp_dat_nme'); chn_nse_grp_dat_nme = repmat({1:numel(tsk)},1,numel(chn_nse_grp)); end
for i = 1:numel(chn_nse_grp)
    if ~isempty(chn_nse_grp{i})
        cfg             = [];
        cfg.data_name   = chn_nse_grp_dat_nme{1};
        cfg.chn_nse_grp = chn_nse_grp(i);
        sll_dat = ft_func(@ft_remove_common_noise,cfg,sll_dat);
    end
end

%% Examine Data for Noise
if plt_spc(2) == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'after' '/' sbj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Epoch Data
cfg = [];
cfg.trialfun = trl_fun{1};
cfg.dataset  = inpath;
cfg.specific = {'dataset' ; 1:numel(inpath)};
cfg.minduration = 0.500;
cfg.pre         = 1.5;
cfg.pos         = 2.5;
cfg.ignore      = 1;
cfg.trg_chn     = 130:136;
trl = ft_func(@ft_definetrial,cfg,trl);

trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end

% - EJK Grab timings & events here

for iTS = 1:numel(trl_hld)
    
    blk_ind = regexpi(inpath{iTS},'\d_\d');
    
    blk = num2str(inpath{iTS}(blk_ind));
    run = num2str(inpath{iTS}(blk_ind+2));
    
    % - EJK Check timings & add events here
    
    lst = mmil_readtext(['/home/ekaestne/Desktop/N439_SL_Behav/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv']);
    trl_hld{iTS}(:,4) = cell2mat(lst(:,3));
    
end

cfg = [];
cfg.specific  = {'trl';[1:numel(inpath) 1:numel(inpath) 1:numel(inpath)]};
cfg.trl       = trl_hld;
sll_dat = ft_func(@ft_redefinetrial,cfg,sll_dat);

%% Combine Data
% combine trials
num_trl       = numel(inpath);
cfg           = [];
cfg.data_name = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2 ; ...
                3:2:num_trl 4:2:num_trl];
cfg.data_new  = 'yes';
cfg.methapp   = 'trials';
sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = [3:2:num_trl 4:2:num_trl];
sll_dat = ft_func([],cfg,sll_dat);

% combine Clin1 & Clin2
cfg           = [];
cfg.data_name = [1:2:numel(sll_dat.data_name) ; 2:2:numel(sll_dat.data_name)];
cfg.data_new  = 'yes';
cfg.methapp   = 'channels';
sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = [2:2:numel(sll_dat.data_name)];
sll_dat = ft_func([],cfg,sll_dat);

%% Create time-frequency data
% parameters for wavelet function
foi =  [2:12,14:2:24,25,30:5:55,70:10:170];    %frequency of interest
sf  =   [1,1,repmat(2,1,9),repmat(3,1,6),repmat(5,1,7),repmat(10,1,11)]; %specific frequency
gwidth = ones(size(foi))*pi; %wavelet

tic      
    cfg         = [];
    cfg.sve_loc = [outpath '/' 'freq' '/' sbj];
    cfg.foi     = foi;
    cfg.sf      = sf;
    cfg.gwidth  = gwidth;
    cfg.new_tme = [-0.5 1.5];
    mmil_calc_freq(cfg,sll_dat);
toc

%% Remove Artifacts
cfg = [];
cfg.dat_loc   = [outpath '/' 'freq' '/' sbj];
cfg.measures  = {'time' 'variance'};
cfg.frq_rng   = {[1 12] [14 24] [30 55] [70 170]};
cfg.thresh    = [0.985 0.985];
cfg.outdir    = [outpath '/' 'freq_rejection' '/' sbj '/' ];
cfg.prefix    = sll_dat.data_name{1};
cfg.pad       = 0;
cfg.plot      = 1;
mmil_reject_freq(cfg)

%% Events
cfg = [];
cfg.dat_loc     = [outpath '/' 'freq' '/' sbj];
cfg.add_eve     = {sll_dat.(sll_dat.data_name{1}).trialinfo};
cfg.alt_eve_lab = {'trialinfo'};
cfg.old_events  = {{[1 2] 3} {[1 2] 4}};
cfg.new_events  = {{101 102} {111 112}};
cfg.crt_alt_eve = {'vis_tot_nse' 'aud_tot_nse'};
mmil_freq_events(cfg);

%% Stats

%% Calculate Phase consistency
cfg = [];
cfg.dat_loc = [outpath '/' 'freq' '/' sbj];
cfg.itc_typ = 'Rayleigh';
cfg.eve_ind = [1 2 3 4];
mmil_calc_itc(cfg);

