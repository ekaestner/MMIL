function SL_initial_load2(fcfg)

%% Initial Variables
subj = fcfg.sbj_nme;
outpath = fcfg.out_pth; 

% In Script Switches 
fprintf('Starting to Load Subject %s \n\n\n',subj)

%% Data Paths
indir       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);  
cln_fld     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_fld');
cln_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir'); %mmil_readtext([fcfg.clr_fld '/cln_dir/' subj]); 
end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]); 
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     
tot_tsk     = unique(tsk);

cfg = [];
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.cln_dir = cln_dir;
cfg.end_dir = end_dir;
cfg.tsk     = tsk;
[inpath,sll_dat,trl] = mmil_find_files(cfg);

%% Load Visual & Auditory data
epc_tme      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
cfg.continuous = 'yes';
cfg.dataset    = inpath;
sll_dat        = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Epoch Data
trialfun         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun');
epc_tme          = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
ignore      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'ignore');

if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
    
    cfg = [];
    if numel(trialfun) == 1; cfg.trialfun = trialfun{1}; else cfg.trialfun = trialfun; end
    cfg.dataset  = inpath;
    cfg.ignore   = ignore; if numel(cfg.ignore)~=numel(cfg.dataset); cfg.ignore = repmat(cfg.ignore,1,numel(cfg.dataset)); end
    if numel(trialfun) == 1; cfg.specific = {'dataset' 'ignore' ; 1:numel(inpath) 1:numel(inpath)}; else cfg.specific = {'dataset' 'trialfun' 'ignore' ; 1:numel(inpath) 1:numel(inpath) 1:numel(trialfun)}; end
    cfg.minduration = 0.500;
    cfg.pre         = abs(epc_tme(1)-1);
    cfg.pos         = epc_tme(2)+1;
    trl = ft_func(@ft_definetrial,cfg,trl);
        
    trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end    
    
    save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
else
    load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
end

cfg = [];
cfg.specific  = {'trl';1:numel(inpath)};
cfg.trl       = trl_hld;
sll_dat = ft_func(@ft_redefinetrial,cfg,sll_dat);

%% First Combine
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb'); cmb = cmb{1};

if isempty(cmb)
    inpath{:}
    disp(['Elaborate on which file(s) go together (cmb): ' ' '])
    pause();
    pause(5);
    cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb'); cmb = cmb{1};
end

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    sll_dat = mmil_combine_data2(cfg,sll_dat);
end

%% Misc

%% Remove identified Unimportant & Bad Channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

sll_dat = mmil_remove_channels(cfg,sll_dat);

%% Examine Data for Noise - Initial
beg_plt_spc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'beg_plt_spc'); beg_plt_spc = beg_plt_spc{1};

if beg_plt_spc %&& ~exist([outpath '/' 'spectrum' '/' subj '/' 'Initial' '/' sll_dat.data_name{1} '_1.png'],'file');
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'Initial' ];
    cfg.prefix = sll_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(sll_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Remove noise through removing common noise
chn_nse_grp = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');
nse_val         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'nse_val'); nse_val = nse_val{1};

% Make temporary ECOG/DEPTH Split
for iD = 1:numel(sll_dat.data_name); sll_dat.(sll_dat.data_name{iD}).cfg.alt_lab.label = sll_dat.(sll_dat.data_name{iD}).label; end
scfg = [];
scfg.dat_nme = sll_dat.data_name;
scfg.clr_fld = fcfg.clr_fld;
scfg.sbj_nme = fcfg.sbj_nme;
scfg.alt_lab = 'label';
scfg.specific = {'dat_nme' ; 1:numel(sll_dat.data_name)};
scfg.add_nme  = '_tmp';
ft_func(@mmil_create_depth2,scfg,sll_dat);

if nse_val   
    
    if strcmpi(chn_nse_grp{1},'split'); chn_nse_grp = repmat(chn_nse_grp,1,numel(sll_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(sll_dat.data_name); pre_fix{iD}{1} = strcat(sll_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(sll_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(sll_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(sll_dat.data_name{iD},'noisy_depth'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(sll_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.sbj_nme     = fcfg.sbj_nme;
    cfg.clr_fld     = fcfg.clr_fld;
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = sll_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(sll_dat.data_name) 1:numel(sll_dat.data_name) 1:numel(sll_dat.data_name)};
    cfg.out_dir     = [fcfg.out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    sll_dat = ft_func(@ft_remove_common_noise2,cfg,sll_dat);    
    
    if nse_val == 1
        cfg.clr_fld = fcfg.clr_fld;
        cfg.sbj_nme = fcfg.sbj_nme;
        cfg.sub_fld_ind = num2cell(1:numel(sll_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(sll_dat.data_name) 1:numel(sll_dat.data_name)};
        sll_dat = ft_func(@mmil_update_cmn_nse2,cfg,sll_dat);
    end
    
end

%% Remove identified noise problems
bse_frq     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';

if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sll_dat      = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Examine Data for Noise
end_plt_spc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_plt_spc'); end_plt_spc = end_plt_spc{1};

if end_plt_spc == 1 %&& ~exist([outpath '/' 'spectrum' '/' subj '/' 'After' '/' sll_dat.data_name{1} '_1.png'],'file');
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'After'];
    cfg.prefix = sll_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(sll_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Combine Clinical/Day if necessary
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;

    sll_dat = mmil_combine_data2(cfg,sll_dat);
end

%% Fix Labels
cfg = [];
cfg.ovr_wrt = 0;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;
cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');

sll_dat = mmil_fix_labels(cfg,sll_dat);

%% Setup events
cfg = [];
cfg.sbj    = subj;
cfg.clr_fld = fcfg.clr_fld;
cfg.inpath = inpath;
trl_hld = SL_fix_trl(cfg,trl_hld);

%% Filter Data for LFP
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 15;
cfg.new_suffix = 'lfp';
sll_dat        = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Filter Data for HGP - EJK
cfg=[];
cfg.data_name  = 1:numel(tsk);
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
cfg.foi    = [70 80 90 100 110 130 140 150 160 170];    %frequency of interest
cfg.sf     = [repmat(10,1,numel(cfg.foi))]; %specific frequency
cfg.width  = cfg.foi./10;
cfg.gwidth = ones(size(cfg.foi))*pi; %wavelet
cfg.keeptrials = 'yes';
cfg.method = 'tfr';
cfg.toi=sll_dat.(sll_dat.data_name{1}).time{1};
cfg.keeptrials = 'yes';
sll_dat        = ft_func(@mmil_hgp_freq_analysis,cfg,sll_dat);

cfg             = [];
cfg.data_name   = numel(tsk)*2+1:numel(tsk)*3;
HGP_smooth_msec = 0.025;
w               = round(HGP_smooth_msec* sll_dat.(sll_dat.data_name{1}).fsample); 
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,round(0.016* sll_dat.(sll_dat.data_name{1}).fsample));
cfg.window      = gauss_y/sum(gauss_y);
sll_dat         = ft_func(@ft_window_data,cfg,sll_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(tsk);
cfg.rmfield   = 'yes';
sll_dat       = ft_func([],cfg,sll_dat);

%% Baseline Data
bse_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_tme'); bse_tme = bse_tme{1};

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = bse_tme;
sll_dat            = ft_func(@ft_baseline,cfg,sll_dat);

%% Automatic rejection & Apply Rejections
epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
rjt_plt = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rjt_plt'); rjt_plt = rjt_plt{1};

cfg = [];
cfg.measures  = {'time' 'time-1' 'variance'};
cfg.thresh    = [0.98 0.98 0.98];
cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
cfg.prefix    = sll_dat.data_name;
cfg.specific  = {'prefix';1:numel(sll_dat.data_name)};
cfg.pad       = sll_dat.(sll_dat.data_name{1}).time{1}(end)-epc_tme(2);
cfg.plot      = rjt_plt;
sll_dat       = ft_func(@auto_rej,cfg,sll_dat);

cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sll_dat     = ft_func(@ft_apply_rejection,cfg,sll_dat);

%% Remove Padding
epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

cfg = [];
cfg.latency = epc_tme;
sll_dat     = ft_func(@ft_selectdata,cfg,sll_dat);

%% Save Data
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sll_dat);

end
