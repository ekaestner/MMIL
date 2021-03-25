function FW_load_depth(fcfg)

%% Initial Variables
subj = fcfg.sbj_nme;
outpath = fcfg.out_pth; 

% In Script Switches 
fprintf('Starting to Load Subject %s \n\n\n',subj)

%% Data Paths
indir       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);  
cln_fld     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_fld');
cln_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir'); %mmil_readtext([fcfg.clr_fld '/cln_dir/' subj]); 
if numel(cln_dir) ~= numel(cln_fld) && numel(cln_fld)>0; cln_dir = repmat(cln_dir(1),1,numel(cln_fld)) ; end
end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]); 
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     

cfg = [];
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.cln_dir = cln_dir;
cfg.end_dir = end_dir;
cfg.tsk     = tsk;
[inpath,fwv_dat,trl] = mmil_find_files(cfg);

%% Load Data
epc_tme      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

for iL = 1:numel(fwv_dat.data_name)
    
    if ~strcmpi(end_dir,'edf')
        
        if strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            ttt = load(inpath{iL});
            fwv_dat.(fwv_dat.data_name{iL}) = ttt.epoch_data;
        elseif strcmpi(end_dir,'eeg') | strcmpi(end_dir,'.eeg')
            ttt = ts_load_data(inpath{iL});
            fwv_dat.(fwv_dat.data_name{iL}) = ttt;
        end
        
        fwv_dat.(fwv_dat.data_name{iL}) = mmil_format_mat(fwv_dat.(fwv_dat.data_name{iL}));
        clear ttt
        
    else
        
        cfg            = [];
        cfg.specific   = {'dataset';1:numel(inpath)};
        cfg.data_new   = 'yes';
        cfg.continuous = 'yes';
        cfg.dataset    = inpath;
        fwv_dat        = ft_func(@ft_preprocessing,cfg,fwv_dat);
        
        if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
            
            trialfun       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);
            
            cfg = [];
            if numel(trialfun) == 1; cfg.trialfun = trialfun{1}; else cfg.trialfun = trialfun; end
            cfg.dataset  = inpath;
            if numel(trialfun) == 1; cfg.specific = {'dataset' ; 1:numel(inpath)}; else cfg.specific = {'dataset' 'trialfun' ; 1:numel(inpath) 1:numel(trialfun)}; end
            cfg.minduration = 0.500;
            cfg.pre         = abs(epc_tme(1)-1);
            cfg.pos         = epc_tme(2)+1;
            cfg.evt         = 1:8;
            trl = ft_func(@ft_definetrial,cfg,trl);
                        
            trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end
                    
            save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
            
        else
            load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
        end
        
        cfg = [];
        cfg.specific  = {'trl';1:numel(inpath)};
        cfg.trl       = trl_hld;
        fwv_dat = ft_func(@ft_redefinetrial,cfg,fwv_dat);
        
    end
    
end

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
    fwv_dat = mmil_combine_data2(cfg,fwv_dat);
end

for iD = 1:numel(fwv_dat.data_name)
    if fwv_dat.(fwv_dat.data_name{iD}).fsample ~= 1./mean(diff(fwv_dat.(fwv_dat.data_name{iD}).time{1}))
        if numel(fwv_dat.(fwv_dat.data_name{iD}).time{1}(1)+1/fwv_dat.(fwv_dat.data_name{iD}).fsample:1/fwv_dat.(fwv_dat.data_name{iD}).fsample:fwv_dat.(fwv_dat.data_name{iD}).time{1}(end)) ~= size(fwv_dat.(fwv_dat.data_name{iD}).trial{1},2)
            dff_num = numel(fwv_dat.(fwv_dat.data_name{iD}).time{1}(1)+1/fwv_dat.(fwv_dat.data_name{iD}).fsample:1/fwv_dat.(fwv_dat.data_name{iD}).fsample:fwv_dat.(fwv_dat.data_name{iD}).time{1}(end)) - size(fwv_dat.(fwv_dat.data_name{iD}).trial{1},2);
            if dff_num<0
                fwv_dat.(fwv_dat.data_name{iD}).time = repmat({fwv_dat.(fwv_dat.data_name{iD}).time{1}(1)+(1/fwv_dat.(fwv_dat.data_name{iD}).fsample)*(1+dff_num):1/fwv_dat.(fwv_dat.data_name{iD}).fsample:fwv_dat.(fwv_dat.data_name{iD}).time{1}(end)},1,numel(fwv_dat.(fwv_dat.data_name{iD}).time));
            elseif dff_num>0
                 fwv_dat.(fwv_dat.data_name{iD}).time = repmat({fwv_dat.(fwv_dat.data_name{iD}).time{1}(1)-(1/fwv_dat.(fwv_dat.data_name{iD}).fsample)*(1+dff_num):1/fwv_dat.(fwv_dat.data_name{iD}).fsample:fwv_dat.(fwv_dat.data_name{iD}).time{1}(end)},1,numel(fwv_dat.(fwv_dat.data_name{iD}).time));
            end
        else
            fwv_dat.(fwv_dat.data_name{iD}).time = repmat({fwv_dat.(fwv_dat.data_name{iD}).time{1}(1)+1/fwv_dat.(fwv_dat.data_name{iD}).fsample:1/fwv_dat.(fwv_dat.data_name{iD}).fsample:fwv_dat.(fwv_dat.data_name{iD}).time{1}(end)},1,numel(fwv_dat.(fwv_dat.data_name{iD}).time));
        end
    end
end

%% Misc
if strcmpi(fcfg.sbj_nme,'NY192_FW')
    cfg = [];
    cfg.trials = [1:581 607:numel(fwv_dat.(fwv_dat.data_name{1}).time)];
    fwv_dat = ft_func(@ft_preprocessing,cfg,fwv_dat);
end

% if strcmpi(fcfg.sbj_nme,'NY240_FW')
%     trl_hld = fwv_dat.(fwv_dat.data_name{iD}).trialinfo;  
% end

%% Remove identified Unimportant & Bad Channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

fwv_dat = mmil_remove_channels(cfg,fwv_dat);

%% Combine Clinical/Day if necessary
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;

    fwv_dat = mmil_combine_data2(cfg,fwv_dat);
end

%% Check Depth Existence
dpt_num = mmil_readtext([fcfg.clr_fld '/' 'ele_idn' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '.csv']);

if sum(cell2mat(dpt_num(2:end,3))==2)~=0
    
%% Fix Labels
cfg = [];
cfg.ovr_wrt = 0;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;
cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');
cfg.dpt = 1;
fwv_dat = mmil_fix_labels(cfg,fwv_dat);

%% Depth Bi-Polar
cfg = [];
cfg.chn_str_cnt = 2;
fwv_dat.(fwv_dat.data_name{1}) = mmil_depth_bipolar(cfg,fwv_dat.(fwv_dat.data_name{1}));
    
%% Remove identified noise problems
bse_frq     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';

if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    fwv_dat      = ft_func(@ft_preprocessing,cfg,fwv_dat);
end

%% Filter Data for LFP
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 15;
cfg.new_suffix = 'lfp';
fwv_dat        = ft_func(@ft_preprocessing,cfg,fwv_dat);

%% Filter Data for HGP - EJK
if fwv_dat.(fwv_dat.data_name{1}).time{1}(end)-fwv_dat.(fwv_dat.data_name{1}).time{1}(2)>1
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
    cfg.toi=fwv_dat.(fwv_dat.data_name{1}).time{1};
    cfg.keeptrials = 'yes';
    fwv_dat        = ft_func(@mmil_hgp_freq_analysis,cfg,fwv_dat);
else
    cfg            = [];
    cfg.data_name  = 1;
    cfg.hilbert    = 'amp';
    cfg.freq_band  = {[70 170]};
    cfg.data_new   = 'yes';
    cfg.new_suffix = 'hgp_hlb';
    fwv_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,fwv_dat);
end

cfg             = [];
cfg.data_name   = numel(tsk)*2+1:numel(tsk)*3;
HGP_smooth_msec = 0.025;
w               = round(HGP_smooth_msec* fwv_dat.(fwv_dat.data_name{1}).fsample); 
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,round(0.016* fwv_dat.(fwv_dat.data_name{1}).fsample));
cfg.window      = gauss_y/sum(gauss_y);
fwv_dat         = ft_func(@ft_window_data,cfg,fwv_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(tsk);
cfg.rmfield   = 'yes';
fwv_dat       = ft_func([],cfg,fwv_dat);

%% Baseline Data
bse_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_tme'); bse_tme = bse_tme{1};

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = bse_tme;
fwv_dat            = ft_func(@ft_baseline,cfg,fwv_dat);

%% Automatic rejection & Apply Rejections
epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
rjt_plt = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rjt_plt'); rjt_plt = rjt_plt{1};

cfg = [];
cfg.measures  = {'time' 'time-1' 'variance'};
cfg.thresh    = [0.95 0.95 0.95];
cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
cfg.prefix    = fwv_dat.data_name;
cfg.specific  = {'prefix';1:numel(fwv_dat.data_name)};
cfg.pad       = fwv_dat.(fwv_dat.data_name{1}).time{1}(end)-epc_tme(2);
cfg.plot      = 0;
fwv_dat       = ft_func(@auto_rej,cfg,fwv_dat);

cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
fwv_dat     = ft_func(@ft_apply_rejection,cfg,fwv_dat);

%% Remove Padding
epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

cfg = [];
cfg.latency = epc_tme;
fwv_dat     = ft_func(@ft_selectdata,cfg,fwv_dat);

%% Save Data
cfg = [];
cfg.str_nme  = 'fwv_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data_depth.mat'];
ft_func([],cfg,fwv_dat);

end

end
