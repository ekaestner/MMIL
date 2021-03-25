%% SL Analysis Script
function SL_initial_load(fcfg)

subj = fcfg.sbj_nme;

fprintf('Starting to Load Subject %s \n\n\n',subj)

indir       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);  
cln_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir'); %mmil_readtext([fcfg.clr_fld '/cln_dir/' subj]); 
end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]); 
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     
rmv_chn     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';
bse_frq     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';
ignore     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'ignore'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';
if any(~cellfun(@isempty,strfind(end_dir,'edf'))); 
    trl_fun = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun'); end %trl_fun = mmil_readtext([fcfg.clr_fld '/trialfun/' subj]); end

chn_nse_grp = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');

% In Script Switches 
plt_spc = 0;
plt_rej = 0;

%% Data Paths
outpath = fcfg.out_pth; 

idn_fle = {};
for iIN = 1:numel(indir)
    if any(~cellfun(@isempty,strfind(end_dir,'eeg')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN}); if numel(inpath_holder.(tsk{iIN})) > 1; inpath_holder.(tsk{iIN})(2) = []; end;
    elseif any(~cellfun(@isempty,strfind(end_dir,'edf')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
        num_fle.(tsk{iIN}) = numel(inpath_holder.(tsk{iIN}));
        for i = 1:num_fle.(tsk{iIN}); idn_fle{end+1} = tsk{iIN}; end
    elseif any(~cellfun(@isempty,strfind(end_dir,'set')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
        num_fle.(tsk{iIN}) = numel(inpath_holder.(tsk{iIN}));
        for i = 1:num_fle.(tsk{iIN}); idn_fle{end+1} = tsk{iIN}; end
    elseif any(~cellfun(@isempty,strfind(end_dir,'mat')));
        inpath(iIN) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
    end
end

if ~isvar('inpath')
    inpath = [];
    for iIN = 1:numel(indir)
        inpath = [inpath inpath_holder.(tsk{iIN})];
    end
end

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
bcc_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

if any(~cellfun(@isempty,strfind(end_dir,'edf'))); inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0); end

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
cfg.continuous = 'yes';
cfg.dataset    = inpath;
bcc_dat        = ft_func(@ft_preprocessing,cfg,bcc_dat);

%% Remove identified Unimportant & Bad Channels
% identified channels
if any(~cellfun(@isempty,rmv_chn))
    if numel(rmv_chn)==1
        cfg = [];
        cfg.channel = ['all',strcat('-',bcc_dat.(bcc_dat.data_name{1}).label(rmv_chn{1}))'];
        bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);
    else
        for iRM = 1:numel(rmv_chn); tmp_rmv{iRM} = ['all',strcat('-',bcc_dat.(bcc_dat.data_name{iRM}).label(rmv_chn{iRM}))']; end
        cfg = [];
        cfg.channel = tmp_rmv;
        cfg.specific = {'channel' ; 1:numel(tmp_rmv)};
        bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);
    end
end

% Remove Unimportant Channels
rmv_chn_uni = find(~cellfun(@isempty,strfind(bcc_dat.(bcc_dat.data_name{1}).label,'EKG')))';
rmv_chn_uni = [rmv_chn_uni find(~cellfun(@isempty,strfind(bcc_dat.(bcc_dat.data_name{1}).label,'DC')))'];

if ~isempty(rmv_chn_uni)
    cfg = [];
    cfg.channel = ['all',strcat('-',bcc_dat.(bcc_dat.data_name{1}).label(rmv_chn_uni))'];
    bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);
end

cfg = [];
cfg.channel = ['all',strcat('-',intersect(bcc_dat.(bcc_dat.data_name{1}).label,bcc_dat.(bcc_dat.data_name{2}).label)')];
bcc_dat = ft_func(@ft_preprocessing,cfg,bcc_dat);

%% Examine Data for Noise - Initial
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'Initial' ];
    cfg.prefix = bcc_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(bcc_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,bcc_dat);
end

%% Remove noise through removing common noise
if ~isempty(chn_nse_grp{1})
    
    for iEJ = 1:numel(chn_nse_grp_nme); for iCN = 1:numel(chn_nse_grp_nme{iEJ}); pre_fix{iEJ}{iCN} = strcat(bcc_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN}); end; end
    
    cfg             = [];
    cfg.data_name   = 1:2;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix = pre_fix;
    cfg.specific  = {'pre_fix' 'chn_nse_grp'; 1:numel(bcc_dat.data_name) 1:numel(bcc_dat.data_name)};
    cfg.out_dir = [fcfg.fle_out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn = 1;
    bcc_dat = ft_func(@ft_remove_common_noise,cfg,bcc_dat);
       
%     cfg.app_all      = [1 1 ; 2 2];
%     mmil_app_all_nse(cfg,bcc_dat);
    
    cfg             = [];
    cfg.data_name   = 3:numel(bcc_dat.data_name);
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix = pre_fix;
    cfg.specific  = {'pre_fix' 'chn_nse_grp'; 3:numel(bcc_dat.data_name) 3:numel(bcc_dat.data_name)};
    cfg.out_dir = [fcfg.fle_out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn = 1;
    bcc_dat = ft_func(@ft_remove_common_noise,cfg,bcc_dat);
    
    cfg.data_name   = 1:numel(bcc_dat.data_name);
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.sub_fld_ind = num2cell(1:numel(bcc_dat.data_name));
    cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(bcc_dat.data_name) 1:numel(bcc_dat.data_name)};
    tmp = ft_func(@mmil_update_cmn_nse,cfg,bcc_dat);
    clear tmp
    
end

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    bcc_dat      = ft_func(@ft_preprocessing,cfg,bcc_dat);
end

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.data_name   = 2;
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'After'];
    cfg.prefix = bcc_dat.data_name;
    cfg.specific  = {'prefix'; 2:numel(bcc_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,bcc_dat);
end

%% Trial Function & Concatenate
if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
    
    cfg = [];
    cfg.trialfun = trl_fun{1};
    cfg.ignore   = ignore;
    cfg.dataset  = inpath;
    cfg.specific = {'dataset' 'ignore'; 1:numel(inpath) 1:numel(inpath)};
    cfg.minduration = 0.500;
    cfg.pre         = 1.5;
    cfg.pos         = 2.5;
    trl = ft_func(@ft_definetrial,cfg,trl);
       
    trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end
    
    for iTS = 1:numel(trl_hld)
        
        if strcmpi(subj,'NY439_SL')
            blk_ind = regexpi(inpath{iTS},'\d_\d');
        elseif strcmpi(subj,'NY537_SL') || strcmpi(subj,'NY540_SL') || strcmpi(subj,'NY523_SL')
            blk_ind = regexpi(inpath{iTS},'\d-\d');
        end
        
        blk = num2str(inpath{iTS}(blk_ind));
        run = num2str(inpath{iTS}(blk_ind+2));
        
        % - EJK Check timings & add events here
        if strcmpi(subj,'NY439_SL')
            lst_loc = ['/home/ekaestne/Desktop/N439_SL_Behav/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
            lst = mmil_readtext(lst_loc);
        elseif strcmpi(subj,'NY537_SL')
            lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY537/NY537_SL/SL_10_13_15_NY537/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
            lst = mmil_readtext(lst_loc);
        elseif strcmpi(subj,'NY540_SL')
            lst_loc = ['/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY537/NY537_SL/SL_10_13_15_NY537/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv'];
            lst = mmil_readtext(lst_loc);
        elseif strcmpi(subj,'NY523_SL')
            lst = mmil_readtext('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_1_1.csv');
            lst2 = mmil_readtext('/space/mdeh5/1/halgdev/incoming/iEEG_NYU/NY523/NY523_SL/SL_2_11_16_NY523/lists/SL_Vis1st_1_2_1.csv');
            lst = [lst ; lst2];
        end
        
        
        trl_hld{iTS}(:,4) = cell2mat(lst(:,3));
        
    end
    
    save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
else
    load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
end

cfg = [];
cfg.specific  = {'trl';1:numel(inpath)};
cfg.trl       = trl_hld;
bcc_dat = ft_func(@ft_redefinetrial,cfg,bcc_dat);

% Combine as necessary     
if any(~cellfun(@isempty,strfind(bcc_dat.data_name,'lin2')))
    
    ttt.clin1 = string_find(bcc_dat.data_name,'lin1');
    ttt.clin2 = string_find(bcc_dat.data_name,'lin2');
    
    for iCL = 1:2
       
        ttt.(['clin' num2str(iCL)]) = string_find(bcc_dat.data_name,['lin' num2str(iCL)]);
        
        cfg           = [];
        cfg.data_name = [repmat(ttt.(['clin' num2str(iCL)])(1),1,numel(ttt.(['clin' num2str(iCL)]))-1) ; ttt.(['clin' num2str(iCL)])(2:end)];
        cfg.data_new  = 'yes';
        cfg.methapp   = 'trials';
        bcc_dat = ft_func(@ft_appenddata,cfg,bcc_dat);
        
        cfg = [];
        cfg.rmfield   = 'yes';
        cfg.data_name =  ttt.(['clin' num2str(iCL)])(2:end);
        bcc_dat = ft_func([],cfg,bcc_dat);
               
    end
        
    cfg           = [];
    cfg.data_name = [1 ; 2 ];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'channels';
    bcc_dat = ft_func(@ft_appenddata,cfg,bcc_dat);
    
    cfg = [];
    cfg.rmfield   = 'yes';
    cfg.data_name = 2;
    bcc_dat = ft_func([],cfg,bcc_dat);
else
    cfg           = [];
    cfg.data_name = [repmat(1,1,numel(bcc_dat.data_name)-1) ; 2:numel(bcc_dat.data_name)];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'trials';
    bcc_dat = ft_func(@ft_appenddata,cfg,bcc_dat);
    
    cfg = [];
    cfg.rmfield   = 'yes';
    cfg.data_name =  2:numel(bcc_dat.data_name);
    bcc_dat = ft_func([],cfg,bcc_dat);
end
        
%% Filter Data for LFP & Baseline
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 15;
cfg.new_suffix = 'lfp';
bcc_dat        = ft_func(@ft_preprocessing,cfg,bcc_dat);

%% Filter Data for HGP & Baseline
cfg            = [];
cfg.data_name  = 1;
cfg.hilbert    = 'amp';
cfg.freq_band  = {[70 170]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
bcc_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,bcc_dat);

cfg             = [];
cfg.data_name   = 1*2+1:1*3;
HGP_smooth_msec = 0.065;
w               = round(HGP_smooth_msec* bcc_dat.(bcc_dat.data_name{1}).fsample); 
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,w/2);
cfg.window      = gauss_y/sum(gauss_y);
bcc_dat         = ft_func(@ft_window_data,cfg,bcc_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1;
cfg.rmfield   = 'yes';
bcc_dat       = ft_func([],cfg,bcc_dat);

%% Baseline Data & Remove Padding
cfg = [];
cfg.latency = [-0.8 1.6];
bcc_dat     = ft_func(@ft_selectdata,cfg,bcc_dat);

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.4 0];
if strcmpi(subj,'NY226_SA_SZ'); cfg.baselinewindow = [-0.125 0]; end
bcc_dat            = ft_func(@ft_preprocessing,cfg,bcc_dat);

%% Automatic rejection & Apply Rejections
cfg = [];
cfg.measures  = {'time' 'variance'};
cfg.thresh    = [0.98 0.98];
cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
cfg.prefix    = bcc_dat.data_name;
cfg.specific  = {'prefix';1:numel(bcc_dat.data_name)};
cfg.pad       = 0.4;
cfg.plot      = plt_rej;
bcc_dat       = ft_func(@auto_rej,cfg,bcc_dat);

cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
bcc_dat     = ft_func(@ft_apply_rejection,cfg,bcc_dat);

%% Remove Padding
cfg = [];
cfg.latency = [-0.4 1.2];
bcc_dat     = ft_func(@ft_selectdata,cfg,bcc_dat);

%% Event Upkeep
% Task Specific Events
cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 3};
cfg.new_events  = {101 102};
cfg.crt_alt_eve = 'vis_tot_nse';
bcc_dat = ft_func(@ft_redefine_events,cfg,bcc_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 4};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_tot_nse';
bcc_dat = ft_func(@ft_redefine_events,cfg,bcc_dat);

% Phoneme Events


%% Fix Labels
if ~exist([fcfg.clr_fld '/' 'ele_idn' '/' subj '/' bcc_dat.data_name{1} '.csv'])
    chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');
    bcc_dat.(bcc_dat.data_name{1}).cfg.alt_lab.label = ft_correct_channel2(chn_loc,bcc_dat.(bcc_dat.data_name{1}).label);
    for iA = 2:numel(bcc_dat.data_name); bcc_dat.(bcc_dat.data_name{iA}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_lab.label; end
else
    for iA = 1:numel(bcc_dat.data_name); 
        ele_idn = mmil_readtext([fcfg.clr_fld '/' 'ele_idn' '/' subj '/' bcc_dat.data_name{iA} '.csv']);
        bcc_dat.(bcc_dat.data_name{iA}).cfg.alt_lab.label = ele_idn(:,1); 
    end
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end
