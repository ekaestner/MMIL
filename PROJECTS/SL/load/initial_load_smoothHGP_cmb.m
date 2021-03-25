%% SL Analysis Script
clear; clc;

sbj_num = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting on Subject %s \n\n\n',subj)

indir       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/indir/' subj]);  
cln_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cln_dir/' subj]); 
end_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/end_dir/' subj]); 
tsk         = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/tsk/' subj]);     
rmv_chn     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/rmv_chn/' subj]),'uni',0)';
bse_frq     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/bse_frq/' subj]),'uni',0)';
trl_fun = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/trialfun/' subj]);

chn_nse_grp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' subj]);
tmp = strfind(chn_nse_grp,';');
if any(~isempty(tmp{:}))
    for i = 1:numel(chn_nse_grp)
        chn_nse_grp_dat_nme(i) = str2num(chn_nse_grp{i}(tmp{1}+1:end))';
        chn_nse_grp{i} = str2num(chn_nse_grp{i}(1:tmp{1}-1))'; 
    end
else
    chn_nse_grp = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/cmn_nse/' subj]),'uni',0)';
end

% In Script Switches 
plt_spc = 0;

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

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sll_dat      = ft_func(@ft_preprocessing,cfg,sll_dat);
end

%% Filter Data for HGP & Baseline
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[70 170]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
sll_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sll_dat);

cfg             = [];
cfg.data_name   = numel(inpath)+1:numel(inpath)*2;
HGP_smooth_msec = 0.065;
w               = round(HGP_smooth_msec* sll_dat.NY439_SL_Day3_Block1_1_Clin1.fsample); 
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,w/2);
cfg.window      = gauss_y/sum(gauss_y);
sll_dat       = ft_func(@ft_window_data,cfg,sll_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(inpath);
cfg.rmfield   = 'yes';
sll_dat       = ft_func([],cfg,sll_dat);

%% Epoch Data
if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/' '/trialfun_output/' subj '.mat'])
    cfg = [];
    cfg.trialfun = trl_fun{1};
    cfg.dataset  = inpath;
    cfg.specific = {'dataset' ; 1:numel(inpath)};
    cfg.minduration = 0.500;
    cfg.pre         = 1.5;
    cfg.pos         = 3;
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
    
    tmp = input('Satisfied with defined trials?: '); clear tmp
    
    save(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
else
    load(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
end

cfg = [];
cfg.specific  = {'trl';[1:numel(inpath) 1:numel(inpath) 1:numel(inpath)]};
cfg.trl       = trl_hld;
sll_dat = ft_func(@ft_redefinetrial,cfg,sll_dat);

%% Combine Data
% combine trials
num_trl       = numel(inpath);
cfg           = [];
cfg.data_name = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2; ...
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
cfg.data_name = [1 ; 2];
cfg.data_new  = 'yes';
cfg.methapp   = 'channels';
sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = [2:2:numel(sll_dat.data_name)];
sll_dat = ft_func([],cfg,sll_dat);

%% Baseline Data
cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.40 0];
sll_dat            = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Automatic rejection
cfg = [];
cfg.measures  = {'time' 'variance'};
cfg.thresh    = [0.98 0.98];
cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
cfg.prefix    = sll_dat.data_name;
cfg.specific  = {'prefix';1};
cfg.pad       = 1;
cfg.plot      = 0;
sll_dat       = ft_func(@auto_rej,cfg,sll_dat);

%% Apply Rejections
cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sll_dat     = ft_func(@ft_apply_rejection,cfg,sll_dat);

%% Remove Padding
cfg = [];
cfg.latency = [-0.5 2];
sll_dat     = ft_func(@ft_selectdata,cfg,sll_dat);

%% Keeping track of events
for iST = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.trialinfo = sll_dat.(sll_dat.data_name{iST}).trialinfo;
end

% Control Events
cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 3};
cfg.new_events  = {101 102};
cfg.crt_alt_eve = 'vis_tot_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 4};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_tot_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

% Phonemes
v_con = []; v_vow = []; a_con = []; a_vow = [];
for iTS = 1:2:numel(inpath) % - EJK need to come back and fix, quick hack
       
        blk_ind = regexpi(inpath{iTS},'\d_\d');
        blk = num2str(inpath{iTS}(blk_ind));
        run = num2str(inpath{iTS}(blk_ind+2));          
        lst = mmil_readtext(['/home/ekaestne/Desktop/N439_SL_Behav/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv']);

        v_con = [v_con ; cell2mat(lst(:,4))];
        v_vow = [v_vow ; cell2mat(lst(:,5))];
        a_con = [a_con ; cell2mat(lst(:,6))];
        a_vow = [a_vow ; cell2mat(lst(:,7))];
        
end

for iST = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.v_con = v_con;
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.v_vow = v_vow;
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.a_con = a_con;
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.a_vow = a_vow;
end

%% Label Upkeep
for iA = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.label = sll_dat.(sll_dat.data_name{iA}).label;
end

% Repair Macro
tru_mac = [16:-1:1 32:-1:17 48:-1:33 64:-1:49];
for iA = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label = sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.label;
    mac_ind = find(~cellfun(@isempty,regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label,'G_')));
    for iC = 1:numel(mac_ind)
        num_ind = regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)},'G_')+2;
        num_ind(2) = num_ind+1;
        rpl_num = tru_mac(str2double(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind)));
        if rpl_num>9
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind) = num2str(rpl_num);
        else
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind) = ['0' num2str(rpl_num)];
        end
    end
end

tru_mac = [6:-1:1];
for iA = 1:numel(sll_dat.data_name)
    mac_ind = find(~cellfun(@isempty,regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label,'LPT_')));
    for iC = 1:numel(mac_ind)
        num_ind = regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)},'LPT_')+4;
        num_ind(2) = num_ind+1;
        rpl_num = tru_mac(str2double(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind)));
        if rpl_num>9
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind) = num2str(rpl_num);
        else
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label{mac_ind(iC)}(num_ind) = ['0' num2str(rpl_num)];
        end
    end
end

% Repair Meso
tru_mes = [16:-1:1 32:-1:17 64:-1:49 48:-1:33];
for iA = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label = sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_label;
    mes_ind = find(~cellfun(@isempty,regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label,'Gr_')));
    for iC = 1:numel(mes_ind)
        num_ind = regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)},'Gr_')+3;
        num_ind(2) = num_ind+1;
        rpl_num = tru_mes(str2double(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind)));
        if rpl_num>9
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind) = num2str(rpl_num);
        else
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind) = ['0' num2str(rpl_num)];
        end
    end
end

tru_mes = [16:-1:1];
for iA = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label = sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label;
    mes_ind = find(~cellfun(@isempty,regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label,'LPTr_')));
    for iC = 1:numel(mes_ind)
        num_ind = regexpi(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)},'LPTr_')+5;
        num_ind(2) = num_ind+1;
        rpl_num = tru_mes(str2double(sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind)));
        if rpl_num>9
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind) = num2str(rpl_num);
        else
            sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label{mes_ind(iC)}(num_ind) = ['0' num2str(rpl_num)];
        end
    end
end

% Get labeling correct
chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
sll_dat.(sll_dat.data_name{1}).cfg.alt_lab.repair_macro_ejk1st_meso_label = ft_correct_channel(chn_loc,sll_dat.(sll_dat.data_name{1}).cfg.alt_lab.repair_macro_ejk1st_meso_label);
for iA = 2:numel(sll_dat.data_name); sll_dat.(sll_dat.data_name{iA}).cfg.alt_lab.repair_macro_ejk1st_meso_label = sll_dat.(sll_dat.data_name{1}).cfg.alt_lab.repair_macro_ejk1st_meso_label; end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename = [outpath '/' subj '_smoothed_hgp.mat'];
ft_func([],cfg,sll_dat);
