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
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';

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

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'before' '/' subj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
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

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' 'after' '/' subj '/' ];
    cfg.prefix = [sll_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sll_dat);
end

%% Filter Data for LFP & Baseline
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 20;
cfg.new_suffix = 'lfp';
sll_dat        = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Filter Data for HGP & Baseline
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[70 170]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
sll_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sll_dat);

cfg           = [];
cfg.data_name = [numel(inpath)*2+1:numel(inpath)*3];
cfg.lpfilter  = 'yes';
cfg.lpfreq    = 20;
sll_dat       = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Filter Data for Theta
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[3 7]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'tht';
sll_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sll_dat);

cfg           = [];
cfg.data_name = numel(inpath)*3+1:numel(inpath)*4;
cfg.lpfilter  = 'yes';
cfg.lpfreq    = 20;
sll_dat       = ft_func(@ft_preprocessing,cfg,sll_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(inpath);
cfg.rmfield   = 'yes';
sll_dat       = ft_func([],cfg,sll_dat);

%% Epoch Data
if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'])
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
        
        % - EJK add auditory specific timings here
        tmp_trl = [trl_hld{iTS}(:,1:2)+(round(0.450*sll_dat.(sll_dat.data_name{iTS}).fsample)) trl_hld{iTS}(:,3) trl_hld{iTS}(:,4)+10];
        trl_hld{iTS} = [trl_hld{iTS} ; tmp_trl];
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
cfg.data_name = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2 ...
                 ones(1,numel(1:2:num_trl)-1)+num_trl ones(1,numel(2:2:num_trl)-1)*2+num_trl ...
                 ones(1,numel(1:2:num_trl)-1)+num_trl*2 ones(1,numel(2:2:num_trl)-1)*2+num_trl*2; ...
                3:2:num_trl 4:2:num_trl ...
                3+num_trl:2:num_trl*2 4+num_trl:2:num_trl*2 ...
                3+num_trl*2:2:num_trl*3 4+num_trl*2:2:num_trl*3];
cfg.data_new  = 'yes';
cfg.methapp   = 'trials';
sll_dat = ft_func(@ft_appenddata,cfg,sll_dat);

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = [3:2:num_trl 4:2:num_trl ...
                3+num_trl:2:num_trl*2 4+num_trl:2:num_trl*2 ...
                3+num_trl*2:2:num_trl*3 4+num_trl*2:2:num_trl*3];
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
cfg.specific  = {'prefix';[1 2 3]};
cfg.pad       = 0.4;
cfg.plot      = 1;
sll_dat       = ft_func(@auto_rej,cfg,sll_dat);

%% Apply Rejections
cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sll_dat     = ft_func(@ft_apply_rejection,cfg,sll_dat);

%% Remove Padding
cfg = [];
cfg.latency = [-0.5 1.5];
sll_dat     = ft_func(@ft_selectdata,cfg,sll_dat);

% %% Add in SNR measures
% cfg = [];
% cfg.nse_tme = [-0.5 0];
% cfg.sig_tme = [0 0.8];
% cfg.eve     = {[1] [2] [3]};
% cfg.snr_lab = {'vis' 'aud' ''};
% sll_dat    = ft_func(@ft_SNR,cfg,sll_dat);

%% Keeping track of events
for iST = 1:numel(sll_dat.data_name)
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.trialinfo = sll_dat.(sll_dat.data_name{iST}).trialinfo;
end

% Control Events
cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 3};
cfg.new_events  = {101 102};
cfg.crt_alt_eve = 'vis_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[11 12] 14};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2 4] 3};
cfg.new_events  = {101 102};
cfg.crt_alt_eve = 'vis_tot_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[11 12 13] 14};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_tot_nse';
sll_dat = ft_func(@ft_redefine_events,cfg,sll_dat);

% Phonemes
v_con = []; v_vow = []; a_con = []; a_vow = [];
for iTS = 1:numel(inpath) % - EJK need to come back and fix, quick hack
       
        blk_ind = regexpi(inpath{iTS},'\d_\d');
        blk = num2str(inpath{iTS}(blk_ind));
        run = num2str(inpath{iTS}(blk_ind+2));          
        lst = mmil_readtext(['/home/ekaestne/Desktop/N439_SL_Behav/lists' '/' 'SL_Vis1st_' blk '_' run '_1.csv']);

        v_con = [v_con ; cell2mat(lst(:,4)) ; cell2mat(lst(:,4))];
        v_vow = [v_vow ; cell2mat(lst(:,5)) ; cell2mat(lst(:,5))];
        a_con = [a_con ; cell2mat(lst(:,6)) ; cell2mat(lst(:,6))];
        a_vow = [a_vow ; cell2mat(lst(:,7)) ; cell2mat(lst(:,7))];
        
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

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sll_dat);
