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
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_freq/';

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

tmp_dat_exs = cellfun(@(x) exist([outpath '/' 'temp_freq' '/' x '_fourier.mat'],'file'),trl.data_name);
if ~all(tmp_dat_exs)
    
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
    
    %% Create Wavelet Dataset
    % parameters for wavelet function
    foi =  [2:12,14:2:24,25,30:5:55,70:10:170];    %frequency of interest
    sf  =   [1,1,repmat(2,1,9),repmat(3,1,6),repmat(5,1,7),repmat(10,1,11)]; %specific frequency
    gwidth = ones(size(foi))*pi; %wavelet
    
    % find relevant channels
    
    tic
    if ~exist([outpath '/' 'temp_freq' '/']); mkdir([outpath '/' 'temp_freq' '/']); end
    for iDT = 1:numel(sll_dat.data_name)
        
        % Calculate
        tmp_frq.data_name{1}             = sll_dat.data_name{iDT};
        tmp_frq.(tmp_frq.data_name{1})   = sll_dat.(sll_dat.data_name{iDT});
        
        cfg = [];
        cfg.data_name  = 1;
        cfg.method     = 'wavelet';
        cfg.output     = 'fourier';
        cfg.foi        = foi;	      % frequency of interest
        cfg.sf         = sf;        % specific frequency
        cfg.width      = foi./sf;
        cfg.gwidth     = gwidth;
        cfg.toi        = tmp_frq.(tmp_frq.data_name{1}).time{1}; %cellfun(@(x) sll_dat.(x).time{1},sll_dat.data_name,'uni',0);
        cfg.keeptrials = 'yes';
        %     cfg.specific   = {'toi'; 1:numel(sll_dat.data_name)};
        tmp_frq        = ft_func(@ft_freqanalysis,cfg,tmp_frq);
        
        tmp_frq.(tmp_frq.data_name{1}).fsample = sll_dat.(sll_dat.data_name{iDT}).fsample;
        
        % Trim Timing
        cfg = [];
        cfg.data_name = 1;
        cfg.latency   = [-0.5 1.5];
        tmp_frq       = ft_func(@ft_selectdata,cfg,tmp_frq);
        
        % Trim Channels
        
        % Create Power
        tmp_frq_dat = tmp_frq.(tmp_frq.data_name{1}).fourierspctrm;
        save([outpath '/' 'temp_freq' '/' tmp_frq.data_name{1} '_fourier.mat'],'tmp_frq_dat','-v7.3')
        clear tmp_frq_dat
       
        tmp_frq.(tmp_frq.data_name{1}) = rmfield(tmp_frq.(tmp_frq.data_name{1}),'fourierspctrm');
        save([outpath '/' 'temp_freq' '/' tmp_frq.data_name{1} '_struct.mat'],'tmp_frq')
        
        tmp_sll.data_name{iDT}         = tmp_frq.data_name{1};
        tmp_sll.(tmp_frq.data_name{1}) = tmp_frq.(tmp_frq.data_name{1});
        
        clear tmp_frq
        
    end
    toc
    
end

%% Combine Power Data
% replace trial with frequency power information
tic
% load & combine trials
num_trl       = numel(inpath);
dat_cmb_ind = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2; ...
    3:2:num_trl 4:2:num_trl];

org_lod = unique(dat_cmb_ind(1,:));
cnt_lod = 1;
for iOR = 1:numel(org_lod)
    
    load([outpath '/' 'temp_freq' '/' sll_dat.data_name{org_lod(iOR)} '_struct.mat']);
    load([outpath '/' 'temp_freq' '/' sll_dat.data_name{org_lod(iOR)} '_fourier.mat']);
    
    sll_dat.(sll_dat.data_name{iOR})       = tmp_frq.(sll_dat.data_name{iOR});
    sll_dat.(sll_dat.data_name{iOR}).power = abs(tmp_frq_dat);
    
    clear tmp_frq tmp_frq_dat
    
    add_lod = unique(dat_cmb_ind(2,dat_cmb_ind(1,:)==org_lod(iOR)));
    for iAD = 1:numel(add_lod);
        
        load([outpath '/' 'temp_freq' '/' sll_dat.data_name{add_lod(iAD)} '_struct.mat']);
        load([outpath '/' 'temp_freq' '/' sll_dat.data_name{add_lod(iAD)} '_fourier.mat']);
        
        sll_dat.(sll_dat.data_name{add_lod(iAD)})       = tmp_frq.(sll_dat.data_name{add_lod(iAD)});
        sll_dat.(sll_dat.data_name{add_lod(iAD)}).power = abs(tmp_frq_dat);
        
        clear tmp_frq tmp_frq_dat
        
        cfg           = [];
        cfg.parameter = 'power';
        cfg.data_name = [dat_cmb_ind(1,cnt_lod) ; dat_cmb_ind(2,cnt_lod)];
        cfg.data_new  = 'yes';
        cfg.methapp   = 'trials';
        sll_dat = ft_func(@ft_appendfreq,cfg,sll_dat);
        
        sll_dat.(sll_dat.data_name{add_lod(iAD)}) = rmfield(sll_dat.(sll_dat.data_name{add_lod(iAD)}),'power');
        
        cnt_lod = cnt_lod + 1;
        
    end
    
end
toc

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

% combine Clin1 & Clin2
dat_cmb_ind = [1:2:numel(sll_dat.data_name) ; ...
    2:2:numel(sll_dat.data_name)];

for iCM = 1:numel(dat_cmb_ind(2,:)) %
    
    cfg           = [];
    cfg.parameter = 'power';
    cfg.data_name = [dat_cmb_ind(1,iCM) ; dat_cmb_ind(2,iCM)];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'channels';
    sll_dat = ft_func(@ft_appendfreq,cfg,sll_dat);
    
    sll_dat.(sll_dat.data_name{dat_cmb_ind(2,iCM)}) = rmfield(sll_dat.(sll_dat.data_name{dat_cmb_ind(2,iCM)}),'power');
    
end

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

%% Automatic rejection
if ~exist([outpath '/' 'data_information' '/' subj '/' 'badtrl_hld.mat'],'file')
    cfg = [];
    cfg.measures  = {'time' 'variance'};
    cfg.frq_rng   = {[1 12] [14 24] [30 55] [70 170]};
    cfg.thresh    = [0.99 0.99];
    cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
    cfg.prefix    = sll_dat.data_name{1};
    cfg.pad       = 0;
    cfg.plot      = 1;
    sll_dat       = ft_func(@auto_rej,cfg,sll_dat);
    
    badtrl_hld      = sll_dat.(sll_dat.data_name{1}).cfg.badtrl;
    badtrl_type_hld = sll_dat.(sll_dat.data_name{1}).cfg.badtrl_type;
    save([outpath '/' 'data_information' '/' subj '/' 'badtrl_hld.mat'],'badtrl_hld','badtrl_type_hld');
else
    load([outpath '/' 'data_information' '/' subj '/' 'badtrl_hld.mat']);
    sll_dat.(sll_dat.data_name{1}).cfg.badtrl      = badtrl_hld;
    sll_dat.(sll_dat.data_name{1}).cfg.badtrl_type = badtrl_type_hld;
end

%% Apply Rejections
cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sll_dat     = ft_func(@ft_apply_rejection,cfg,sll_dat);

%% Keeping track of events
if ~exist([outpath '/' 'data_information' '/' subj '/' 'alt_eve_hld.mat'],'file')
    
    eve_hld = [];
    for iE = 1:2:numel(inpath);
        eve_hld = [eve_hld ; tmp_sll.(tmp_sll.data_name{iE}).trialinfo];
    end
    sll_dat.(sll_dat.data_name{1}).trialinfo = eve_hld;
    
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
    cfg.old_events  = {[1 2] 4};
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
    cfg.old_events  = {[1 2 3] 4};
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
    
    alt_eve_hld = sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve;
    save([outpath '/' 'data_information' '/' subj '/' 'alt_eve_hld.mat'],'alt_eve_hld');
    
else
    load([outpath '/' 'data_information' '/' subj '/' 'alt_eve_hld.mat']);
    sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve = alt_eve_hld;
end

%% Label Upkeep
if ~exist([outpath '/' 'data_information' '/' subj '/' 'alt_lab_hld.mat'],'file')
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
    
    alt_lab_hld = sll_dat.(sll_dat.data_name{iST}).cfg.alt_lab;
    save([outpath '/' 'data_information' '/' subj '/' 'alt_lab_hld.mat'],'alt_lab_hld');
else
    load([outpath '/' 'data_information' '/' subj '/' 'alt_lab_hld.mat'])
    sll_dat.(sll_dat.data_name{1}).cfg.alt_lab = alt_lab_hld;
end

%% Save Power Data
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename = [outpath '/' subj '_overall_power_data.mat'];
ft_func([],cfg,sll_dat);

clear sll_dat

%% Combine Phase Data
% replace trial with frequency phase information
sll_dat = tmp_sll;

tic
% load & combine trials
num_trl       = numel(inpath);
dat_cmb_ind = [ones(1,numel(1:2:num_trl)-1) ones(1,numel(2:2:num_trl)-1)*2; ...
    3:2:num_trl 4:2:num_trl];

org_lod = unique(dat_cmb_ind(1,:));
cnt_lod = 1;
for iOR = 1:numel(org_lod)
    
    load([outpath '/' 'temp_freq' '/' sll_dat.data_name{org_lod(iOR)} '_struct.mat']);
    load([outpath '/' 'temp_freq' '/' sll_dat.data_name{org_lod(iOR)} '_fourier.mat']);
    
    sll_dat.(sll_dat.data_name{iOR})       = tmp_frq.(sll_dat.data_name{iOR});
    sll_dat.(sll_dat.data_name{iOR}).phase = angle(tmp_frq_dat);
    
    clear tmp_frq tmp_frq_dat
    
    add_lod = unique(dat_cmb_ind(2,dat_cmb_ind(1,:)==org_lod(iOR)));
    for iAD = 1:numel(add_lod);
        
        load([outpath '/' 'temp_freq' '/' sll_dat.data_name{add_lod(iAD)} '_struct.mat']);
        load([outpath '/' 'temp_freq' '/' sll_dat.data_name{add_lod(iAD)} '_fourier.mat']);
        
        sll_dat.(sll_dat.data_name{add_lod(iAD)})       = tmp_frq.(sll_dat.data_name{add_lod(iAD)});
        sll_dat.(sll_dat.data_name{add_lod(iAD)}).phase = angle(tmp_frq_dat);
        
        clear tmp_frq tmp_frq_dat
        
        cfg           = [];
        cfg.parameter = 'phase';
        cfg.data_name = [dat_cmb_ind(1,cnt_lod) ; dat_cmb_ind(2,cnt_lod)];
        cfg.data_new  = 'yes';
        cfg.methapp   = 'trials';
        sll_dat = ft_func(@ft_appendfreq,cfg,sll_dat);
        
        sll_dat.(sll_dat.data_name{add_lod(iAD)}) = rmfield(sll_dat.(sll_dat.data_name{add_lod(iAD)}),'phase');
        
        cnt_lod = cnt_lod + 1;
        
    end
    
end
toc

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

% combine Clin1 & Clin2
dat_cmb_ind = [1:2:numel(sll_dat.data_name) ; ...
    2:2:numel(sll_dat.data_name)];

for iCM = 1:numel(dat_cmb_ind(2,:)) %
    
    cfg           = [];
    cfg.parameter = 'phase';
    cfg.data_name = [dat_cmb_ind(1,iCM) ; dat_cmb_ind(2,iCM)];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'channels';
    sll_dat = ft_func(@ft_appendfreq,cfg,sll_dat);
    
    sll_dat.(sll_dat.data_name{dat_cmb_ind(2,iCM)}) = rmfield(sll_dat.(sll_dat.data_name{dat_cmb_ind(2,iCM)}),'phase');
    
end

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = dat_cmb_ind(2,:);
sll_dat = ft_func([],cfg,sll_dat);

% Add in rejection, trial, & label information
sll_dat.(sll_dat.data_name{1}).cfg.alt_eve      = alt_eve_hld;
sll_dat.(sll_dat.data_name{1}).cfg.alt_lab      = alt_lab_hld;
sll_dat.(sll_dat.data_name{1}).cfg.badtrl_type  = badtrl_type_hld;
sll_dat.(sll_dat.data_name{1}).cfg.badtrl       = badtrl_hld;

%% Apply Rejections
cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sll_dat     = ft_func(@ft_apply_rejection,cfg,sll_dat);

%% Save Phase Data
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename = [outpath '/' subj '_overall_phase_data.mat'];
ft_func([],cfg,sll_dat);

clear sll_dat


