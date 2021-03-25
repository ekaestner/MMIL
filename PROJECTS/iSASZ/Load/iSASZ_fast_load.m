function iSASZ_fast_load(fcfg)

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

cfg = [];
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.cln_dir = cln_dir;
cfg.end_dir = end_dir;
cfg.tsk     = tsk;
[inpath,sem_dat,trl] = mmil_find_files(cfg);

%% Load Visual & Auditory data
epc_tme      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

for iL = 1:numel(sem_dat.data_name)
    
    if ~isempty(string_find(end_dir(iL),{'mat'}))
        ttt = load(inpath{iL});
        sem_dat.(sem_dat.data_name{iL}) = ttt.epoch_data;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
    elseif ~isempty(string_find(end_dir(iL),{'eeg'}))
        ttt = ts_load_data(inpath{iL},'accept_all_flag',1);
        sem_dat.(sem_dat.data_name{iL}) = ttt;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
    elseif ~isempty(string_find(end_dir(iL),{'set'}))
        ttt = ts_load_data(inpath{iL},'accept_all_flag',1);
        sem_dat.(sem_dat.data_name{iL}) = ttt;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
        
    elseif ~isempty(string_find(end_dir(iL),{'edf'}))
        
        cfg            = [];
        cfg.specific   = {'dataset';1:numel(inpath)};
        cfg.data_new   = 'yes';
        cfg.continuous = 'yes';
        cfg.dataset    = inpath;
        sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);
        
        if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
            
            trialfun       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);  
            
            cfg = [];
            if numel(trialfun) == 1; cfg.trialfun = trialfun{1}; else cfg.trialfun = trialfun; end
            cfg.dataset  = inpath;
            if numel(trialfun) == 1; cfg.specific = {'dataset' ; 1:numel(inpath)}; else cfg.specific = {'dataset' 'trialfun' ; 1:numel(inpath) 1:numel(trialfun)}; end
            cfg.minduration = 0.500;
            cfg.pre         = epc_tme(1)-1;
            cfg.pos         = epc_tme(2)+1;
            cfg.evt         = 1:8;
            trl = ft_func(@ft_definetrial,cfg,trl);
            
            if any(strcmpi(subj,{'NY496_SA_SZ' 'NY503_SA_SZ' 'NY511_SA_SZ'}))
                trl = sasz_fix_trl(subj,trl);
            end
            
            trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end
            
            if any(strcmpi(subj,{'NY442_SA_SZ'})); trl_hld{5}(341:end,:) = []; end
            
            for i = 1:numel(tsk); if strcmpi(tsk{i},'SA'); trl_hld{i}(:,4) = trl_hld{i}(:,4)+10; end; end
            
            save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
        else
            load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
            if any(strcmpi(subj,{'NY511_SA_SZ'})); 
                trl_hld{3}(145:end,:) = [];
            end            
        end
        
        cfg = [];
        cfg.specific  = {'trl';1:numel(inpath)};
        cfg.trl       = trl_hld;
        sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
        
    end
    
end

%% Misc
for iL = 1:numel(sem_dat.data_name)
    if strcmpi(tsk{iL},'SA') && all(sem_dat.(sem_dat.data_name{iL}).trialinfo < 10)
        sem_dat.(sem_dat.data_name{iL}).trialinfo = sem_dat.(sem_dat.data_name{iL}).trialinfo + 10;
    end
end

if strcmpi(subj,'NY008_SA_SZ')
    sem_dat.(sem_dat.data_name{2}).label = sem_dat.(sem_dat.data_name{1}).label;
elseif strcmpi(subj,'NY011_SA_SZ')
    sem_dat.(sem_dat.data_name{2}).label = sem_dat.(sem_dat.data_name{1}).label;
    
    cfg = [];
    cfg.latency = [-1.4 2.4];
    sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
elseif strcmpi(subj,'NY068_SA')
    sem_dat.NY68_SAo.trialinfo = sem_dat.NY68_SAo.trialinfo';
    sem_dat.NY68_SAo.trial = sem_dat.NY68_SAo.trial';
elseif strcmpi(subj,'NY017_SA_SZ')
    sem_dat.NY17_SA.cfg.trl = [sem_dat.NY17_SA.sampleinfo sem_dat.NY17_SA.sampleinfo(:,2)];
    sem_dat.NY17_SZ.cfg.trl = [sem_dat.NY17_SZ.sampleinfo sem_dat.NY17_SZ.sampleinfo(:,2)];
elseif strcmpi(subj,'NY226_SA_SZ')
    cfg = [];
    cfg.latency = [-0.2 0.8];
    sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
    
    sem_dat.(sem_dat.data_name{1}).fsample = upper(sem_dat.(sem_dat.data_name{1}).fsample);
    
    sem_dat.(sem_dat.data_name{1}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{1}).time));
    sem_dat.(sem_dat.data_name{2}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{2}).time));
    sem_dat.(sem_dat.data_name{3}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{3}).time));
    
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
    sem_dat = mmil_combine_data2(cfg,sem_dat);
end

%% Remove identified Unimportant & Bad Channels
% identified channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

sem_dat = mmil_remove_channels(cfg,sem_dat);

%% Remove noise through removing common noise
chn_nse_grp = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');
nse_val         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'nse_val'); nse_val = nse_val{1};

% Make temporary ECOG/DEPTH Split
for iD = 1:numel(sem_dat.data_name); sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{iD}).label; end
scfg = [];
scfg.dat_nme = sem_dat.data_name;
scfg.clr_fld = fcfg.clr_fld;
scfg.sbj_nme = fcfg.sbj_nme;
scfg.alt_lab = 'label';
scfg.specific = {'dat_nme' ; 1:numel(sem_dat.data_name)};
scfg.add_nme  = '_tmp';
ft_func(@mmil_create_depth2,scfg,sem_dat);

if nse_val
    
    if strcmpi(chn_nse_grp{1},'split'); chn_nse_grp = repmat(chn_nse_grp,1,numel(sem_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(sem_dat.data_name); pre_fix{iD}{1} = strcat(sem_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(sem_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(sem_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(sem_dat.data_name{iD},'noisy_depth'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(sem_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.sbj_nme     = fcfg.sbj_nme;
    cfg.clr_fld     = fcfg.clr_fld;
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = sem_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name)};
    cfg.out_dir     = [fcfg.out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    sem_dat = ft_func(@ft_remove_common_noise2,cfg,sem_dat);   
    
    if nse_val == 1
        cfg.clr_fld = fcfg.clr_fld;
        cfg.sbj_nme = fcfg.sbj_nme;
        cfg.sub_fld_ind = num2cell(1:numel(sem_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@mmil_update_cmn_nse2,cfg,sem_dat);
    end
    
end

%% Setup events
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.clr_fld = fcfg.clr_fld;
cfg.inpath  = inpath;
cfg.end_dir = end_dir;
cfg.tsk     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk');
sasz_sve_eve(cfg,sem_dat)

%% Combine Clinical/Day if necessary
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;

    sem_dat = mmil_combine_data2(cfg,sem_dat);
end

%% Fix Labels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;
cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');

sem_dat = mmil_fix_labels(cfg,sem_dat);

% Load codes
if ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif ~isempty(string_find(end_dir,'mat'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
elseif ~isempty(string_find(end_dir,'edf'))
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes.csv']);
end

% Load real events
stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv']);
eve_col = 3;

%% Check data
wrg_num = 1;
wrg_mtc = 1;

while wrg_num || wrg_mtc
    
    if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])
        if  ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],{});
        elseif ~isempty(string_find(end_dir,'mat'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],{});
        elseif ~isempty(string_find(end_dir,'edf'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv'],{});
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],{})
    else
        if  ~isempty(string_find(end_dir,'eeg')) | ~isempty(string_find(end_dir,'set'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'mat'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'edf'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv']);
        end
        stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
    if size(eve_cde,1) ~= size(stm_cde,1)
        disp([fcfg.sbj_nme ': number of trials do no match' ])
        pause()
        pause(5)
    elseif ~all(cell2mat(eve_cde(:,1)) == cell2mat(stm_cde(:,eve_col)))
        wrg_num = 0;
        
        disp([fcfg.sbj_nme ': events do no match' ])
        find(cell2mat(eve_cde(:,1)) ~= cell2mat(stm_cde(:,eve_col)))
        pause()
        pause(5)
    else
        wrg_num = 0;
        wrg_mtc = 0;
    end
    
    
    try isempty(mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])); catch; emp_txt = 1; end
    if ~wrg_num && ~wrg_mtc && isvar('emp_txt') && emp_txt
        if ~isempty(string_find(end_dir,'eeg'))  | ~isempty(string_find(end_dir,'set'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],eve_cde);
        elseif ~isempty(string_find(end_dir,'mat'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],eve_cde);
            elseif ~isempty(string_find(end_dir,'edf'))
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv'],eve_cde);
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],stm_cde)
    else
        if ~isempty(string_find(end_dir,'eeg'))  | ~isempty(string_find(end_dir,'set'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif ~isempty(string_find(end_dir,'mat'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
            elseif ~isempty(string_find(end_dir,'edf'))
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_edf_codes_corrected.csv']);
        end
        stm_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
end

end