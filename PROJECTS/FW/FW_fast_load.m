function FW_fast_load(fcfg)

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

%% Remove identified Unimportant & Bad Channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

fwv_dat = mmil_remove_channels(cfg,fwv_dat);

%% Remove noise through removing common noise
chn_nse_grp = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');
nse_val         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'nse_val'); nse_val = nse_val{1};

% Make temporary ECOG/DEPTH Split
for iD = 1:numel(fwv_dat.data_name); fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.label = fwv_dat.(fwv_dat.data_name{iD}).label; end
scfg = [];
scfg.dat_nme = fwv_dat.data_name;
scfg.clr_fld = fcfg.clr_fld;
scfg.sbj_nme = fcfg.sbj_nme;
scfg.alt_lab = 'label';
scfg.specific = {'dat_nme' ; 1:numel(fwv_dat.data_name)};
scfg.add_nme  = '_tmp';
ft_func(@mmil_create_depth2,scfg,fwv_dat);

if nse_val==1   
    
    if strcmpi(chn_nse_grp{1},'split'); chn_nse_grp = repmat(chn_nse_grp,1,numel(fwv_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(fwv_dat.data_name); pre_fix{iD}{1} = strcat(fwv_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(fwv_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(fwv_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(fwv_dat.data_name{iD},'noisy_depth'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(fwv_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.sbj_nme     = fcfg.sbj_nme;
    cfg.clr_fld     = fcfg.clr_fld;
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = fwv_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(fwv_dat.data_name) 1:numel(fwv_dat.data_name) 1:numel(fwv_dat.data_name)};
    cfg.out_dir     = [fcfg.out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    fwv_dat = ft_func(@ft_remove_common_noise2,cfg,fwv_dat);    
    
    if nse_val == 1
        cfg.clr_fld = fcfg.clr_fld;
        cfg.sbj_nme = fcfg.sbj_nme;
        cfg.sub_fld_ind = num2cell(1:numel(fwv_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(fwv_dat.data_name) 1:numel(fwv_dat.data_name)};
        fwv_dat = ft_func(@mmil_update_cmn_nse2,cfg,fwv_dat);
    end
    
end

%% Combine Clinical/Day if necessary
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;

    fwv_dat = mmil_combine_data2(cfg,fwv_dat);
end

%% Fix Labels
cfg = [];
cfg.ovr_wrt = 0;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;
cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');

fwv_dat = mmil_fix_labels(cfg,fwv_dat);

%% Setup events
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.clr_fld = fcfg.clr_fld;
cfg.inpath  = inpath;
cfg.end_dir = end_dir;
cfg.tsk     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk');
fw_sve_eve(cfg)

%%
end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]);
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     

if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes.csv']);
elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
    eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes.csv']);
end

eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events.csv']);
if strcmpi(tsk,'fw_1')
    idn_col = 2;
    eve_col = 4;
else
    idn_col = 4;
    eve_col = 3;
end

%% Check data
wrg_num = 1;
wrg_mtc = 1;

while wrg_num || wrg_mtc
    
    if ~exist([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],{});
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],{});
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],{})
    else
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        end
        eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
        
    if size(eve_cde,1) ~= size(eve_fwv,1)
        disp([fcfg.sbj_nme ': number of trials do no match' ])
        pause()
        pause(5)       
    elseif ~all(cell2mat(eve_cde(:,1)) == cell2mat(eve_fwv(:,4)))
        wrg_num = 0;
        
        disp([fcfg.sbj_nme ': events do no match' ])
        find(cell2mat(eve_cde(:,1)) ~= cell2mat(eve_fwv(:,4)))
        pause()
        pause(5)
    else
        wrg_num = 0;
        wrg_mtc = 0;
    end

    
    try isempty(mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'])); catch; emp_txt = 1; end
    if ~wrg_num && ~wrg_mtc && isvar('emp_txt') && emp_txt
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv'],eve_cde);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv'],eve_cde);
        end
        cell2csv([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv'],eve_fwv)
    else
        if strcmpi(end_dir,'eeg') | strcmpi(end_dir,'edf')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_eeg_codes_corrected.csv']);
        elseif strcmpi(end_dir,'mat') | strcmpi(end_dir,'.mat')
            eve_cde = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_dio_codes_corrected.csv']);
        end
        eve_fwv = mmil_readtext([fcfg.clr_fld '/' 'task' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_events_corrected.csv']);
    end
    
end

end