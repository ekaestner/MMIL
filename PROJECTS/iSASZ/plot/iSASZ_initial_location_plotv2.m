clear; clc;

sbj_num = 2;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

cfg = [];
cfg.load    = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sem_dat     = ft_func([],cfg);

cln_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';

try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
catch
    has_loc = 0;
end
if ~isvar('has_loc'); has_loc = 1; end

if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo > 9); has_aud = 1; else has_aud = 0; end
if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo < 9); has_vis = 1; else has_vis = 0; end

fprintf('Starting work on %s \n',subj)

if has_loc; for iA = 2:numel(sem_dat.data_name); sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label; end; end

for iD = 1:numel(sem_dat.data_name)
    
    cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}];
    if exist([cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}],'file');
        fle = 1;
        sig_chn = mmil_readtext(cfg.sig_fle,[',']);
        tmp = find(strcmpi(sig_chn(1,:),'stars'));
        sig_chn = cell2mat(sig_chn(2:end,2:tmp-2));
    else
        fle = 0;
    end
   
    if has_aud && has_vis
        
        rep_sel_lbl = {'Visual Lexical Overall Activity' 'Visual Lexical Repetition Effect' 'Auditory Lexical Overall Activity' 'Auditory Lexical Repetition Effect'};
        
        rep_lab_pos = [1 2 6 7];
        
        if fle
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn(:,[1 3 7 9]),1),'Uni',0),'Uni',0);
        end
        
        rep_col     = {rgb('reddish grey') rgb('neon red') rgb('greyish blue') rgb('cyan')};
        
    elseif has_aud
        
        rep_sel_lbl = {'Auditory Lexical Overall Activity' 'Auditory Lexical Repetition Effect'};
        
        rep_lab_pos = [1 2];
        
        if fle
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn(:,[1 3]),1),'Uni',0),'Uni',0);
        end
        
        rep_col     = {rgb('greyish blue') rgb('cyan')};
        
    elseif has_vis
        
        rep_sel_lbl = {'Visual Lexical Overall Activity' 'Visual Lexical Repetition Effect'};
        
        rep_lab_pos = [1 2];
        
        if fle
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn(:,[1 3]),1),'Uni',0),'Uni',0);
        end
        
        rep_col     = {rgb('reddish grey') rgb('neon red')};
        
    end
    
    %% Electrode all condition locations
    if has_loc      
            
        %% Electrode repetition locations
        cfg = [];
        
        cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location' '/' subj]);
        cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/pial' '/' subj]);
        cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi' '/' subj]);
        
        cfg.sel_lbl   = rep_sel_lbl;
        cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi/' subj]);
        cfg.lab_pos   = rep_lab_pos;
        cfg.sel_ele   = rep_sel_ele;
        cfg.col       = rep_col;
        
        cfg.sve_img   = 'jpg';
        cfg.sve_loc   = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' 'ThreeStat' '/' subj '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.sve_pre   = [subj '_repetition_and_overall'];
        
        mmil_ieeg_sensor_location_plot_v1(cfg);
        
    end
    
end
