clear; clc;

sbj_num = 10;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

cfg = [];
cfg.load    = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sem_dat     = ft_func([],cfg);

cln_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';

try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
catch
    has_loc = 0;
end
if ~isvar('has_loc'); has_loc = 1; end

if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo > 9); has_aud = 1; else has_aud = 0; end
if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo < 9); has_vis = 1; else has_vis = 0; end

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
        
        all_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl;
        rep_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl([2 5]);
        sem_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl([3 6]);
        
        all_lab_pos = [1 2 3 6 7 8];
        rep_lab_pos = [1 6];
        sem_lab_pos = [1 6];
        
        if fle
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn(:,[2 5]),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn(:,[3 6]),1),'Uni',0),'Uni',0);
        else
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,[2 5]),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,[3 6]),1),'Uni',0),'Uni',0);
        end
        
        all_col     = {ft_stt_col(rgb('red')) rgb('red') rgb('green') ft_stt_col(rgb('blue')) rgb('blue') rgb('yellow')};
        rep_col     = {rgb('red') rgb('blue')};
        sem_col     = {rgb('green') rgb('yellow')};
        
    elseif has_aud
        
        all_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl;
        rep_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl(2);
        sem_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl(3);
        
        all_lab_pos = [1 2 3 6 7 8];
        rep_lab_pos = [1 6];
        sem_lab_pos = [1 6];
        
        if fle
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,2),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,3),1),'Uni',0),'Uni',0);
        else
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,2),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,3),1),'Uni',0),'Uni',0);
        end
        
        all_col     = {ft_stt_col(rgb('blue')) rgb('blue') rgb('yellow')};
        rep_col     = {rgb('blue')};
        sem_col     = {rgb('yellow')};
        
    elseif has_vis
        
        all_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl;
        rep_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl(2);
        sem_sel_lbl = sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.lbl(3);
        
        all_lab_pos = [1 2 3];
        rep_lab_pos = [1];
        sem_lab_pos = [1];
        
        if fle
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,2),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,3),1),'Uni',0),'Uni',0);
        else
            all_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn,1),'Uni',0),'Uni',0);
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,2),1),'Uni',0),'Uni',0);
            sem_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sem_dat.(sem_dat.data_name{iD}).cfg.alt_sig_chn.chn(:,3),1),'Uni',0),'Uni',0);
        end
        
        all_col     = {ft_stt_col(rgb('red')) rgb('red') rgb('green')};
        rep_col     = {rgb('red')};
        sem_col     = {rgb('green')};
        
    end
    
    %% Electrode all condition locations
    if has_loc
        cfg = [];
        
        cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location' '/' subj]);
        cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/pial' '/' subj]);
        cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi' '/' subj]);
        
        cfg.sel_lbl   = all_sel_lbl;
        cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi/' subj]);
        cfg.lab_pos   = all_lab_pos;
        cfg.sel_ele   = all_sel_ele;
        cfg.col       = all_col;
        
        cfg.sve_img   = 'jpg';
        cfg.sve_loc   = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'initial_location' '/' sem_dat.data_name{iD}];
        cfg.sve_pre   = [subj '_all'];
        
        mmil_ieeg_sensor_location_plot_v1(cfg);
        
        
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
        cfg.sve_loc   = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'initial_location' '/' sem_dat.data_name{iD}];
        cfg.sve_pre   = [subj '_rep'];
        
        mmil_ieeg_sensor_location_plot_v1(cfg);
        
        %% Electrode semantic locations
        cfg = [];
        
        cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location' '/' subj]);
        cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/pial' '/' subj]);
        cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi' '/' subj]);
        
        cfg.sel_lbl   = sem_sel_lbl;
        cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi/' subj]);
        cfg.lab_pos   = sem_lab_pos;
        cfg.sel_ele   = sem_sel_ele;
        cfg.col       = sem_col;
        
        cfg.sve_img   = 'jpg';
        cfg.sve_loc   = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'initial_location' '/' sem_dat.data_name{iD}];
        cfg.sve_pre   = [subj '_sem'];
        
        mmil_ieeg_sensor_location_plot_v1(cfg);
        
    end
    
end