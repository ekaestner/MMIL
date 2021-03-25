clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/lfp';
indir   = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

indir_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

ovr_ele_loc.lhs = {};
ovr_ele_loc.rhs = {};
ovr_ele_act = cell(2,4);

for iS = [1:3 5:16]
    
    sbj = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_lfp_data.mat'];
    sem_dat  = ft_func([],cfg);
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}]);
        
    try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' 'electrode_location' '/' sbj]);
    catch
        has_loc = 0;
    end
    if ~isvar('has_loc'); has_loc = 1; end
    
    cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}];
    if exist([clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}],'file');
        fle = 1;
        sig_chn = mmil_readtext(cfg.sig_fle,['\t']);
        tmp = find(strcmpi(sig_chn(1,:),'stars'));
        sig_chn = cell2mat(sig_chn(2:end,2:tmp-2));
    else
        fle = 0;
    end
    
    if has_loc
                
        cfg = [];
        cfg.load = 'yes';
        cfg.file = [indir_old '/' sbj '_overall_data.mat'];
        sem_dat_old  = ft_func([],cfg);
        sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label = sem_dat_old.(sem_dat_old.data_name{1}).cfg.alt_lab.label;
        clear sem_dat_old
        
        new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],['\t']);
        has_vis = any(strcmpi('visual_early',new_sig(1,:)));
        has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
        
        cur_ele_fle = mmil_readtext([clr_fld_old '/' 'electrode_location_mni' '/' sbj]);
        cur_ele_hms = mmil_readtext([clr_fld_old '/' 'hemi' '/' sbj]);
        
        for iF = 1:numel(cur_ele_fle)
            
            cur_ele_loc = mmil_readtext(cur_ele_fle{1},[' ']);
            ovr_ele_loc.(cur_ele_hms{iF}) = [ovr_ele_loc.(cur_ele_hms{iF}) ; [strcat(sbj(1:5),' ',cur_ele_loc(:,1)) cur_ele_loc(:,2:4)]];
            
            if has_aud && has_vis
                
                if fle; sig_chn(find(sig_chn(:,1)),2) = 0; sig_chn(find(sig_chn(:,3)),4) = 0; end
                rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
                ind = 1;
                for iA = 1:4
                    ovr_ele_act{1,iA} = [ovr_ele_act{1,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                    ind = ind+1;
                end
                
            elseif has_aud
                
                if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
                rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
                ind = 1;
                for iA = 3:4
                    ovr_ele_act{1,iA} = [ovr_ele_act{1,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                    ind = ind+1;
                end
                
            elseif has_vis
                
                if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
                rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
                ind = 1;
                for iA = 1:2
                    ovr_ele_act{1,iA} = [ovr_ele_act{1,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                    ind = ind+1;
                end
                
            end
            
            clear has_loc
            
        end
        
    else
        
        new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],['\t']);
        has_vis = any(strcmpi('visual_early',new_sig(1,:)));
        has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
        
        if has_aud && has_vis
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; sig_chn(find(sig_chn(:,3)),4) = 0; end
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            ind = 1;
            for iA = 1:4
                ovr_ele_act{2,iA} = [ovr_ele_act{2,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                ind = ind+1;
            end
            
        elseif has_aud
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            ind = 1;
            for iA = 3:4
                ovr_ele_act{2,iA} = [ovr_ele_act{2,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                ind = ind+1;
            end
            
        elseif has_vis
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            ind = 1;
            for iA = 1:2
                ovr_ele_act{2,iA} = [ovr_ele_act{2,iA} ; strcat(sbj(1:5),' ',rep_sel_ele{ind})];
                ind = ind+1;
            end
            
        end
        
    end
    
    clear has_loc
    
end

cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location_mni/' '/'  'lhs_combined_mni.csv'],ovr_ele_loc.lhs,[' ']);
cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location_mni/' '/'  'rhs_combined_mni.csv'],ovr_ele_loc.rhs,[' ']);

cfg = [];
cfg.elec_text = {['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location_mni/' '/'  'lhs_combined_mni.csv']  ... 
    ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location_mni/' '/'  'rhs_combined_mni.csv']};
cfg.pial_mat  = { '/space/mdeh1/5/halgdev/projects/nyuproj/loc/mni_template/template_surf/ch2_template_mni_lh_pial.mat'  ... 
    '/space/mdeh1/5/halgdev/projects/nyuproj/loc/mni_template/template_surf/ch2_template_mni_rh_pial.mat'};
cfg.hms       = {'lhs' 'rhs'};
cfg.sel_lbl   = {'visual_early' 'visual_late' 'auditory_early' 'auditory_late'};
cfg.hem       = {'lhs' 'rhs'};
cfg.lab_pos   = [1 2 6 7];
cfg.sel_ele   = ovr_ele_act;
cfg.col       = {rgb('light red') rgb('dark red') rgb('light blue') rgb('dark blue')};
cfg.sve_img   = 'jpg';
cfg.sve_loc   = [out_dir '/'];
cfg.sve_pre   = ['lfp_rep_loc'];
mmil_ieeg_sensor_location_plot_v1(cfg);

%% Make Side by Side Plots of Relevant Channels






