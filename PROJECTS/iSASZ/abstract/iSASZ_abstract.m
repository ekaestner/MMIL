%% Move HGP to it's own folder
clear; clc;

clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
sbj_nme = mmil_readtext([clr_fld '/' 'subjects']);

indir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
abs_out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

for iS = 1:16
    sbj     = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_overall_data.mat'];
    sem_dat  = ft_func([],cfg);
       
    cfg = [];
    cfg.data_name = [1 3 4];
    cfg.rmfield   = 'yes';
    sem_dat       = ft_func([],cfg,sem_dat);
    
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.filename = [abs_out_dir '/' sbj '_hgp_data.mat'];
    ft_func([],cfg,sem_dat);

end

%% Create Early & Late HGP divisions
clear; clc;

clr_fld    = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_in = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
sbj_nme    = mmil_readtext([clr_fld_in '/' 'subjects']);

indir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

abs_out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

for iS = 1:16
    sbj     = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    if has_vis && has_aud
        
        cfg = [];
        cfg.alt_stt = {'vis_dif_rep_stt' ...
                       'aud_dif_rep_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('light red')) ft_stt_col(rgb('dark red')) ...
            ft_stt_col(rgb('light blue')) ft_stt_col(rgb('dark blue'))};
        cfg.cmp_stt = [1 1 2 2];
        cfg.cmp_trl = {'rep_all_eve'};
        cfg.cmp     = {'111!112' '111!112'...
                       '113!114' '113!114'};
        cfg.cmp_nme = {'visual_early' 'visual_late' ...
                       'auditory_early' 'auditory_late'};
        cfg.tme_win = {[0.100 0.500] [0.500 1.000] ...
            [0.100 0.500] [0.500 1.000]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = sem_dat.data_name;
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    elseif has_vis
        
        cfg = [];
        cfg.alt_stt     = {'vis_dif_rep_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('light red')) ft_stt_col(rgb('dark red'))};
        cfg.cmp_stt = [1 1];
        cfg.cmp_trl = {'rep_all_eve'};
        cfg.cmp     = {'111!112' '111!112'};
        cfg.cmp_nme = {'visual_early' 'visual_late'};
        cfg.tme_win = {[0.100 0.500] [0.500 1.000]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = sem_dat.data_name;
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    elseif has_aud
        
        cfg = [];
        cfg.alt_stt     = {'aud_dif_rep_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('light blue')) ft_stt_col(rgb('dark blue'))};
        cfg.cmp_stt = [1 1];
        cfg.cmp_trl = {'rep_all_eve'};
        cfg.cmp     = {'113!114'  '113!114' };
        cfg.cmp_nme = {'auditory_early' 'auditory_late'};
        cfg.tme_win = {[0.100 0.500] [0.500 1.000]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = sem_dat.data_name;
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    end
    
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.filename = [abs_out_dir '/' sbj '_hgp_data.mat'];
    ft_func([],cfg,sem_dat);
    
end

%% check against existing corrections
clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

indir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

abs_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

for iS = 1:16
    sbj     = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    old_sig = mmil_readtext([clr_fld_old '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}]);
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}]);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    if has_vis && has_aud
        
        new_sig(2:end,2) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,2),old_sig(2:end,3)));
        new_sig(2:end,3) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,3),old_sig(2:end,3)));
        
        new_sig(2:end,4) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,4),old_sig(2:end,6)));
        new_sig(2:end,5) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,5),old_sig(2:end,6)));
        
        cell2csv([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],new_sig);
        
    elseif has_vis
        
        new_sig(2:end,2) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,2),old_sig(2:end,3)));
        new_sig(2:end,3) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,3),old_sig(2:end,3)));      
        
        cell2csv([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],new_sig);
        
    elseif has_aud
        
        new_sig(2:end,2) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,2),old_sig(2:end,3)));
        new_sig(2:end,3) = num2cell(cellfun(@(x,y) double(all([x y])),new_sig(2:end,3),old_sig(2:end,3)));
        
        cell2csv([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],new_sig);
        
    end
    
end

%% Plot subjects for checking
clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/hgp';
indir   = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

for iS = 1:16
    sbj     = sbj_nme{iS};
        
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}]);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    if has_vis && has_aud
        %
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmb = 1;
        cfg.cmp_ind = 1:2;
        cfg.cmb_typ = {[1 2]};
        cfg.cmb_lbl = {'visual_early'};
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}]; 
        [grp_chn{1},pre_fix{1},grp_stt{1},grp_stt_col{1}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{1}));
        
        cfg.chn_grp         = grp_chn{1};
        cfg.alt_eve         = {'rep_all_eve'};
        cfg.eve             = [112 111];
        cfg.lnstyle.col_ord = {rgb('green') rgb('red')};
        cfg.cnd_nme         = {'V Nov' 'V Rep'};
        cfg.alt_lbl = 'vis_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{1};
        cfg.stt_col = grp_stt_col{1};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [out_dir '/' sbj '/' 'Visual' '/'];
        cfg.prefix     = strcat(sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
        %
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmb = 1; 
        cfg.cmb_typ = {[3 4]};
        cfg.cmp_ind = 3:4;
        cfg.cmb_lbl = {'auditory_early'};
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}]; 
        [grp_chn{1},pre_fix{1},grp_stt{1},grp_stt_col{1}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{1}));
        
        cfg.chn_grp         = grp_chn{1};
        cfg.alt_eve         = {'rep_all_eve'};
        cfg.eve             = [113 114];
        cfg.lnstyle.col_ord = {rgb('yellow') rgb('blue')};
        cfg.cnd_nme         = {'A Nov' 'A Rep'};
        cfg.alt_lbl = 'aud_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{1};
        cfg.stt_col = grp_stt_col{1};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [out_dir '/' sbj '/' 'Auditory' '/'];
        cfg.prefix     = strcat(sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_vis
    
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmb = 1;
        cfg.cmb_typ = {[1 2]};
        cfg.cmb_lbl = {'visual_early'};
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}]; 
        [grp_chn{1},pre_fix{1},grp_stt{1},grp_stt_col{1}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{1}));
        
        cfg.chn_grp         = grp_chn{1};
        cfg.alt_eve         = {'rep_all_eve'};
        cfg.eve             = [112 111];
        cfg.lnstyle.col_ord = {rgb('green') rgb('red')};
        cfg.cnd_nme         = {'V Nov' 'V Rep'};
        cfg.alt_lbl = 'vis_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{1};
        cfg.stt_col = grp_stt_col{1};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [out_dir '/' sbj '/' 'Visual' '/'];
        cfg.prefix     = strcat(sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_aud
    
                cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmb = 1; 
        cfg.cmb_typ = {[1 2]};
        cfg.cmb_lbl = {'auditory_early'};
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' sem_dat.data_name{1}]; 
        [grp_chn{1},pre_fix{1},grp_stt{1},grp_stt_col{1}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{1}));
        
        cfg.chn_grp         = grp_chn{1};
        cfg.alt_eve         = {'rep_all_eve'};
        cfg.eve             = [113 114];
        cfg.lnstyle.col_ord = {rgb('yellow') rgb('blue')};
        cfg.cnd_nme         = {'A Nov' 'A Rep'};
        cfg.alt_lbl = 'aud_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{1};
        cfg.stt_col = grp_stt_col{1};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [out_dir '/' sbj '/' 'Auditory' '/'];
        cfg.prefix     = strcat(sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    end
    
end

%% Caculate total numbers & effects of electrodes
clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/hgp';
indir   = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

ele_hld = {'subject','v_tot','v_rep_ely','v_rep_lte','a_tot','a_rep_ely','a_rep_lte','ely_ovr_lap','lte_ovr_lap'};

for iS = [1:3 5:16]
    
    sbj = sbj_nme{iS};
    
    ele_hld{iS+1,1} = sbj;
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],['\t']);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    num_rmv = sum(strcmpi(new_sig(~cellfun(@isempty,new_sig(:,end)),end),'remove'));
    
    ind_icd = 2:size(new_sig,1);  
    
    dpt_ind = find(~cellfun(@isempty,strfind(new_sig(:,1),'D')) | ~cellfun(@isempty,strfind(new_sig(:,1),'d'))); 
    if any(~cellfun(@isempty,strfind(new_sig(dpt_ind,1),'>')))
        str_ind = strfind(new_sig(dpt_ind,1),'>');
    else
        str_ind = strfind(new_sig(dpt_ind,1),'-');
    end
    if ~isempty(dpt_ind); dpt_ind = dpt_ind(cellfun(@(x,y) strcmpi(x(y+1),'d'),new_sig(dpt_ind,1),str_ind)); end
    ind_icd(dpt_ind-1) = [];
    
    if has_vis && has_aud

        ele_hld{iS+1,2} = numel(ind_icd)-num_rmv;
        ele_hld{iS+1,3} = sum(cell2mat(new_sig(ind_icd,2)));
        ele_hld{iS+1,4} = sum(cell2mat(new_sig(ind_icd,3)));
        
        ele_hld{iS+1,5} = numel(ind_icd)-num_rmv;
        ele_hld{iS+1,6} = sum(cell2mat(new_sig(ind_icd,4)));
        ele_hld{iS+1,7} = sum(cell2mat(new_sig(ind_icd,5)));
        
        new_sig(find(cell2mat(new_sig(ind_icd,2)))+1,3) = num2cell(0); new_sig(find(cell2mat(new_sig(2:end,4)))+1,5) = num2cell(0);
        ele_hld{iS+1,8} = sum(numel(intersect(find(cell2mat(new_sig(ind_icd,2))),find(cell2mat(new_sig(ind_icd,4)))))); 
        ele_hld{iS+1,9} = sum(numel(intersect(find(cell2mat(new_sig(ind_icd,3))),find(cell2mat(new_sig(ind_icd,5)))))); 
        
    elseif has_vis
        
        ele_hld{iS+1,2} = numel(ind_icd)-num_rmv;
        ele_hld{iS+1,3} = sum(cell2mat(new_sig(ind_icd,2)));
        ele_hld{iS+1,4} = sum(cell2mat(new_sig(ind_icd,3)));
        
        ele_hld{iS+1,5} = 0;
        ele_hld{iS+1,6} = 0;
        ele_hld{iS+1,7} = 0;
        
        ele_hld{iS+1,8} = 0; 
        ele_hld{iS+1,9} = 0; 
        
    elseif has_aud
        
        ele_hld{iS+1,2} = 0;
        ele_hld{iS+1,3} = 0;
        ele_hld{iS+1,4} = 0;
        
        ele_hld{iS+1,5} = numel(ind_icd)-num_rmv;
        ele_hld{iS+1,6} = sum(cell2mat(new_sig(ind_icd,2)));
        ele_hld{iS+1,7} = sum(cell2mat(new_sig(ind_icd,3)));       
        
        ele_hld{iS+1,8} = 0; 
        ele_hld{iS+1,9} = 0;
    end
    
end

v_tot     = sum(cell2mat(ele_hld(2:end,2)));
v_rep_ely = sum(cell2mat(ele_hld(2:end,3)));
v_rep_lte = sum(cell2mat(ele_hld(2:end,4)));

a_tot     = sum(cell2mat(ele_hld(2:end,5)));
a_rep_ely = sum(cell2mat(ele_hld(2:end,6)));
a_rep_lte = sum(cell2mat(ele_hld(2:end,7)));

ovr_lap_ely = sum(cell2mat(ele_hld(2:end,8)));
ovr_lap_lte = sum(cell2mat(ele_hld(2:end,9)));

%% Plot location plots
clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/hgp';
indir   = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

indir_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

ele_hld = {'subject','v_tot','v_rep_ely','v_rep_lte','a_tot','a_rep_ely','a_rep_lte'};

for iS = 13:16
    
    sbj = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
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
        sig_chn = mmil_readtext(cfg.sig_fle,[',']);
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
        
        new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}]);
        has_vis = any(strcmpi('visual_early',new_sig(1,:)));
        has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
        
        if has_aud && has_vis
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; sig_chn(find(sig_chn(:,3)),4) = 0; end
            rep_sel_lbl = sem_dat.(sem_dat.data_name{1}).cfg.alt_sig_chn.lbl;
            rep_lab_pos = [1 2 6 7];
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            rep_col     = {rgb('light red') rgb('dark red') rgb('light blue') rgb('dark blue')};
            
        elseif has_aud
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
            rep_sel_lbl = sem_dat.(sem_dat.data_name{1}).cfg.alt_sig_chn.lbl;
            rep_lab_pos = [1 2];
            rep_sel_ele = cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            rep_col     = {rgb('light blue') rgb('dark blue')};
            
        elseif has_vis
            
            if fle; sig_chn(find(sig_chn(:,1)),2) = 0; end
            rep_sel_lbl = sem_dat.(sem_dat.data_name{1}).cfg.alt_sig_chn.lbl;
            rep_lab_pos = [1 2];
            cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label(x),cellfun(@find,num2cell(sig_chn,1),'Uni',0),'Uni',0);
            rep_col     = {rgb('light red') rgb('dark red')};
            
        end
        
        cfg = [];
        cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/electrode_location' '/' sbj]);
        cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/pial' '/' sbj]);
        cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi' '/' sbj]);
        cfg.sel_lbl   = rep_sel_lbl;
        cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/hemi/' sbj]);
        cfg.lab_pos   = rep_lab_pos;
        cfg.sel_ele   = rep_sel_ele;
        cfg.col       = rep_col;
        cfg.sve_img   = 'jpg';
        cfg.sve_loc   = [out_dir '/' sbj];
        cfg.sve_pre   = [sbj 'hgp_rep'];
        mmil_ieeg_sensor_location_plot_v1(cfg);       
        
    end
    
    clear has_loc
    
end
















