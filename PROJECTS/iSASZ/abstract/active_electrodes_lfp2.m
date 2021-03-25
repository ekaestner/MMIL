%% Create Early HGP divisions
clear; clc;

clr_fld    = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_in = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
sbj_nme    = mmil_readtext([clr_fld_in '/' 'subjects']);

indir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

abs_out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

for iS =  9:16 %[1:3 5:16]
    sbj     = sbj_nme{iS};
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    if has_vis && has_aud
        
        cfg = [];
        cfg.alt_stt = {'vis_ovr_all_stt' ...
                       'aud_ovr_all_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('reddish grey')) ...
            ft_stt_col(rgb('bluish grey'))};
        cfg.cmp_stt = [1 2];
        cfg.cmp_trl = 'ovr_all_evt';
        cfg.cmp     = {'101>999' ...
                       '102>999'};
        cfg.cmp_nme = {'visual_early' ...
                       'auditory_early'};
        cfg.tme_win = {[0.100 0.500] ...
            [0.100 0.500]};
        cfg.bse_tme = {[-0.300 0] ...
            [-0.300 0]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = strcat('ovr_all_',sem_dat.data_name);
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    elseif has_vis
        
        cfg = [];
        cfg.alt_trl = 'ovr_all_evt';
        cfg.alt_stt = {'vis_ovr_all_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('reddish grey'))};
        cfg.cmp_stt = [1];
        cfg.cmp_trl = 'ovr_all_evt';
        cfg.cmp     = {'101>999'};
        cfg.cmp_nme = {'visual_early'};
        cfg.tme_win = {[0.100 0.500]};
                cfg.bse_tme = {[-0.300 0]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = strcat('ovr_all_',sem_dat.data_name);
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    elseif has_aud
        
        cfg = [];
        cfg.alt_trl = 'ovr_all_evt';
        cfg.alt_stt = {'aud_ovr_all_stt'};
        cfg.alt_stt_col = {ft_stt_col(rgb('bluish grey'))};
        cfg.cmp_stt = [1];
        cfg.cmp_trl = 'ovr_all_evt';
        cfg.cmp     = {'102>999'};
        cfg.cmp_nme = {'auditory_early'};
        cfg.tme_win = {[0.100 0.500]};
        cfg.bse_tme = {[-0.300 0]};
        cfg.clr_fld  = clr_fld;
        cfg.sbj_nme  = sbj;
        cfg.typ      = strcat('ovr_all_',sem_dat.data_name);
        cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
        
    end
    
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.filename = [abs_out_dir '/' sbj '_hgp_data.mat'];
    ft_func([],cfg,sem_dat);
    
end

%% Plot subjects for checking
clear; clc;

clr_fld_new  = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
clr_fld_old = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/hgp';
indir   = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

sbj_nme    = mmil_readtext([clr_fld_old '/' 'subjects']);

for iS = [1:3 5:16]
    sbj     = sbj_nme{iS};
        
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})]);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    if has_vis && has_aud
        %
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
         
        cfg.cmp_ind = 1;
        cfg.cmb     = 1;
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})]; 
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
        cfg.prefix     = strcat('ovr_all_',sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
        %
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmp_ind = 2;
        cfg.cmb     = 1;
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})]; 
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
        cfg.prefix     = strcat('ovr_all_',sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_vis
    
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
         
        cfg.cmp_ind = 1;
        cfg.cmb     = 1;
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})]; 
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
        cfg.prefix     = strcat('ovr_all_',sem_dat.data_name{1},{'_'},pre_fix{1});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_aud
    
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             
        cfg.cmp_ind = 1;
        cfg.cmb     = 1;
        cfg.sig_fle = [clr_fld_new '/' 'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})]; 
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
        cfg.prefix     = strcat('ovr_all_',sem_dat.data_name{1},{'_'},pre_fix{1});
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

ele_hld = {'subject','v_act_ely','a_act_ely'};

for iS = [1:3 5:16]
    
    sbj = sbj_nme{iS};
    
    ele_hld{iS+1,1} = sbj;
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld_new '/'  'sig_chn' '/' sbj '/' strcat('ovr_all_',sem_dat.data_name{1})],['\t']);
    
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

        ele_hld{iS+1,2} = sum(cell2mat(new_sig(ind_icd,2)));
        ele_hld{iS+1,3} = sum(cell2mat(new_sig(ind_icd,3)));
                
    elseif has_vis
        
        ele_hld{iS+1,2} = sum(cell2mat(new_sig(ind_icd,2)));
               
    elseif has_aud
        
        ele_hld{iS+1,3} = sum(cell2mat(new_sig(ind_icd,3)));
        
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