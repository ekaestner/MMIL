function FW_figure2

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.tsk     = 'FW';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = {'pap_wrd_600' 'pap_con_600'};
fcfg.eff_clm = { [ 1 3 ] [1]};
fcfg.eff_col = { rgb('purple') rgb('bright red') rgb('reddish grey') };
fcfg.eff_lbl = { 'Letters' 'Words' 'FalseFont' };
fcfg.nsl_col = {rgb('white')};

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = 'pap_rsp_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 1; % 1

%% Middle Portion
% Selective Electrodes
eff_txt_one = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{1} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{1} '_plt']);
eff_txt_two = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{2} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{2} '_plt']);

fcfg.sel_ele{1} = eff_txt_one(find(cell2mat(eff_txt_one(2:end,fcfg.eff_clm{1}(1)+2)))+1,2);
fcfg.sel_ele{2} = eff_txt_one(find(cell2mat(eff_txt_one(2:end,fcfg.eff_clm{1}(2)+2)))+1,2);
fcfg.sel_ele{3} = eff_txt_two(find(cell2mat(eff_txt_two(2:end,fcfg.eff_clm{2}(1)+2)))+1,2);

% Overall Electrodes
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);

fcfg.all_ele = intersect(ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2),unique([fcfg.sel_ele{1} ; fcfg.sel_ele{2} ; fcfg.sel_ele{3}]));

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_ele = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.inc_prc = { [fcfg.clr_fld '/' 'electrode_location_files' '/' 'fsaverage' '/'  'label' '/' 'lh.aparc.split.annot'] ...
                [fcfg.clr_fld '/' 'electrode_location_files' '/' 'fsaverage' '/'  'label' '/' 'rh.aparc.split.annot'] };
            
mmil_ieeg_sensor_location_plot_v4(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_ele;
mmil_loc_dot(pcfg)

%% Middle Portion ALTERNATE
% Letter
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_wrd_600';
cfg.eff_clm = 1;
cfg.top_pct = 0.25;
cfg.eff_col = {{'dark purple' 'purple' 'neon purple'}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

% Word
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_wrd_600';

cfg.eff_clm = 3;
cfg.top_pct = 0.25;
cfg.eff_col = {{'dark red' 'red' 'neon red'}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

% False-Font
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_con_600';
cfg.eff_clm = 1;
cfg.top_pct = 0.25;
cfg.eff_col = {{[rgb('reddish grey') - [0.16 0.06 0.06]] rgb('reddish grey') [rgb('reddish grey') + [0.24 0.16 0.16]]}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

% Colorbar
pcfg = [];
pcfg.col_map = {'dark purple' 'purple' 'neon purple'};
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure2';
pcfg.sve_pre = 'letter';
mmil_color_bar(pcfg)

pcfg = [];
pcfg.col_map = {'dark red' 'red' 'neon red'};
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure2';
pcfg.sve_pre = 'word';
mmil_color_bar(pcfg)

pcfg = [];
pcfg.col_map = {[rgb('reddish grey') - [0.16 0.06 0.06]] rgb('reddish grey') [rgb('reddish grey') + [0.24 0.16 0.16]]};
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure2';
pcfg.sve_pre = 'falsefont';
mmil_color_bar(pcfg)

%% Cirlces
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'left'  '/' ])
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'right' '/' ])

tot_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/'  'table1.csv']);
tot_col = [4 11];

crc_lft_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/'  'table2.csv']);
crc_rgh_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/'  'table3.csv']);
col = [3 4 5 11];

for iT = 1:size(crc_lft_tbl,1)
    if isnumeric(crc_lft_tbl{iT,col(1)}) && ~isempty(crc_lft_tbl{iT,col(1)})
        lft_tot_scl(iT,1) = crc_lft_tbl{iT,col(1)}/tot_tbl{iT,tot_col(1)};
    end
end
for iT = 1:size(crc_rgh_tbl,1)
    if isnumeric(crc_rgh_tbl{iT,col(1)}) && ~isempty(crc_rgh_tbl{iT,col(1)})
        rgh_tot_scl(iT,1) = crc_rgh_tbl{iT,col(1)}/tot_tbl{iT,tot_col(2)};
    end
end

eff_col = { rgb('purple') rgb('bright red') rgb('reddish grey') };

ttt = [lft_tot_scl ; rgh_tot_scl];
max_num = max(ttt);

for iR = 3:size(crc_lft_tbl,1)
    
    if ~strcmpi(crc_lft_tbl{iR,col(2)},'-') && ~strcmpi(crc_lft_tbl{iR,col(3)},'-') && ~strcmpi(crc_lft_tbl{iR,col(4)},'-') ...
            && sum([crc_lft_tbl{iR,col(2)} crc_lft_tbl{iR,col(3)} crc_lft_tbl{iR,col(4)}])>0
        
        pfg = figure();
        pch = pie([crc_lft_tbl{iR,col(2)} crc_lft_tbl{iR,col(3)} crc_lft_tbl{iR,col(4)}],{'' '' ''});
        
        txt_axe = axes('Parent',pfg);
        text(0.5,1.05,['Left ' crc_lft_tbl{iR,1}],'FontSize',24,'HorizontalAlignment','center','Parent',txt_axe)
        text(0.5,0.995,['Language-Sensitive Electrodes : ' num2str(crc_lft_tbl{iR,col(1)}) ' / ' num2str(tot_tbl{iR,tot_col(1)})],'FontSize',16,'HorizontalAlignment','center','Parent',txt_axe)
        axis off
        
        pch_axe = get(pch(1),'Parent');
        pch_axe_pos = get(pch_axe,'Position');
        new_rad = sqrt(lft_tot_scl(iR,1)/max_num);
        new_rad = pch_axe_pos .* [1 1 new_rad new_rad ];
        new_rad(1) = pch_axe_pos(1)+((pch_axe_pos(3)-new_rad(3))/2);
        new_rad(2) = pch_axe_pos(2)+((pch_axe_pos(4)-new_rad(4))/2);
        set(pch_axe,'Position',new_rad)
        
        set(pch,'EdgeColor','none')
        hp = flipud(findobj(pfg,'Type','patch'));
        cnt = 1;
        for iC = 1:3
            if crc_lft_tbl{iR,col(iC+1)}>0
                set(hp(cnt),'FaceColor',eff_col{iC})
                cnt = cnt + 1;
            end
        end
        set(gcf, 'color','white','InvertHardCopy', 'off');
        axis off;
        
        cfg = [];
        cfg.jpg = 1;
        cfg.prn_typ = 'eps';
        cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'left'  '/' ['left_' crc_rgh_tbl{iR,1}] ];
        mmil_print_plot(cfg)
        close all
        
    end
    
    if ~strcmpi(crc_rgh_tbl{iR,col(2)},'-') && ~strcmpi(crc_rgh_tbl{iR,col(3)},'-') && ~strcmpi(crc_rgh_tbl{iR,col(4)},'-') ...
            && sum([crc_rgh_tbl{iR,col(2)} crc_rgh_tbl{iR,col(3)} crc_rgh_tbl{iR,col(4)}])>0
        
        pfg = figure();
        pch = pie([crc_rgh_tbl{iR,col(2)} crc_rgh_tbl{iR,col(3)} crc_rgh_tbl{iR,col(4)}],{'' '' ''});
        
        txt_axe = axes('Parent',pfg);
        text(0.5,1.05,['Right ' crc_rgh_tbl{iR,1}],'FontSize',24,'HorizontalAlignment','center','Parent',txt_axe)
        text(0.5,0.995,['Language-Sensitive Electrodes : ' num2str(crc_rgh_tbl{iR,col(1)}) ' / ' num2str(tot_tbl{iR,tot_col(2)})],'FontSize',16,'HorizontalAlignment','center','Parent',txt_axe)
        axis off
        
        pch_axe = get(pch(1),'Parent');
        pch_axe_pos = get(pch_axe,'Position');
        new_rad = sqrt(rgh_tot_scl(iR,1)/max_num);
        new_rad = pch_axe_pos .* [1 1 new_rad new_rad ];
        new_rad(1) = pch_axe_pos(1)+((pch_axe_pos(3)-new_rad(3))/2);
        new_rad(2) = pch_axe_pos(2)+((pch_axe_pos(4)-new_rad(4))/2);
        set(pch_axe,'Position',new_rad)
        
        set(pch,'EdgeColor','none')
        hp = flipud(findobj(pfg,'Type','patch'));
        cnt = 1;
        for iC = 1:3
            if crc_rgh_tbl{iR,col(iC+1)}>0
                set(hp(cnt),'FaceColor',eff_col{iC})
                cnt = cnt + 1;
            end
        end
        set(gcf, 'color','white','InvertHardCopy', 'off');
        axis off;
        
        cfg = [];
        cfg.jpg = 1;
        cfg.prn_typ = 'eps';
        cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'right'  '/' ['right_' crc_rgh_tbl{iR,1}]];
        mmil_print_plot(cfg)
        close all
        
    end
    
end

%% Cirlces ALTERNATE
%% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg_plt_nme = {'Lateral Occipital' ...
               'Fusiform' ...
               'Inferior Temporal Gyrus' ...
               'Middle Temporal Gyrus' ...
               'Superior Temporal Gyrus' ...
               'Lateral Parietal' ...
               'Supramarginal' ...
               'Precentral' ...
               'Postcentral' ...
               'Pars Opercularis' ...
               'Pars Triangularis' ...
               'Middle Frontal Gyrus' };           
           
reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
% reg.cmb_reg.orb = { 'Pars Orbitalis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
reg_nme = fieldnames(reg);
cmb_nme = fieldnames(reg.(reg_nme{1}));

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600' 'pap_rsp_600'};
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]         [2 3]};
fcfg.plt     = 0;
hms_sig = mmil_fisher_exact_region_test(fcfg);

% %%%%%%%%%%%%%%%%%%%%%%%
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'left_combined'  '/' ])
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'right_combined' '/' ])

for iT = 1:numel(cmb_nme)-1; lft_tot_scl(iT,1) = hms_sig.cmb_reg.Specific{iT,4}(1,1) / sum( hms_sig.cmb_reg.Specific{iT,4}(1,1:2)); end
for iT = 1:numel(cmb_nme)-1; rgh_tot_scl(iT,1) = hms_sig.cmb_reg.Specific{iT,4}(2,1) / sum( hms_sig.cmb_reg.Specific{iT,4}(2,1:2)); end

eff_col = { rgb('purple') rgb('bright red') rgb('reddish grey') };

ttt = [lft_tot_scl ; rgh_tot_scl];
max_num = max(ttt);

for iR = 1:numel(cmb_nme)-1
    
    % Left Hemisphere
    pfg = figure();
    pch = pie([hms_sig.cmb_reg.VisualOrthographic{iR,4}(1,1) hms_sig.cmb_reg.VisualTotal{iR,4}(1,1) hms_sig.cmb_reg.VisualSensory{iR,4}(1,1)],{'' '' ''});
    
    txt_axe = axes('Parent',pfg);
    text(0.5,1.05,['Left ' reg_plt_nme{iR}],'FontSize',24,'HorizontalAlignment','center','Parent',txt_axe)
    text(0.5,0.995,['Language-Sensitive Electrodes : ' num2str(hms_sig.cmb_reg.Specific{iR,4}(1,1)) ' / ' num2str(sum( hms_sig.cmb_reg.Specific{iR,4}(1,1:2)))],'FontSize',16,'HorizontalAlignment','center','Parent',txt_axe)
    axis off
    
    pch_axe = get(pch(1),'Parent');
    pch_axe_pos = get(pch_axe,'Position');
    new_rad = sqrt(lft_tot_scl(iR,1)/max_num);
    new_rad = pch_axe_pos .* [1 1 new_rad new_rad ];
    new_rad(1) = pch_axe_pos(1)+((pch_axe_pos(3)-new_rad(3))/2);
    new_rad(2) = pch_axe_pos(2)+((pch_axe_pos(4)-new_rad(4))/2);
    set(pch_axe,'Position',new_rad)
    
    set(pch,'EdgeColor','none')
    hp = flipud(findobj(pfg,'Type','patch'));
    cnt = 1;
    num_hld = [hms_sig.cmb_reg.VisualOrthographic{iR,4}(1,1) hms_sig.cmb_reg.VisualTotal{iR,4}(1,1) hms_sig.cmb_reg.VisualSensory{iR,4}(1,1)];
    for iC = 1:3
        if num_hld(iC)>0
            set(hp(cnt),'FaceColor',eff_col{iC})
            cnt = cnt + 1;
        end
    end
    set(gcf, 'color','white','InvertHardCopy', 'off');
    axis off;
    
    cfg = [];
    cfg.jpg = 1;
    cfg.prn_typ = 'eps';
    cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'left_combined'  '/' ['left_' reg_plt_nme{iR}] ];
    mmil_print_plot(cfg)
    close all    
    
    % Right Hemisphere
    pfg = figure();
    pch = pie([hms_sig.cmb_reg.VisualOrthographic{iR,4}(2,1) hms_sig.cmb_reg.VisualTotal{iR,4}(2,1) hms_sig.cmb_reg.VisualSensory{iR,4}(2,1)],{'' '' ''});
    
    txt_axe = axes('Parent',pfg);
    text(0.5,1.05,['Right ' reg_plt_nme{iR}],'FontSize',24,'HorizontalAlignment','center','Parent',txt_axe)
    text(0.5,0.995,['Language-Sensitive Electrodes : ' num2str(hms_sig.cmb_reg.Specific{iR,4}(2,1)) ' / ' num2str(sum( hms_sig.cmb_reg.Specific{iR,4}(2,1:2)))],'FontSize',16,'HorizontalAlignment','center','Parent',txt_axe)
    axis off
    
    pch_axe = get(pch(1),'Parent');
    pch_axe_pos = get(pch_axe,'Position');
    new_rad = sqrt(rgh_tot_scl(iR,1)/max_num);
    new_rad = pch_axe_pos .* [1 1 new_rad new_rad ];
    new_rad(1) = pch_axe_pos(1)+((pch_axe_pos(3)-new_rad(3))/2);
    new_rad(2) = pch_axe_pos(2)+((pch_axe_pos(4)-new_rad(4))/2);
    set(pch_axe,'Position',new_rad)
    
    set(pch,'EdgeColor','none')
    hp = flipud(findobj(pfg,'Type','patch'));
    cnt = 1;
    num_hld = [hms_sig.cmb_reg.VisualOrthographic{iR,4}(2,1) hms_sig.cmb_reg.VisualTotal{iR,4}(2,1) hms_sig.cmb_reg.VisualSensory{iR,4}(2,1)];
    for iC = 1:3
        if num_hld(iC)>0
            set(hp(cnt),'FaceColor',eff_col{iC})
            cnt = cnt + 1;
        end
    end
    set(gcf, 'color','white','InvertHardCopy', 'off');
    axis off;
    
    cfg = [];
    cfg.jpg = 1;
    cfg.prn_typ = 'eps';
    cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'right_combined'  '/' ['right_' reg_plt_nme{iR}] ];
    mmil_print_plot(cfg)
    close all  
    
end

end















