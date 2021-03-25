% open iSASZ_plv_stt 

clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

%%
vis_cnt = [];
aud_cnt = [];
bim_cnt = [];

fst_nme = fieldnames(plv_hld);
for iN = 1:numel(fst_nme)
    
    scd_nme = fieldnames(plv_hld.(fst_nme{iN})); scd_nme(end) = [];
    for iS = 1:numel(scd_nme)
        
        vis_cnt = [vis_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){1}];
        aud_cnt = [aud_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){2}];
        bim_cnt = [bim_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){3}];
        
    end    
end

%% Make VISUAL Plot
%
vis_cnt = unique(vis_cnt);

vis_cnt_nme = cellfun(@(x) x(1:5),vis_cnt,'uni',0);

vis_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),vis_cnt,'uni',0),'uni',0);
for iE = 1:numel(vis_cnt_par); vis_cnt_par{iE}{1} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{1}]; vis_cnt_par{iE}{2} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{2}]; end

vis_cnt_ele = unique([vis_cnt_par{:}]);

%
% ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_vis_rep' '/' 'subjects' '/' 'total' '/' 'pap_vis_rep' '_plt']);
% ovr_txt = ovr_txt(2:end,:);
% ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_vis_act' '/' 'subjects' '/' 'total' '/' 'pap_vis_act' '_plt']);
ovr_txt = ovr_txt(2:end,:);
ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

% PLOT
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {vis_cnt_ele};
cfg.par_ele   = {vis_cnt_par};
cfg.all_ele   = ['NY503_SA_SZ_RPT06' vis_cnt_ele]; %ovr_txt;

cfg.sel_lbl = {'Visual'};
cfg.col     = {rgb('red')};
cfg.nsl_col = {rgb('dark red')};

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'visual_repetition'];
cfg.sve_pre   = ['middle_pic_visual'];
cfg.sep_str   = [','];

cfg.rad       = 2.35;

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Make AUDITORY Plot
%
aud_cnt = unique(aud_cnt);

aud_cnt_nme = cellfun(@(x) x(1:5),aud_cnt,'uni',0);

aud_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),aud_cnt,'uni',0),'uni',0);
for iE = 1:numel(aud_cnt_par); aud_cnt_par{iE}{1} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{1}]; aud_cnt_par{iE}{2} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{2}]; end

aud_cnt_ele = unique([aud_cnt_par{:}]);

%
% ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_aud_rep' '/' 'subjects' '/' 'total' '/' 'pap_aud_rep' '_plt']);
% ovr_txt = ovr_txt(2:end,:);
% ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_aud_act' '/' 'subjects' '/' 'total' '/' 'pap_aud_act' '_plt']);
ovr_txt = ovr_txt(2:end,:);
ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

% PLOT
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {aud_cnt_ele};
cfg.par_ele   = {aud_cnt_par};
cfg.all_ele   = [aud_cnt_ele 'NY503_SA_SZ_RPT06']; % ovr_txt;

cfg.sel_lbl = {'Auditory'};
cfg.col     = {rgb('blue')};
cfg.nsl_col = {rgb('dark blue')};

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'auditory_repetition'];
cfg.sve_pre   = ['middle_pic_auditory'];
cfg.sep_str   = [','];

cfg.rad       = 2.35;

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Make BIMODAL Plot
%
bim_cnt = unique(bim_cnt);

bim_cnt_nme = cellfun(@(x) x(1:5),bim_cnt,'uni',0);

bim_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),bim_cnt,'uni',0),'uni',0);
for iE = 1:numel(bim_cnt_par); bim_cnt_par{iE}{1} = [bim_cnt_nme{iE} '_SA_SZ_' bim_cnt_par{iE}{1}]; bim_cnt_par{iE}{2} = [bim_cnt_nme{iE} '_SA_SZ_' bim_cnt_par{iE}{2}]; end

bim_cnt_ele = unique([bim_cnt_par{:}]);

%
% ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_bim_rep' '/' 'subjects' '/' 'total' '/' 'pap_bim_rep' '_plt']);
% ovr_txt = ovr_txt(2:end,:);
% ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_bim_act' '/' 'subjects' '/' 'total' '/' 'pap_bim_act' '_plt']);
ovr_txt = ovr_txt(2:end,:);
ovr_txt = ovr_txt(logical(cell2mat(ovr_txt(:,3))),2);

% PLOT
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {bim_cnt_ele};
cfg.par_ele   = {bim_cnt_par};
cfg.all_ele   = [bim_cnt_ele 'NY503_SA_SZ_RPT06']; %ovr_txt;

cfg.sel_lbl = {'Bimodal'};
cfg.col     = {rgb('purple')};
cfg.nsl_col = {rgb('dark purple')};

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'bimodal_repetition'];
cfg.sve_pre   = ['middle_pic_bimodal'];
cfg.sep_str   = [','];

cfg.rad       = 2.35;

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Overlapping Bi-Modal
% VISUAL
vis_cnt = unique(vis_cnt);

vis_cnt_nme = cellfun(@(x) x(1:5),vis_cnt,'uni',0);

vis_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),vis_cnt,'uni',0),'uni',0);
for iE = 1:numel(vis_cnt_par); vis_cnt_par{iE}{1} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{1}]; vis_cnt_par{iE}{2} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{2}]; end

vis_cnt_ele = unique([vis_cnt_par{:}]);

vis_ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_vis_act' '/' 'subjects' '/' 'total' '/' 'pap_vis_act' '_plt']);
vis_ovr_txt = vis_ovr_txt(2:end,:);
vis_ovr_txt = vis_ovr_txt(logical(cell2mat(vis_ovr_txt(:,3))),2);

% AUDITORY
aud_cnt = unique(aud_cnt);

aud_cnt_nme = cellfun(@(x) x(1:5),aud_cnt,'uni',0);

aud_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),aud_cnt,'uni',0),'uni',0);
for iE = 1:numel(aud_cnt_par); aud_cnt_par{iE}{1} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{1}]; aud_cnt_par{iE}{2} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{2}]; end

aud_cnt_ele = unique([aud_cnt_par{:}]);

aud_ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_aud_act' '/' 'subjects' '/' 'total' '/' 'pap_aud_act' '_plt']);
aud_ovr_txt = aud_ovr_txt(2:end,:);
aud_ovr_txt = aud_ovr_txt(logical(cell2mat(aud_ovr_txt(:,3))),2);

% Put together Overlap
ovr_lap_cnt = intersect(vis_cnt_ele,aud_cnt_ele);

% Put together Vis
rmv_ind = [];
for iV = 1:numel(vis_cnt_par)
    if isempty(intersect(ovr_lap_cnt,vis_cnt_par{iV}))
        rmv_ind = [rmv_ind iV];
    end    
end
vis_cnt_par(rmv_ind) = [];

% Put together Aud
rmv_ind = [];
for iA = 1:numel(aud_cnt_par)
    if isempty(intersect(ovr_lap_cnt,aud_cnt_par{iA}))
        rmv_ind = [rmv_ind iA];
    end
end
aud_cnt_par(rmv_ind) = [];

ovr_txt = unique([vis_ovr_txt ; aud_ovr_txt]);

% Tidy up
vis_cnt_ele = unique(cat(2,vis_cnt_par{:}));
vis_cnt_ele = setxor(vis_cnt_ele,ovr_lap_cnt);

aud_cnt_ele = unique(cat(2,aud_cnt_par{:}));
aud_cnt_ele = setxor(aud_cnt_ele,ovr_lap_cnt);

% PLOT
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {vis_cnt_ele aud_cnt_ele ovr_lap_cnt};
cfg.sel_lbl = {'Visual'   'Auditory'  'Overlap'};
cfg.col     = {rgb('red') rgb('blue') rgb('bright purple')};
cfg.nsl_col = {rgb('dark purple')};

cfg.par_ele   = {vis_cnt_par aud_cnt_par};
cfg.par_ele_col = {rgb('reddish grey')+[0.15 0.075 0.075] rgb('bluish grey')+[0.1 0.1 0.2]};
cfg.par_ele_lne = [0.75                0.75];

cfg.all_ele   = [ovr_txt ; 'NY503_SA_SZ_RPT06']; %ovr_txt;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'overlap_repetition'];
cfg.sve_pre   = ['middle_pic_overlap'];
cfg.sep_str   = [','];

cfg.rad       = [1.5 1.5 2.5];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Alternate Bi-Modal
% %
% sbj_nme_hld = unique([vis_cnt_nme(:)' aud_cnt_nme(:)']);
% 
% for iSB = 1:numel(sbj_nme_hld)
% 
% vis_cnt = [];
% aud_cnt = [];
% bim_cnt = [];
% 
% fst_nme = fieldnames(plv_hld);
% for iN = 1:numel(fst_nme)
%     
%     scd_nme = fieldnames(plv_hld.(fst_nme{iN})); scd_nme(end) = [];
%     for iS = 1:numel(scd_nme)
%         
%         vis_cnt = [vis_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){1}];
%         aud_cnt = [aud_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){2}];
%         bim_cnt = [bim_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){3}];
%         
%     end    
% end
% 
% vis_cnt = unique(vis_cnt);
% aud_cnt = unique(aud_cnt);
% 
% vis_cnt_nme = cellfun(@(x) x(1:5),vis_cnt,'uni',0);
% aud_cnt_nme = cellfun(@(x) x(1:5),aud_cnt,'uni',0);
% 
% vis_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),vis_cnt,'uni',0),'uni',0);
% for iE = 1:numel(vis_cnt_par); vis_cnt_par{iE}{1} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{1}]; vis_cnt_par{iE}{2} = [vis_cnt_nme{iE} '_SA_SZ_' vis_cnt_par{iE}{2}]; end
% aud_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),aud_cnt,'uni',0),'uni',0);
% for iE = 1:numel(aud_cnt_par); aud_cnt_par{iE}{1} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{1}]; aud_cnt_par{iE}{2} = [aud_cnt_nme{iE} '_SA_SZ_' aud_cnt_par{iE}{2}]; end
% 
% vis_cnt_ele = unique([vis_cnt_par{:}]);
% aud_cnt_ele = unique([aud_cnt_par{:}]);
% 
% bim_ele = intersect(vis_cnt_ele,aud_cnt_ele);
% 
% bim_vis_cnt_par = [];
% for iV = 1:numel(vis_cnt_par)
%     if ~isempty(intersect(vis_cnt_par{iV},bim_ele))
%         bim_vis_cnt_par = [bim_vis_cnt_par vis_cnt_par(iV)];
%     end
% end
% bim_aud_cnt_par = [];
% for iA = 1:numel(aud_cnt_par)
%     if ~isempty(intersect(aud_cnt_par{iA},bim_ele))
%         bim_aud_cnt_par = [bim_aud_cnt_par aud_cnt_par(iA)];
%     end
% end
% 
% % Specific Subject
% sbj_nme_exc = sbj_nme_hld{iSB};
% 
% vis_cnt_ele = vis_cnt_ele(string_find(vis_cnt_ele,sbj_nme_exc));
% aud_cnt_ele = aud_cnt_ele(string_find(aud_cnt_ele,sbj_nme_exc));
% 
% bim_vis_cnt_par_sbj = [];
% for iV = 1:numel(bim_vis_cnt_par)
%     if ~isempty(string_find(bim_vis_cnt_par{iV},sbj_nme_exc))
%         bim_vis_cnt_par_sbj = [bim_vis_cnt_par_sbj bim_vis_cnt_par(iV)];
%     end
% end
% bim_aud_cnt_par_sbj = [];
% for iA = 1:numel(bim_aud_cnt_par)
%     if ~isempty(string_find(bim_aud_cnt_par{iA},sbj_nme_exc))
%         bim_aud_cnt_par_sbj = [bim_aud_cnt_par_sbj bim_aud_cnt_par(iA)];
%     end
% end
% 
% % PLOT
% cfg = [];
% 
% cfg.hms = {'lhs' 'rhs'};
% cfg.hem = {'lhs' 'rhs'};
% 
% cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
% cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};
% 
% cfg.sml_vew = 1;
% 
% cfg.sel_ele   = {vis_cnt_ele         aud_cnt_ele};
% cfg.par_ele   = {bim_vis_cnt_par_sbj bim_aud_cnt_par_sbj};
% cfg.all_ele   = {vis_cnt_ele{:} aud_cnt_ele{:} 'NY503_SA_SZ_RPT06'}; %ovr_txt;
% 
% cfg.sel_lbl = {'Visual' 'Auditory'};
% cfg.col     = {rgb('red') rgb('blue')};
% cfg.nsl_col = {rgb('black')};
% 
% cfg.sep_str         = ',';
% 
% cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
%                 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
%                 'Supramarginal' ...
%                 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
%                 'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
% 
% cfg.sve_img   = 'eps';
% cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'bimodal_liberal_repetition'];
% cfg.sve_pre   = [sbj_nme_exc '_' 'middle_pic_bimodal'];
% cfg.sep_str   = [','];
% 
% cfg.rad       = 2.35;
% 
% cfg.sve_img   = 'png';
% mmil_ieeg_sensor_location_plot_v4(cfg);
% 
% end

%% MAKE LINE PLOTS
% % Find Highlight Plots
% plv_hld.Precentral.tot
% plv_hld.Precentral.STG{1}
% plv_hld.Precentral.STG{2}
% 
% plv_hld.Precentral.Fusiform{1}
% plv_hld.Precentral.Fusiform{2}
% 
% % VISUAL
% 
% % AUDITORY
% 
% % BIMODAL






































