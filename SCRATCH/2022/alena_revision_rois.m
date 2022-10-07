clear; clc;

%% Inputs
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Collaborations/Alena/bil_net_epd/Revisions/ROIs/';

roi_sve = { 'Language' 'EF'  };
roi_nme = { {{ 'parsopercularis' 'parstriangularis' 'supramarginal' 'middletemporal' 'superiortemporal' } {'superiortemporal'}} ...
            {{ 'caudalmiddlefrontal' 'caudalanteriorcingulate' 'rostralanteriorcingulate' } {'caudalmiddlefrontal' 'caudalanteriorcingulate' 'rostralanteriorcingulate'}} };

%% Constants
cfg.sph = { 'lh' 'rh' };
cfg.sph_vew = { 'lat' 'ven' 'med' };

%% surf loading
% .pial files
srf_brn{1} = fs_read_surf('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/surf/lh.pial');
srf_brn{1}.coords = srf_brn{1}.vertices;
srf_brn{1}.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/surf/rh.pial');
srf_brn{2}.coords = srf_brn{2}.vertices;
srf_brn{2}.faces = srf_brn{2}.faces;

% Desikan
[lhs_roi_num, lhs_roi_lbl, lhs_col_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/lh.aparc.annot');
[rhs_roi_num, rhs_roi_lbl, rhs_col_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/rh.aparc.annot');

% Destrieaux
[lhs_roi_dst_num, lhs_roi_dst_lbl, lhs_col_dst_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/lh.aparc.a2009s.annot');
[rhs_roi_dst_num, rhs_roi_dst_lbl, rhs_col_dst_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/rh.aparc.a2009s.annot');

% HCP
[lhs_roi_hcp_num, lhs_roi_hcp_lbl, lhs_col_hcp_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/lh.HCPMMP1.annot');
[rhs_roi_hcp_num, rhs_roi_hcp_lbl, rhs_col_hcp_tab] = fs_read_annotation('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/fsaverage_xhemi/label/rh.HCPMMP1.annot');


%% Desikan to Destrieaux Overlap
for iR = 1:numel(roi_sve)
    
    [ ~, lhs_roi_des_ind ] = intersect(lhs_roi_lbl,roi_nme{iR}{1});
    lhs_roi_des_ind = find(ismember(lhs_roi_num,lhs_roi_des_ind));
    
    [ ~, rhs_roi_des_ind ] = intersect(rhs_roi_lbl,roi_nme{iR}{2});
    rhs_roi_des_ind = find(ismember(rhs_roi_num,rhs_roi_des_ind));
    
    % Plot Desikan
    cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
    cfg.prc_nme = '.aparc.annot';
    
    fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
    
    plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.0 0.4 0.2] [0.0 0.2 0.4 0.4]} {[0.4 0.6 0.4 0.4] [0.4 0.0 0.4 0.2] [0.4 0.2 0.4 0.4]}};
    
    col_hld{1} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{1}(lhs_roi_des_ind,:) = repmat(rgb('light green'),numel(lhs_roi_des_ind),1);
    col_hld{2} = repmat(rgb('light grey'), size(srf_brn{2}.vertices,1), 1 );
    col_hld{2}(rhs_roi_des_ind,:) = repmat(rgb('light green'),numel(rhs_roi_des_ind),1);
    
    
    for iH = 1:numel(cfg.sph)
        for iSP = 1:numel(cfg.sph_vew)
            
            pcfg = [];
            
            pcfg.fig_hdl = fig_hld(1);
            
            pcfg.surf_brain  = srf_brn{iH};
            pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
            
            pcfg.sph         = cfg.sph{iH};
            pcfg.sph_vew     = cfg.sph_vew{iSP};
            
            pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
            
            pcfg.col_map = col_hld{iH};
            
            pcfg.top_pct = 1;
            
            if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
            
            pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
            
            pcfg.sve_img = 1;
            
            ejk_surf_plot(pcfg);
            
        end
    end
    
    print( [ out_dir '/' 'Desikan' '_' roi_sve{iR} '.png'] ,'-dpng')
    close all  
    
    % Figure Out Destrieaux
    roi_dst_nme{iR}{1} = unique(lhs_roi_dst_num(lhs_roi_des_ind));
    roi_dst_nme{iR}{1} = lhs_roi_dst_lbl(roi_dst_nme{iR}{1});
    roi_dst_nme{iR}{2} = unique(rhs_roi_dst_num(rhs_roi_des_ind));
    roi_dst_nme{iR}{2} = rhs_roi_dst_lbl(roi_dst_nme{iR}{2});
    
    [ ~, lhs_roi_dst_ind ] = intersect(lhs_roi_dst_lbl,roi_dst_nme{iR}{1});
    lhs_roi_dst_ind = find(ismember(lhs_roi_dst_num,lhs_roi_dst_ind));
    
    [ ~, rhs_roi_dst_ind ] = intersect(rhs_roi_dst_lbl,roi_dst_nme{iR}{2});
    rhs_roi_dst_ind = find(ismember(rhs_roi_dst_num,rhs_roi_dst_ind));
    
    % Plot Destrieaux
    cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
    cfg.prc_nme = '.aparc.annot';
    
    fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
    
    plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.0 0.4 0.2] [0.0 0.2 0.4 0.4]} {[0.4 0.6 0.4 0.4] [0.4 0.0 0.4 0.2] [0.4 0.2 0.4 0.4]}};
    
    lhs_ovr_lap_ind = intersect(lhs_roi_dst_ind,lhs_roi_des_ind);
    lhs_old_mss_ind = setxor(lhs_roi_des_ind,lhs_ovr_lap_ind);
    lhs_new_add_ind = setxor(lhs_roi_dst_ind,lhs_ovr_lap_ind);
    
    col_hld{1} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{1}(lhs_ovr_lap_ind,:) = repmat(rgb('light green'),numel(lhs_ovr_lap_ind),1);
    col_hld{1}(lhs_old_mss_ind,:) = repmat(rgb('light red'),numel(lhs_old_mss_ind),1);
    col_hld{1}(lhs_new_add_ind,:) = repmat(rgb('light purple'),numel(lhs_new_add_ind),1);
    
    rhs_ovr_lap_ind = intersect(rhs_roi_dst_ind,rhs_roi_des_ind);
    rhs_old_mss_ind = setxor(rhs_roi_des_ind,rhs_ovr_lap_ind);
    rhs_new_add_ind = setxor(rhs_roi_dst_ind,rhs_ovr_lap_ind);
    
    col_hld{2} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{2}(rhs_ovr_lap_ind,:) = repmat(rgb('light green'),numel(rhs_ovr_lap_ind),1);
    col_hld{2}(rhs_old_mss_ind,:) = repmat(rgb('light red'),numel(rhs_old_mss_ind),1);
    col_hld{2}(rhs_new_add_ind,:) = repmat(rgb('light purple'),numel(rhs_new_add_ind),1);
    
    
    for iH = 1:numel(cfg.sph)
        for iSP = 1:numel(cfg.sph_vew)
            
            pcfg = [];
            
            pcfg.fig_hdl = fig_hld(1);
            
            pcfg.surf_brain  = srf_brn{iH};
            pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
            
            pcfg.sph         = cfg.sph{iH};
            pcfg.sph_vew     = cfg.sph_vew{iSP};
            
            pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
            
            pcfg.col_map = col_hld{iH};
            
            pcfg.top_pct = 1;
            
            if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
            
            pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
            
            pcfg.sve_img = 1;
            
            ejk_surf_plot(pcfg);
            
        end
    end
    
    print( [ out_dir '/' 'Destrieaux' '_' roi_sve{iR} '.png'] ,'-dpng')
    close all 
    
    % Put out report on ROIs
    out_csv = cell(numel(cat(1,roi_dst_nme{iR}{:})),5);
    out_csv_lbl = { 'roi' '% Overlap' '% Novel' '# Overlap' '# Novel' };
    roi_cnt = 1;
    for iH = 1:numel(cfg.sph)
        for iRO = 1:numel(roi_dst_nme{iR}{iH})
            out_csv{roi_cnt,1} = [ cfg.sph{iH} '_' roi_dst_nme{iR}{iH}{iRO}];
            
            if iH==1
                roi_ind_use = find(ismember( lhs_roi_dst_num, find(strcmpi(lhs_roi_dst_lbl,roi_dst_nme{iR}{iH}{iRO}))));
                out_csv{roi_cnt,4} = sum(ismember(lhs_ovr_lap_ind,roi_ind_use));
                out_csv{roi_cnt,2} = round((out_csv{roi_cnt,4} / numel(roi_ind_use) )*100);
                out_csv{roi_cnt,5} = sum(ismember(lhs_new_add_ind,roi_ind_use));
                out_csv{roi_cnt,3} = round((out_csv{roi_cnt,5} / numel(roi_ind_use) )*100);
            elseif iH==2
                roi_ind_use = find(ismember( rhs_roi_dst_num, find(strcmpi(rhs_roi_dst_lbl,roi_dst_nme{iR}{iH}{iRO}))));
                out_csv{roi_cnt,4} = sum(ismember(rhs_ovr_lap_ind,roi_ind_use));
                out_csv{roi_cnt,2} = round((out_csv{roi_cnt,4} / numel(roi_ind_use) )*100);
                out_csv{roi_cnt,5} = sum(ismember(rhs_new_add_ind,roi_ind_use));
                out_csv{roi_cnt,3} = round((out_csv{roi_cnt,5} / numel(roi_ind_use) )*100);
            end
            
            roi_cnt = roi_cnt+1;
        end
    end
    cell2csv([ out_dir '/' 'Destrieaux' '_' roi_sve{iR} '.csv'],[out_csv_lbl; out_csv])
end

%% Desikan to HCP Overlap
for iR = 1:numel(roi_sve)
    
    [ ~, lhs_roi_des_ind ] = intersect(lhs_roi_lbl,roi_nme{iR}{1});
    lhs_roi_des_ind = find(ismember(lhs_roi_num,lhs_roi_des_ind));
    
    [ ~, rhs_roi_des_ind ] = intersect(rhs_roi_lbl,roi_nme{iR}{2});
    rhs_roi_des_ind = find(ismember(rhs_roi_num,rhs_roi_des_ind));
    
    % Plot Desikan
    cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
    cfg.prc_nme = '.aparc.annot';
    
    fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
    
    plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.0 0.4 0.2] [0.0 0.2 0.4 0.4]} {[0.4 0.6 0.4 0.4] [0.4 0.0 0.4 0.2] [0.4 0.2 0.4 0.4]}};
    
    col_hld{1} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{1}(lhs_roi_des_ind,:) = repmat(rgb('light green'),numel(lhs_roi_des_ind),1);
    col_hld{2} = repmat(rgb('light grey'), size(srf_brn{2}.vertices,1), 1 );
    col_hld{2}(rhs_roi_des_ind,:) = repmat(rgb('light green'),numel(rhs_roi_des_ind),1);
    
    
    for iH = 1:numel(cfg.sph)
        for iSP = 1:numel(cfg.sph_vew)
            
            pcfg = [];
            
            pcfg.fig_hdl = fig_hld(1);
            
            pcfg.surf_brain  = srf_brn{iH};
            pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
            
            pcfg.sph         = cfg.sph{iH};
            pcfg.sph_vew     = cfg.sph_vew{iSP};
            
            pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
            
            pcfg.col_map = col_hld{iH};
            
            pcfg.top_pct = 1;
            
            if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
            
            pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
            
            pcfg.sve_img = 1;
            
            ejk_surf_plot(pcfg);
            
        end
    end
    
    print( [ out_dir '/' 'Desikan' '_' roi_sve{iR} '.png'] ,'-dpng')
    close all  
    
    % Figure Out HCP
    roi_hcp_nme{iR}{1} = unique(lhs_roi_hcp_num(lhs_roi_des_ind)); roi_hcp_nme{iR}{1}(roi_hcp_nme{iR}{1}==0) = [];
    roi_hcp_nme{iR}{1} = lhs_roi_hcp_lbl(roi_hcp_nme{iR}{1});
    roi_hcp_nme{iR}{2} = unique(rhs_roi_hcp_num(rhs_roi_des_ind)); roi_hcp_nme{iR}{2}(roi_hcp_nme{iR}{2}==0) = [];
    roi_hcp_nme{iR}{2} = rhs_roi_hcp_lbl(roi_hcp_nme{iR}{2});
    
    [ ~, lhs_roi_hcp_ind ] = intersect(lhs_roi_hcp_lbl,roi_hcp_nme{iR}{1});
    lhs_roi_hcp_ind = find(ismember(lhs_roi_hcp_num,lhs_roi_hcp_ind));
    
    [ ~, rhs_roi_hcp_ind ] = intersect(rhs_roi_hcp_lbl,roi_hcp_nme{iR}{2});
    rhs_roi_hcp_ind = find(ismember(rhs_roi_hcp_num,rhs_roi_hcp_ind));
    
    % Plot Destrieaux
    cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
    cfg.prc_nme = '.aparc.annot';
    
    fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
    
    plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.0 0.4 0.2] [0.0 0.2 0.4 0.4]} {[0.4 0.6 0.4 0.4] [0.4 0.0 0.4 0.2] [0.4 0.2 0.4 0.4]}};
    
    lhs_ovr_lap_ind = intersect(lhs_roi_hcp_ind,lhs_roi_des_ind);
    lhs_old_mss_ind = setxor(lhs_roi_des_ind,lhs_ovr_lap_ind);
    lhs_new_add_ind = setxor(lhs_roi_hcp_ind,lhs_ovr_lap_ind);
    
    col_hld{1} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{1}(lhs_ovr_lap_ind,:) = repmat(rgb('light green'),numel(lhs_ovr_lap_ind),1);
    col_hld{1}(lhs_old_mss_ind,:) = repmat(rgb('light red'),numel(lhs_old_mss_ind),1);
    col_hld{1}(lhs_new_add_ind,:) = repmat(rgb('light purple'),numel(lhs_new_add_ind),1);
    
    rhs_ovr_lap_ind = intersect(rhs_roi_hcp_ind,rhs_roi_des_ind);
    rhs_old_mss_ind = setxor(rhs_roi_des_ind,rhs_ovr_lap_ind);
    rhs_new_add_ind = setxor(rhs_roi_hcp_ind,rhs_ovr_lap_ind);
    
    col_hld{2} = repmat(rgb('light grey'), size(srf_brn{1}.vertices,1), 1 );
    col_hld{2}(rhs_ovr_lap_ind,:) = repmat(rgb('light green'),numel(rhs_ovr_lap_ind),1);
    col_hld{2}(rhs_old_mss_ind,:) = repmat(rgb('light red'),numel(rhs_old_mss_ind),1);
    col_hld{2}(rhs_new_add_ind,:) = repmat(rgb('light purple'),numel(rhs_new_add_ind),1);
    
    
    for iH = 1:numel(cfg.sph)
        for iSP = 1:numel(cfg.sph_vew)
            
            pcfg = [];
            
            pcfg.fig_hdl = fig_hld(1);
            
            pcfg.surf_brain  = srf_brn{iH};
            pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
            
            pcfg.sph         = cfg.sph{iH};
            pcfg.sph_vew     = cfg.sph_vew{iSP};
            
            pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
            
            pcfg.col_map = col_hld{iH};
            
            pcfg.top_pct = 1;
            
            if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
            
            pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
            
            pcfg.sve_img = 1;
            
            ejk_surf_plot(pcfg);
            
        end
    end
    
    print( [ out_dir '/' 'HCP' '_' roi_sve{iR} '.png'] ,'-dpng')
    close all 
    
    % Put out report on ROIs
    out_csv = cell(numel(cat(1,roi_hcp_nme{iR}{:})),5);
    out_csv_lbl = { 'roi' '% Overlap' '% Novel' '# Overlap' '# Novel' };
    roi_cnt = 1;
    for iH = 1:numel(cfg.sph)
        for iRO = 1:numel(roi_hcp_nme{iR}{iH})
            out_csv{roi_cnt,1} = [ cfg.sph{iH} '_' roi_hcp_nme{iR}{iH}{iRO}];
            
            if iH==1
                roi_ind_use = find(ismember( lhs_roi_hcp_num, find(strcmpi(lhs_roi_hcp_lbl,roi_hcp_nme{iR}{iH}{iRO}))));
                out_csv{roi_cnt,4} = sum(ismember(lhs_ovr_lap_ind,roi_ind_use));
                out_csv{roi_cnt,2} = round((out_csv{roi_cnt,4} / numel(roi_ind_use) )*100);
                out_csv{roi_cnt,5} = sum(ismember(lhs_new_add_ind,roi_ind_use));
                out_csv{roi_cnt,3} = round((out_csv{roi_cnt,5} / numel(roi_ind_use) )*100);
            elseif iH==2
                roi_ind_use = find(ismember( rhs_roi_hcp_num, find(strcmpi(rhs_roi_hcp_lbl,roi_hcp_nme{iR}{iH}{iRO}))));
                out_csv{roi_cnt,4} = sum(ismember(rhs_ovr_lap_ind,roi_ind_use));
                out_csv{roi_cnt,2} = round((out_csv{roi_cnt,4} / numel(roi_ind_use) )*100);
                out_csv{roi_cnt,5} = sum(ismember(rhs_new_add_ind,roi_ind_use));
                out_csv{roi_cnt,3} = round((out_csv{roi_cnt,5} / numel(roi_ind_use) )*100);
            end
            
            roi_cnt = roi_cnt+1;
        end
    end
    cell2csv([ out_dir '/' 'HCP' '_' roi_sve{iR} '.csv'],[out_csv_lbl; out_csv])
end















