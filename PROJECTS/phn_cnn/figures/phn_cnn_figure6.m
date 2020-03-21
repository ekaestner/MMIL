
clear; clc;

lbl_nme = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Data_all3T_temporal_toall_labels.csv');

%% PCA importances
pca_imp = cell2mat(mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_importances.csv'));

pca_com = cell2mat(mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_pca_components.csv'));

[ gdd_pcs , gdd_pcs_ind] = sort(pca_imp);
    gdd_pcs     = flipud(gdd_pcs);
    gdd_pcs_ind = flipud(gdd_pcs_ind);
    [ gdd_pcs_ind gdd_pcs ];
    str_pcs = gdd_pcs_ind(gdd_pcs>0.1);
    oth_pcs = gdd_pcs_ind(gdd_pcs>0.0);

% Top 4, Get top 27 regions
for iPC = 1:numel(str_pcs)
    
    pcs_hld = pca_com(str_pcs(iPC),:);
    
    [ imp_pcs , imp_pcs_ind ] = sort(pcs_hld);
        imp_pcs     = flipud(imp_pcs');
        imp_pcs_ind = flipud(imp_pcs_ind');
        
    pca_con_hld_tp4{iPC} = [ lbl_nme(imp_pcs_ind(1:27)) num2cell(imp_pcs(1:27)) ];
    
end

% Top 9, Get top 27 regions
for iPC = 1:numel(oth_pcs)
    
    pcs_hld = pca_com(oth_pcs(iPC),:);
    
    [ imp_pcs , imp_pcs_ind ] = sort(pcs_hld);
        imp_pcs     = flipud(imp_pcs');
        imp_pcs_ind = flipud(imp_pcs_ind');
        
    pca_con_hld_tp9{iPC} = [ lbl_nme(imp_pcs_ind(1:27)) num2cell(imp_pcs(1:27)) ];
    
end

%% Get Overall Laterality Distributions
% Top 4
lat_out_hld_tp4 = zeros(numel(pca_con_hld_tp4),3);

for iPC = 1:numel(pca_con_hld_tp4)

    lhs_cnt = 0;
    rhs_cnt = 0;
    bil_cnt = 0;

    for iRW = 1:size(pca_con_hld_tp4{iPC},1)
        nme_hld = regexp(pca_con_hld_tp4{iPC}{iRW,1},'==','split');
        if strcmpi(nme_hld{1}(1:2),'rh') && strcmpi(nme_hld{2}(1:2),'rh')
            rhs_cnt = rhs_cnt + 1;
        elseif strcmpi(nme_hld{1}(1:2),'lh') && strcmpi(nme_hld{2}(1:2),'lh')
            lhs_cnt = lhs_cnt + 1;
        else
            bil_cnt = bil_cnt + 1;
        end
    end
    
    lat_out_hld_tp4(iPC,:) = [ lhs_cnt bil_cnt rhs_cnt ];
    
end

sum( lat_out_hld_tp4 )

% Top 9
lat_out_hld_tp9 = zeros(numel(pca_con_hld_tp9),3);

for iPC = 1:numel(pca_con_hld_tp9)

    lhs_cnt = 0;
    rhs_cnt = 0;
    bil_cnt = 0;

    for iRW = 1:size(pca_con_hld_tp9{iPC},1)
        nme_hld = regexp(pca_con_hld_tp9{iPC}{iRW,1},'==','split');
        if strcmpi(nme_hld{1}(1:2),'rh') && strcmpi(nme_hld{2}(1:2),'rh')
            rhs_cnt = rhs_cnt + 1;
        elseif strcmpi(nme_hld{1}(1:2),'lh') && strcmpi(nme_hld{2}(1:2),'lh')
            lhs_cnt = lhs_cnt + 1;
        else
            bil_cnt = bil_cnt + 1;
        end
    end
    
    lat_out_hld_tp9(iPC,:) = [ lhs_cnt bil_cnt rhs_cnt ];
    
end

sum( lat_out_hld_tp9 )

[ sum(lat_out_hld_tp4(:,1))/sum(sum(lat_out_hld_tp4)) sum(lat_out_hld_tp9(:,1))/sum(sum(lat_out_hld_tp9)) ]
[ sum(lat_out_hld_tp4(:,2))/sum(sum(lat_out_hld_tp4)) sum(lat_out_hld_tp9(:,2))/sum(sum(lat_out_hld_tp9)) ]
[ sum(lat_out_hld_tp4(:,3))/sum(sum(lat_out_hld_tp4)) sum(lat_out_hld_tp9(:,3))/sum(sum(lat_out_hld_tp9)) ]

% Top4 BAR PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.ydt     = { sum(lat_out_hld_tp4(:,1)) sum(lat_out_hld_tp4(:,2)) sum(lat_out_hld_tp4(:,3)) };
fcfg.xdt     = { 1                       2                       3 };

fcfg.fce_col = { rgb('black')            rgb('black')            rgb('black') };

fcfg.xlb     = { 'Left'                  'Bilateral'             'Right'};
fcfg.ylb = 'Edge Count';

fcfg.ylm = [0 45];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure6/';
fcfg.out_nme = 'LateralityCount_Top4';

ejk_bar(fcfg)

% Top9 BAR PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.ydt     = { sum(lat_out_hld_tp9(:,1)) sum(lat_out_hld_tp9(:,2)) sum(lat_out_hld_tp9(:,3)) };
fcfg.xdt     = { 1                       2                       3 };

fcfg.fce_col = { rgb('black')            rgb('black')            rgb('black') };

fcfg.xlb     = { 'Left'                  'Bilateral'             'Right'};
fcfg.ylb = 'Edge Count';

fcfg.ylm = [0 120];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure6/';
fcfg.out_nme = 'LateralityCount_Top9';

ejk_bar(fcfg)

%% Get Overall Region Distributions
lat_tmp = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' };
        
ven_med_tmp = { 'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
                'parahippocampal' 'entorhinal' };
  
lat_frn     = { 'inferior-precentral' 'middle-precentral' ...
                'parsorbitalis'       'parstriangularis'  'parsopercularis' ...
                'caudalmiddlefrontal' 'middle-middlefrontal' 'rostral-middlefrontal'};

occ_ctx     = { 'lateraloccipital' 'cuneus' 'lingual' 'pericalcarine' 'calcarine'};  
            
lat_par_ctx = { 'inferior-postcentral' 'middle-postcentral' 'supramarginal' ...
                 'inferiorparietal' 'superiorparietal' };  
             
sup_ctx     = { 'caudal-superiorfrontal' 'middle-superiorfrontal' 'rostral-superiorfrontal' ...
                'superior-precentral'    'superior-postcentral' ...
                'paracentral' };

% Top 4, Get top 27 regions
tmp_out_hld_tp4 = zeros(numel(pca_con_hld_tp4),7);
            
for iPC = 1:numel(pca_con_hld_tp4)
    
    lat_tmp_cnt = 0;
    med_tmp_cnt = 0;
    lat_frn_cnt = 0;
    occ_ctx_cnt = 0;
    lat_par_cnt = 0;
    sup_ctx_cnt = 0;
    
    oth_cnt = 0;
    oth_hld = {};
    
    for iRW = 1:size(pca_con_hld_tp4{iPC},1)
        
        nme_hld = regexp(pca_con_hld_tp4{iPC}{iRW,1},'==','split');
        
        if  any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end)))
            lat_tmp_cnt = lat_tmp_cnt + sum( [ any(strcmpi(lat_tmp,nme_hld{1}(4:end))) any(strcmpi(lat_tmp,nme_hld{2}(4:end))) ] );
        end
        if  any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end)))
            med_tmp_cnt = med_tmp_cnt + sum( [ any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) any(strcmpi(ven_med_tmp,nme_hld{2}(4:end))) ] );
        end
        if  any(strcmpi(lat_frn,nme_hld{1}(4:end)))     || any(strcmpi(lat_frn,nme_hld{2}(4:end)))
            lat_frn_cnt = lat_frn_cnt + 1;
        end
        if  any(strcmpi(occ_ctx,nme_hld{1}(4:end)))     || any(strcmpi(occ_ctx,nme_hld{2}(4:end)))
            occ_ctx_cnt = occ_ctx_cnt + 1;
        end
        if  any(strcmpi(lat_par_ctx,nme_hld{1}(4:end))) || any(strcmpi(lat_par_ctx,nme_hld{2}(4:end)))
            lat_par_cnt = lat_par_cnt + 1;
        end
        if  any(strcmpi(sup_ctx,nme_hld{1}(4:end)))     || any(strcmpi(sup_ctx,nme_hld{2}(4:end)))
            sup_ctx_cnt = sup_ctx_cnt + 1;
        end
        
        if ~(any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end))))     && ...
                ~(any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end))))
            
            oth_cnt = oth_cnt + 1;
            oth_hld = { oth_hld{:} [nme_hld{1} '--' nme_hld{2}] };
            
        end
                        
    end
    
    tmp_out_hld_tp4(iPC,:) = [ lat_tmp_cnt med_tmp_cnt lat_frn_cnt occ_ctx_cnt lat_par_cnt sup_ctx_cnt oth_cnt ];
    
end

sum( tmp_out_hld_tp4 )

% Top 9, Get top 27 regions
tmp_out_hld_tp9 = zeros(numel(pca_con_hld_tp9),7);
            
for iPC = 1:numel(pca_con_hld_tp9)
    
    lat_tmp_cnt = 0;
    med_tmp_cnt = 0;
    lat_frn_cnt = 0;
    occ_ctx_cnt = 0;
    lat_par_cnt = 0;
    sup_ctx_cnt = 0;
    
    oth_cnt = 0;
    oth_hld = {};
    
    for iRW = 1:size(pca_con_hld_tp9{iPC},1)
        
        nme_hld = regexp(pca_con_hld_tp9{iPC}{iRW,1},'==','split');
        
        if  any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end)))
            lat_tmp_cnt = lat_tmp_cnt + sum( [ any(strcmpi(lat_tmp,nme_hld{1}(4:end))) any(strcmpi(lat_tmp,nme_hld{2}(4:end))) ] );
        end
        if  any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end)))
            med_tmp_cnt = med_tmp_cnt + sum( [ any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) any(strcmpi(ven_med_tmp,nme_hld{2}(4:end))) ] );
        end
        if  any(strcmpi(lat_frn,nme_hld{1}(4:end)))     || any(strcmpi(lat_frn,nme_hld{2}(4:end)))
            lat_frn_cnt = lat_frn_cnt + 1;
        end
        if  any(strcmpi(occ_ctx,nme_hld{1}(4:end)))     || any(strcmpi(occ_ctx,nme_hld{2}(4:end)))
            occ_ctx_cnt = occ_ctx_cnt + 1;
        end
        if  any(strcmpi(lat_par_ctx,nme_hld{1}(4:end))) || any(strcmpi(lat_par_ctx,nme_hld{2}(4:end)))
            lat_par_cnt = lat_par_cnt + 1;
        end
        if  any(strcmpi(sup_ctx,nme_hld{1}(4:end)))     || any(strcmpi(sup_ctx,nme_hld{2}(4:end)))
            sup_ctx_cnt = sup_ctx_cnt + 1;
        end
        
        if ~(any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end))))     && ...
                ~(any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end))))
            
            oth_cnt = oth_cnt + 1;
            oth_hld = { oth_hld{:} [nme_hld{1} '--' nme_hld{2}] };
            
        end
                        
    end
    
    tmp_out_hld_tp9(iPC,:) = [ lat_tmp_cnt med_tmp_cnt lat_frn_cnt occ_ctx_cnt lat_par_cnt sup_ctx_cnt oth_cnt ];
    
end

sum( tmp_out_hld_tp9 )

[ sum(tmp_out_hld_tp4(:,1))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,1))/sum(sum(tmp_out_hld_tp9)) ]
[ sum(tmp_out_hld_tp4(:,2))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,2))/sum(sum(tmp_out_hld_tp9)) ]
[ sum(tmp_out_hld_tp4(:,3))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,3))/sum(sum(tmp_out_hld_tp9)) ]
[ sum(tmp_out_hld_tp4(:,4))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,2))/sum(sum(tmp_out_hld_tp9)) ]
[ sum(tmp_out_hld_tp4(:,5))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,5))/sum(sum(tmp_out_hld_tp9)) ]
[ sum(tmp_out_hld_tp4(:,6))/sum(sum(tmp_out_hld_tp4)) sum(tmp_out_hld_tp9(:,6))/sum(sum(tmp_out_hld_tp9)) ]

% Colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dst_col = distinguishable_colors(10);

col{1} = rgb('grey');
col{2} = dst_col(1,:);
col{3} = dst_col(2,:);
col{4} = dst_col(3,:);
col{5} = dst_col(4,:);
col{6} = dst_col(5,:);
col{7} = dst_col(6,:);

col_map = [];
top_pct = 1;
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
end
top_pct = numel(col)-1;

% BAR PLOT top 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.ydt     = { sum(tmp_out_hld_tp4(:,1)) sum(tmp_out_hld_tp4(:,2)) sum(tmp_out_hld_tp4(:,3)) sum(tmp_out_hld_tp4(:,4)) sum(tmp_out_hld_tp4(:,5)) sum(tmp_out_hld_tp4(:,6)) };
fcfg.xdt     = { 1                       2                       3                       4                       5                       6 };

fcfg.fce_col = { dst_col(1,:)            dst_col(2,:)            dst_col(3,:)            dst_col(4,:)            dst_col(5,:)            dst_col(6,:)};

fcfg.xlb     = { 'Lateral Temporal'      'Ventral Temporal'      'Lateral Frontal'       'Occipital'             'Lateral Parietal'     'Superior Cortex' };
fcfg.ylb     = 'Edge Count';

fcfg.ylm     = [0 80];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure6/';
fcfg.out_nme = 'RegionCount_top4';

ejk_bar(fcfg)

% BAR PLOT top 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.ydt     = { sum(tmp_out_hld_tp9(:,1)) sum(tmp_out_hld_tp9(:,2)) sum(tmp_out_hld_tp9(:,3)) sum(tmp_out_hld_tp9(:,4)) sum(tmp_out_hld_tp9(:,5)) sum(tmp_out_hld_tp9(:,6)) };
fcfg.xdt     = { 1                       2                       3                       4                       5                       6 };

fcfg.fce_col = { dst_col(1,:)            dst_col(2,:)            dst_col(3,:)            dst_col(4,:)            dst_col(5,:)            dst_col(6,:)};

fcfg.xlb     = { 'Lateral Temporal'      'Ventral Temporal'      'Lateral Frontal'       'Occipital'             'Lateral Parietal'     'Superior Cortex' };
fcfg.ylb     = 'Edge Count';

fcfg.ylm     = [0 180];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure6/';
fcfg.out_nme = 'RegionCount_top9';

ejk_bar(fcfg)

% BRAIN PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
sph = { 'lh' 'rh' };
sph_vew = { 'lat' 'ven' 'med' };

% Make Data
[~,albl,~]=fs_read_annotation(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' sph{1} '.aparc.split.annot']);
pct_hld = [albl num2cell(zeros(size(albl,1),1))];

pct_hld(ismember(pct_hld(:,1),lat_tmp),2)     = {1};
pct_hld(ismember(pct_hld(:,1),ven_med_tmp),2) = {2};
pct_hld(ismember(pct_hld(:,1),lat_frn),2)     = {3};
pct_hld(ismember(pct_hld(:,1),occ_ctx),2)     = {4};
pct_hld(ismember(pct_hld(:,1),lat_par_ctx),2) = {5};
pct_hld(ismember(pct_hld(:,1),sup_ctx),2)     = {6};

% Load Surfaces
srf_brn{1} = fs_read_surf(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

% Make Plot
for iH = 1:2
    for iSP = 1:numel(sph_vew)
        
        fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
        
        pcfg = [];
        
        pcfg.surf_brain  = srf_brn{iH};
        pcfg.aparc       = ['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' sph{iH} '.aparc.split.annot'];
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iSP};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        if ~strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 1],'visible','off','Parent',fig_hld(1));
        elseif strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 0.6],'visible','off','Parent',fig_hld(1));
        end
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = col_map;
        
        pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
        pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
        
        pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
        print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure6/' '/' 'ROI_guide' '_' sph{iH} '_' sph_vew{iSP} '.png'],'-dpng','-r200')
        close all
        
    end
end
