
clear; clc;

lbl_nme = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Data_all3T_temporal_toall_labels.csv');

%% PCA importances
pca_imp = cell2mat(mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_importances.csv'));

% pca_com = abs(cell2mat(mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_pca_components.csv')));
pca_com = cell2mat(mmil_readtext('/space/syn09/1/data/MMILDB/ekaestne/CONNECTOME/output/001_CONNECTOME_Outcomes_all3T_temporal_toall_pca_components.csv'));

[ gdd_pcs , gdd_pcs_ind] = sort(pca_imp);
    gdd_pcs     = flipud(gdd_pcs);
    gdd_pcs_ind = flipud(gdd_pcs_ind);
    [ gdd_pcs_ind gdd_pcs ]
    str_pcs = gdd_pcs_ind(gdd_pcs>0.1);
    oth_pcs = gdd_pcs_ind(gdd_pcs>0.0);
    
%% Investigate overall 
% Summed, original, Top 4
men_pca_com_tp4 = sum(pca_com(str_pcs,:));

[ ~ , max_ind_tp4 ] = sort(men_pca_com_tp4);
    max_ind_tp4 = fliplr(max_ind_tp4);

[ lbl_nme(max_ind_tp4(1:15)) num2cell(men_pca_com_tp4(max_ind_tp4(1:15)))' ]

% Summed, original, All
men_pca_com_tp9 = sum(pca_com(oth_pcs,:));

[ ~ , max_ind_tp9 ] = sort(men_pca_com_tp9);
    max_ind_tp9 = fliplr(max_ind_tp9);

[ lbl_nme(max_ind_tp9(1:15)) num2cell(men_pca_com_tp9(max_ind_tp9(1:15)))' ]

%% Laterality
% Summed, original, Top 4
lhs_cnt = 0;
rhs_cnt = 0;
bil_cnt = 0;

lbl_nme_hld = lbl_nme(max_ind_tp4(1:27));

for iRW = 1:size(lbl_nme_hld,1)
    nme_hld = regexp(lbl_nme_hld{iRW,1},'==','split');
    if strcmpi(nme_hld{1}(1:2),'rh') && strcmpi(nme_hld{2}(1:2),'rh')
        rhs_cnt = rhs_cnt + 1;
    elseif strcmpi(nme_hld{1}(1:2),'lh') && strcmpi(nme_hld{2}(1:2),'lh')
        lhs_cnt = lhs_cnt + 1;
    else
        bil_cnt = bil_cnt + 1;
    end
end

[ lhs_cnt rhs_cnt bil_cnt ];

[ lhs_cnt/27 bil_cnt/27 rhs_cnt/27 ]

% Summed, original, All
lhs_cnt = 0;
rhs_cnt = 0;
bil_cnt = 0;

lbl_nme_hld = lbl_nme(max_ind_tp9(1:27));

for iRW = 1:size(lbl_nme_hld,1)
    nme_hld = regexp(lbl_nme_hld{iRW,1},'==','split');
    if strcmpi(nme_hld{1}(1:2),'rh') && strcmpi(nme_hld{2}(1:2),'rh')
        rhs_cnt = rhs_cnt + 1;
    elseif strcmpi(nme_hld{1}(1:2),'lh') && strcmpi(nme_hld{2}(1:2),'lh')
        lhs_cnt = lhs_cnt + 1;
    else
        bil_cnt = bil_cnt + 1;
    end
end

[ lhs_cnt rhs_cnt bil_cnt ];

[ lhs_cnt/27 bil_cnt/27 rhs_cnt/27 ]

%% Location
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

% Summed, original, Top 4
lbl_nme_hld = lbl_nme(max_ind_tp4(1:27));

lat_tmp_cnt = 0;
med_tmp_cnt = 0;
lat_frn_cnt = 0;
occ_ctx_cnt = 0;
lat_par_cnt = 0;
sup_ctx_cnt = 0;

oth_cnt = 0;
oth_hld = {};

for iRW = 1:size(lbl_nme_hld,1)
    
    nme_hld = regexp(lbl_nme_hld{iRW,1},'==','split');
    
    if  any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end)))
        lat_tmp_cnt = lat_tmp_cnt + sum( [ any(strcmpi(lat_tmp,nme_hld{1}(4:end))) any(strcmpi(lat_tmp,nme_hld{2}(4:end))) ] );
    end
    if  any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end)))
        med_tmp_cnt = med_tmp_cnt +  sum( [ any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) any(strcmpi(ven_med_tmp,nme_hld{2}(4:end))) ] );
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

tot_num = sum([ lat_tmp_cnt med_tmp_cnt lat_frn_cnt occ_ctx_cnt lat_par_cnt sup_ctx_cnt oth_cnt ]);
[ lat_tmp_cnt/tot_num med_tmp_cnt/tot_num lat_frn_cnt/tot_num occ_ctx_cnt/tot_num lat_par_cnt/tot_num sup_ctx_cnt/tot_num oth_cnt/tot_num ]

% Summed, original, Top 9
lbl_nme_hld = lbl_nme(max_ind_tp9(1:27));

lat_tmp_cnt = 0;
med_tmp_cnt = 0;
lat_frn_cnt = 0;
occ_ctx_cnt = 0;
lat_par_cnt = 0;
sup_ctx_cnt = 0;

oth_cnt = 0;
oth_hld = {};

for iRW = 1:size(lbl_nme_hld,1)
    
    nme_hld = regexp(lbl_nme_hld{iRW,1},'==','split');
    
    if  any(strcmpi(lat_tmp,nme_hld{1}(4:end)))     || any(strcmpi(lat_tmp,nme_hld{2}(4:end)))
        lat_tmp_cnt = lat_tmp_cnt + 1;
    end
    if  any(strcmpi(ven_med_tmp,nme_hld{1}(4:end))) || any(strcmpi(ven_med_tmp,nme_hld{2}(4:end)))
        med_tmp_cnt = med_tmp_cnt + 1;
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

tot_num = sum([ lat_tmp_cnt med_tmp_cnt lat_frn_cnt occ_ctx_cnt lat_par_cnt sup_ctx_cnt oth_cnt ]);
[ lat_tmp_cnt/tot_num med_tmp_cnt/tot_num lat_frn_cnt/tot_num occ_ctx_cnt/tot_num lat_par_cnt/tot_num sup_ctx_cnt/tot_num oth_cnt/tot_num ]

%% Individual PC Plots - Original Approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best - 1
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Number1Connection';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(1))    };
fcfg.edg_wgh = { repmat(3,1,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

% Best - 2
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Number2Connection';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(2))    };
fcfg.edg_wgh = { repmat(3,1,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

% Best - 3
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Number3Connection';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(3))    };
fcfg.edg_wgh = { repmat(3,1,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Top 3
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Top3Connections';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(1:3))    };
fcfg.edg_wgh = { repmat(3,3,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

% Top 10
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Top10Connections';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(1:10))    };
fcfg.edg_wgh = { repmat(3,10,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

% Top 15
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Top15Connections';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(1:15))    };
fcfg.edg_wgh = { repmat(3,15,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

% Top 27
fcfg = [];

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/output/figure7/';
fcfg.out_nme = 'Top27Connections';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_nme(max_ind(1:27))    };
fcfg.edg_wgh = { repmat(3,27,1) };
fcfg.edg_col = { rgb('royal purple')      };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)

%% Individual PC Plots - Leo Approach
connectome_link_visualization_ejk( '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/example_data/fc011_roi.nii' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/example_data/phn_cnn_table.xlsx' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/tractogram_data/rMNI152_T1_1mm_brain.nii' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/tractogram_data/HCP842_100k.trk' )
                               
connectome_link_visualization_ejk( '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/example_data/aal.nii' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/example_data/example_table.xlsx' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/tractogram_data/rMNI152_T1_1mm_brain.nii' , ...
                                   '/home/ekaestne/gitrep/MMIL/OUTSIDE/LeoConnectomeVisualization/tractogram_data/HCP842_100k.trk' )
                               
                               
                               
                               
                               
                               
                               
                               
                               

