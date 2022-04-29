load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

% Load Wada/Handedness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = cog_tst_inc
    sbj_nme_hld = [ sbj_nme_hld ; sbj_cog.sbj_nme(~isnan(sbj_cog.(cog_tst_nme{iT}))) ];
end
sbj_nme = unique(sbj_nme_hld);

clear sbj_cog

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Data/Wada.csv';
[ sbj_dem, ~ , ~ , ~, ~, sbj_srg] = mmil_load_redcap(fcfg);

sbj_wda_hld = [ sbj_dem.sbj_nme sbj_srg.wda_lng cell( numel(sbj_srg.wda_lng), 1) sbj_dem.sbj_hnd ];

% Load fMRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fMRI_laterality.csv'];
fcfg.dta_col = 3;
[ fmr_aln_dta, fmr_aln_dta_sbj, fmr_aln_dta_col] = ejk_dta_frm( fcfg );

for iS = 1:size(sbj_wda_hld,1)
    sbj_ind = find(strcmpi(fmr_aln_dta_sbj,sbj_wda_hld{iS,1}));
    if isempty(sbj_ind)
       sbj_wda_hld{iS,3} = '';
    else
        sbj_wda_hld{iS,3} = fmr_aln_dta{sbj_ind,end};
    end    
end

% Categorize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_dom = cell(size(sbj_wda_hld,1),2);
for iS = 1:size(sbj_wda_hld,1)
   if ~isempty(sbj_wda_hld{iS,2})
       if strcmpi(sbj_wda_hld{iS,2},'left')
           sbj_dom{iS,1} = 'T';
       else
           sbj_dom{iS,1} = 'A';
       end
       sbj_dom{iS,2} = 'Wada'; 
   elseif ~isempty(sbj_wda_hld{iS,3})
       sbj_dom{iS,1} = sbj_wda_hld{iS,3};
       sbj_dom{iS,2} = 'fMRI'; 
   elseif  ~isempty(sbj_wda_hld{iS,4})
       if strcmpi(sbj_wda_hld{iS,4},'R')
           sbj_dom{iS,1} = 'T';
       else
           sbj_dom{iS,1} = 'A';
       end 
       sbj_dom{iS,2} = 'Hand'; 
   else
       sprintf('%s : missing\n',sbj_wda_hld{iS,1})
   end   
end

sbj_wda_hld = [ sbj_wda_hld sbj_dom ];

% cell2csv('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Data/LanguageDominance.csv', sbj_wda_hld)

%% Load Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
% fcfg.all_num = 1;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical_wdominance.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'subcort_vol_ICV_cor_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ mri_dta, mri_dta_sbj, mri_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ fib_dta, fib_dta_sbj, fib_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_QC.csv'];
fcfg.dta_col = 5;
% fcfg.all_num = 1;
[ wmp_dta, wmp_dta_sbj, wmp_dta_col] = ejk_dta_frm( fcfg );

prd_dta_sbj = [ cog_dta_sbj ];

prd_dta_col = [ cog_dta_col cln_dta_col mri_dta_col fib_dta_col wmp_dta_col {'lang_dom'} ];

prd_dta     = [ cog_dta     cln_dta     mri_dta     fib_dta     wmp_dta     sbj_wda_hld(:,end-1)];

%% Categorize
pre_grp_use = { 'controls_pre_3T_allSurg_all' 'tle_pre_3T_allSurg_left' 'tle_pre_3T_allSurg_right' };

pst_grp_use = { 'tle_post_3T_ATLonly_left' 'tle_post_3T_ATLonly_right' };

cat_cut = -1.28;

% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pre_dom_tbl = cell(3,4);
pst_dom_tbl = cell(9,3);

cog_cat_bnt_pst = repmat({''},numel(cog_dta(:,3)),1);
cog_cat_bnt_pst( cell2mat(cog_dta(:,3)) <= cat_cut) = {'Impaired'};
cog_cat_bnt_pst( cell2mat(cog_dta(:,3)) > cat_cut) = {'NoChange'};

cog_cat_ant_pst = repmat({''},numel(cog_dta(:,4)),1);
cog_cat_ant_pst( cell2mat(cog_dta(:,4)) <= cat_cut) = {'Impaired'};
cog_cat_ant_pst( cell2mat(cog_dta(:,4)) > cat_cut) = {'NoChange'};

% Pre-operative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(pre_grp_use)
    
    % Total
    tot_hld = tabulate( sbj_wda_hld( grp.(pre_grp_use{iG}), end-1) );
    pre_dom_tbl{ 2, iG+1} = [num2str(tot_hld{strcmpi(tot_hld(:,1),'T'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'A'),2})];
    
    % Type
    tot_hld = tabulate( sbj_wda_hld( grp.(pre_grp_use{iG}), end) );
    try pre_dom_tbl{ 3, iG+1} = [num2str(tot_hld{strcmpi(tot_hld(:,1),'Wada'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'fMRI'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'Hand'),2})];
    catch pre_dom_tbl{ 3, iG+1} = ['0' ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'fMRI'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'Hand'),2})]; end
    
end

pre_dom_tbl{1,2} = 'HC'; pre_dom_tbl{1,3} = 'L-TLE'; pre_dom_tbl{1,4} = 'R-TLE';
pre_dom_tbl{2,1} = 'Dominance (Typical/Atypical)'; pre_dom_tbl{3,1} = 'Type (Wada/fMRI/Hand)';

% cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/scratch/lang_lat/presurgical_table.csv',pre_dom_tbl);

% Post-operative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(pst_grp_use)

% Total
tot_hld = tabulate( sbj_wda_hld( grp.(pst_grp_use{iG}), end-1) );
pst_dom_tbl{ 2, iG+1} = [num2str(tot_hld{strcmpi(tot_hld(:,1),'T'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'A'),2})];

% Type
tot_hld = tabulate( sbj_wda_hld( grp.(pst_grp_use{iG}), end) );
try pst_dom_tbl{ 3, iG+1} = [num2str(tot_hld{strcmpi(tot_hld(:,1),'Wada'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'fMRI'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'Hand'),2})];
catch pst_dom_tbl{ 3, iG+1} = ['0' ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'fMRI'),2}) ' / ' num2str(tot_hld{strcmpi(tot_hld(:,1),'Hand'),2})]; end

% Dom-x-Cognition
[ crs_tbl, ~, ~, crs_lbl ] = crosstab( sbj_wda_hld( grp.(pst_grp_use{iG}), end-1), cog_cat_bnt_pst(grp.(pst_grp_use{iG})));
pst_dom_tbl{ 5, iG+1} = [num2str(crs_tbl(strcmpi(crs_lbl(:,1),'T'),strcmpi(crs_lbl(:,2),'NoChange'))) ' / ' num2str(crs_tbl(strcmpi(crs_lbl(:,1),'T'),strcmpi(crs_lbl(:,2),'Impaired'))) ];
pst_dom_tbl{ 6, iG+1} = [num2str(crs_tbl(strcmpi(crs_lbl(:,1),'A'),strcmpi(crs_lbl(:,2),'NoChange'))) ' / ' num2str(crs_tbl(strcmpi(crs_lbl(:,1),'A'),strcmpi(crs_lbl(:,2),'Impaired'))) ];

[ crs_tbl, ~, ~, crs_lbl ] = crosstab( sbj_wda_hld( grp.(pst_grp_use{iG}), end-1), cog_cat_ant_pst(grp.(pst_grp_use{iG})));
pst_dom_tbl{ 8, iG+1} = [num2str(crs_tbl(strcmpi(crs_lbl(:,1),'T'),strcmpi(crs_lbl(:,2),'NoChange'))) ' / ' num2str(crs_tbl(strcmpi(crs_lbl(:,1),'T'),strcmpi(crs_lbl(:,2),'Impaired'))) ];
pst_dom_tbl{ 9, iG+1} = [num2str(crs_tbl(strcmpi(crs_lbl(:,1),'A'),strcmpi(crs_lbl(:,2),'NoChange'))) ' / ' num2str(crs_tbl(strcmpi(crs_lbl(:,1),'A'),strcmpi(crs_lbl(:,2),'Impaired'))) ];

end

pst_dom_tbl{1,2} = 'L-TLE'; pst_dom_tbl{1,3} = 'R-TLE';
pst_dom_tbl{2,1} = 'Dominance (Typical/Atypical)'; pst_dom_tbl{3,1} = 'Type (Wada/fMRI/Hand)';
pst_dom_tbl{5,1} = 'BNT: Typical (Decline/None)'; pst_dom_tbl{6,1} = 'BNT: Atypical (Decline/None)';
pst_dom_tbl{8,1} = 'ANT: Typical (Decline/None)'; pst_dom_tbl{9,1} = 'ANT: Atypical (Decline/None)';

% cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/atl_nme/scratch/lang_lat/post-surgical_table.csv',pst_dom_tbl);

%% 
dom_col = find(strcmpi( prd_dta_col, 'lang_dom' ));
typ_dom_num = find( strcmpi(prd_dta( :, dom_col),'T'));
abn_dom_num = find( strcmpi(prd_dta( :, dom_col),'A'));
emp_dom_num = find( cellfun(@isempty,prd_dta( :, dom_col)));

bnt_pre_scr_col = find(strcmpi( prd_dta_col, 'bnt_raw_scr' ));
ant_pre_scr_col = find(strcmpi( prd_dta_col, 'ant_mem_raw_scr' ));

bnt_pst_scr_col = find(strcmpi( prd_dta_col, 'bnt_raw_scr_pst' ));
ant_pst_scr_col = find(strcmpi( prd_dta_col, 'ant_mem_raw_scr_pst' ));

edu_col     = find(strcmpi( prd_dta_col, 'Educ' ));
age_sze_col = find(strcmpi( prd_dta_col, 'AgeOfSeizureOnset' ));
asm_col     = find(strcmpi( prd_dta_col, 'NumAEDs' ));

lft_hip_col = find(strcmpi( prd_dta_col, 'Left_Hippocampus' ));
rgh_hip_col = find(strcmpi( prd_dta_col, 'Right_Hippocampus' ));

lft_ifo_col = find(strcmpi( prd_dta_col, 'L_IFO' ));
rgh_ifo_col = find(strcmpi( prd_dta_col, 'R_IFO' ));

lft_ilf_col = find(strcmpi( prd_dta_col, 'L_ILF' ));
rgh_ilf_col = find(strcmpi( prd_dta_col, 'R_ILF' ));

lft_fus_col = find(strcmpi( prd_dta_col, 'lh_fusiform' ));
rgh_fus_col = find(strcmpi( prd_dta_col, 'rh_fusiform' ));

cog_nme = { { 'BNT'           'ANT' }           { 'BNT'           'ANT' } };
cog_col = { [ bnt_pre_scr_col ant_pre_scr_col ] [ bnt_pst_scr_col ant_pst_scr_col ] };

mes_nme = { { 'Education' 'AgeOfOnset' 'ASM' ...
              'L_Hip'     'R_Hip'     'L_ILF'     'R_ILF'     'L_IFOF'    'R_IFOF'    'L_Fusiform' 'R_Fusiform'  } ...
            { 'BNT_pre' 'ANT_pre' 'Education' 'AgeOfOnset' 'ASM' ...
              'L_Hip'     'R_Hip'     'L_ILF'     'R_ILF'     'L_IFOF'    'R_IFOF'    'L_Fusiform' 'R_Fusiform'  }  };
mes_col = { [ edu_col age_sze_col asm_col ...
              lft_hip_col rgh_hip_col lft_ilf_col rgh_ilf_col lft_ifo_col rgh_ifo_col lft_fus_col  rgh_fus_col ] ...
            [ bnt_pre_scr_col ant_pre_scr_col edu_col age_sze_col asm_col ...
              lft_hip_col rgh_hip_col lft_ilf_col rgh_ilf_col lft_ifo_col rgh_ifo_col lft_fus_col  rgh_fus_col ] };

fdr_grp = { { 1:3 4:11 } { 1:5 6:13 } };
          
%% PLOTTING
plt_grp = { 'tle_post_3T_ATLonly_left' };
cog_typ = [ 2 ];

grp_col = { rgb('Blue') rgb('teal') rgb('dark grey')    };

for iP = 1:numel(plt_grp)
    
    ejk_chk_dir([ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/NoTrendLine' '/' ])
    
    grp_num{1} = intersect( grp.(plt_grp{iP}), typ_dom_num );
    grp_num{2} = intersect( grp.(plt_grp{iP}), abn_dom_num );
    grp_num{3} = intersect( grp.(plt_grp{iP}), emp_dom_num );
    if isempty(grp_num{3}); grp_num(3) = []; end
    
    for iT = 1:numel(cog_nme{cog_typ(iP)})
        for iN = 1:numel(mes_nme{cog_typ(iP)})
            
            % Data Gather
            for iG = 1:numel(grp_num)
                ydt{iG} = cell2mat(prd_dta( grp_num{iG}, mes_col{cog_typ(iP)}(iN) ));
                xdt{iG} = cell2mat(prd_dta( grp_num{iG}, cog_col{cog_typ(iP)}(iT) ));
            end
            
            % Data Plot
            fcfg = [];
            
            fcfg.xdt     = xdt;
            fcfg.ydt     = ydt;
            
            fcfg.fce_col = grp_col;
            fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
            
            fcfg.ylb = { mes_nme{cog_typ(iP)}{iN}  };
            fcfg.xlb = { cog_nme{cog_typ(iP)}{iT} };
            
            fcfg.trd_lne = [ 1 0 ];
            
            fcfg.out_dir = [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/Revisions/NoTrendLine' ];
            fcfg.out_nme = [ cog_nme{cog_typ(iP)}{iT} '_' mes_nme{cog_typ(iP)}{iN} '_left_only' ];
            
            ejk_scatter(fcfg)
            
            clear xdt ydt
            
        end
    end
    
end


