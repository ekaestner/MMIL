clear; clc
 
cov_rmv_ind = [ 173 207 212 218 222 ];
srf_rmv_ind = cov_rmv_ind-1; 

out_dta_loc = '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/';

%%
cov_dta     = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/epd_age_covariates_sm313_no_nan.csv');
    cov_dta(cov_rmv_ind,1)

cov_dta(cov_rmv_ind,:) = [];

cell2csv( [ out_dta_loc '/' 'epd_age_covariates_no_nan_no_mass.csv' ] , cov_dta)

% lhs_176_dta
load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_sm176_no_nan.mat');
    srf_dta_sbj(srf_rmv_ind)
    
srf_dta_sbj(srf_rmv_ind) = [];    
srf_dta(srf_rmv_ind,:)   = [];   

save( [ out_dta_loc '/' 'surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat' ] , 'srf_dta_sbj' , 'srf_dta' )

% rhs_176_dta
load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_rhs_sm176_no_nan.mat');
    srf_dta_sbj(srf_rmv_ind)
    
srf_dta_sbj(srf_rmv_ind) = [];    
srf_dta(srf_rmv_ind,:)   = [];   

save( [ out_dta_loc '/' 'surf_aMRI_thickness_rhs_sm176_no_nan_no_mass.mat' ] , 'srf_dta_sbj' , 'srf_dta' )

% lhs_313_dta
load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_sm313_no_nan.mat');
    srf_dta_sbj(srf_rmv_ind)
    
srf_dta_sbj(srf_rmv_ind) = [];    
srf_dta(srf_rmv_ind,:)   = [];   

save( [ out_dta_loc '/' 'surf_aMRI_thickness_lhs_sm313_no_nan_no_mass.mat' ] , 'srf_dta_sbj' , 'srf_dta' )

% rhs_313_dta
load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_rhs_sm313_no_nan.mat');
    srf_dta_sbj(srf_rmv_ind)
    
srf_dta_sbj(srf_rmv_ind) = [];    
srf_dta(srf_rmv_ind,:)   = [];   

save( [ out_dta_loc '/' 'surf_aMRI_thickness_rhs_sm313_no_nan_no_mass.mat' ] , 'srf_dta_sbj' , 'srf_dta' )

%% Check Which Subjects Are Removed
cov_rmv_ind = [ 173 207 212 218 222 ];
srf_rmv_ind = cov_rmv_ind-1; 

lft_org_srf_dta = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_sm176_no_nan.mat');
    % lft_org_srf_dta.srf_dta_sbj(srf_rmv_ind)

cov_dta     = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/epd_age_covariates_sm313_no_nan.csv');    
    ely_ind = find(strcmpi( cov_dta(:,5) , 'Early' ))-1;
        ely_rmv_ind = intersect( srf_rmv_ind , ely_ind );
        ely_kep_ind = setxor( ely_ind , ely_rmv_ind );
    lte_ind = find(strcmpi( cov_dta(:,5) , 'Late' ))-1;
        lte_rmv_ind = intersect( srf_rmv_ind , lte_ind );
        lte_kep_ind = setxor( lte_ind , lte_rmv_ind );
        
[col_loc,albl,actbl] = fs_read_annotation(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' 'lh' '.aparc.split.annot']);

mid_pre_gry = find(col_loc==find(strcmpi(albl,'middle-precentral')));
ent_rhn_gry = find(col_loc==find(strcmpi(albl,'entorhinal')));

% middle precentral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ely_kep_thk = mean( lft_org_srf_dta.srf_dta(ely_kep_ind,mid_pre_gry) , 2);
ely_rmv_thk = mean( lft_org_srf_dta.srf_dta(ely_rmv_ind,mid_pre_gry) , 2);

lte_kep_thk = mean( lft_org_srf_dta.srf_dta(lte_kep_ind,mid_pre_gry) , 2);
lte_rmv_thk = mean( lft_org_srf_dta.srf_dta(lte_rmv_ind,mid_pre_gry) , 2);

%
fcfg = [];

fcfg.ydt     = { ely_kep_thk  ely_rmv_thk lte_kep_thk  lte_rmv_thk };
fcfg.xdt     = { 1            1           2            2           };

fcfg.fce_col = { rgb('dark grey') rgb('red')  rgb('dark grey') rgb('red') };
fcfg.edg_col = { rgb('black')     rgb('red')  rgb('black')     rgb('black') };

fcfg.xlb = { 'Early' 'Late' };
fcfg.ylb = 'Thickness';

fcfg.ylm = [1 4];

ejk_scatter(fcfg)

%
men_ely_all_thk = mean(mean( lft_org_srf_dta.srf_dta(ely_ind,mid_pre_gry) , 2));
men_ely_rmv_thk = mean(mean( lft_org_srf_dta.srf_dta(ely_kep_ind,mid_pre_gry) , 2));

men_lte_all_thk = mean(mean( lft_org_srf_dta.srf_dta(lte_ind,mid_pre_gry) , 2));
men_lte_rmv_thk = mean(mean( lft_org_srf_dta.srf_dta(lte_kep_ind,mid_pre_gry) , 2));

% entorhinal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ely_kep_thk = mean( lft_org_srf_dta.srf_dta(ely_kep_ind,ent_rhn_gry) , 2);
ely_rmv_thk = mean( lft_org_srf_dta.srf_dta(ely_rmv_ind,ent_rhn_gry) , 2);

lte_kep_thk = mean( lft_org_srf_dta.srf_dta(lte_kep_ind,ent_rhn_gry) , 2);
lte_rmv_thk = mean( lft_org_srf_dta.srf_dta(lte_rmv_ind,ent_rhn_gry) , 2);

%
fcfg = [];

fcfg.ydt     = { ely_kep_thk  ely_rmv_thk lte_kep_thk  lte_rmv_thk };
fcfg.xdt     = { 1            1           2            2           };

fcfg.fce_col = { rgb('dark grey') rgb('red')  rgb('dark grey') rgb('red') };
fcfg.edg_col = { rgb('black')     rgb('red')  rgb('black')     rgb('black') };

fcfg.xlb = { 'Early' 'Late' };
fcfg.ylb = 'Thickness';

fcfg.ylm = [1 4];

ejk_scatter(fcfg)

%
men_ely_all_thk = mean(mean( lft_org_srf_dta.srf_dta(ely_ind,ent_rhn_gry) , 2));
men_ely_rmv_thk = mean(mean( lft_org_srf_dta.srf_dta(ely_kep_ind,ent_rhn_gry) , 2));

men_lte_all_thk = mean(mean( lft_org_srf_dta.srf_dta(lte_ind,ent_rhn_gry) , 2));
men_lte_rmv_thk = mean(mean( lft_org_srf_dta.srf_dta(lte_kep_ind,ent_rhn_gry) , 2));






