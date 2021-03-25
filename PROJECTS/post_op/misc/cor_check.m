
%% Original
thk_age = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/ROI/MRI_thickness_aparc__Epilepsy_and_Aging_updated_withTransverse.csv');
srf_age = load('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/NEWDATA/surf_aMRI_thickness_lhs_sm176_no_nan_no_mass.mat');
rvl_age = load( '/home/ekaestner/Downloads/surface_correlations/lhs_fusiform/rvalues_lhs_dep_var_TLE.mat' );

ttt = find(rvl_age.rvalues>0.86);

scatter( cell2mat(thk_age(2:end,7)) , srf_age.srf_dta(:,ttt(1)) )
corrcoef( cell2mat(thk_age(2:end,7)) , srf_age.srf_dta(:,ttt(1)) )

%% Post-op Fusiform ALL
lod_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData/T1.csv');
        lod_dta( cellfun(@(x) strcmpi(x,'NA'),lod_dta) ) = {NaN};
srf_pst = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_aMRI_thickness_lhs_sm176.mat');

non_nan_lod_dta = find(~isnan(cell2mat(lod_dta(3:end,28))))+2;
non_nan_srf_dta = find(~isnan(srf_pst.srf_dta(:,ttt(1))))+2;
non_nan = intersect( non_nan_lod_dta, non_nan_srf_dta);

scatter( cell2mat(lod_dta(non_nan,28)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )
corrcoef( cell2mat(lod_dta(non_nan,28)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )

%% Post-op CVLT Total ALL
lod_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    lod_dta( cellfun(@isempty,lod_dta) ) = {NaN}; 
srf_pst = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_aMRI_thickness_lhs_sm176.mat');

non_nan_lod_dta = find(~isnan(cell2mat(lod_dta(3:end,5))))+2;
non_nan_srf_dta = find(~isnan(srf_pst.srf_dta(:,ttt(1))))+2;
non_nan = intersect( non_nan_lod_dta, non_nan_srf_dta);

scatter( cell2mat(lod_dta(non_nan,5)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )
corrcoef( cell2mat(lod_dta(non_nan,5)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )    

%% Post-op Left/3T only
lod_dta     = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/AllData/T1.csv');
        lod_dta( cellfun(@(x) strcmpi(x,'NA'),lod_dta) ) = {NaN};
grp_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    grp_dta( cellfun(@isempty,grp_dta) ) = {NaN};
srf_pst = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_aMRI_thickness_lhs_sm176.mat');

crr_sbj = find(strcmpi(grp_dta(3:end,end-1),'3_left'))+2;

non_nan_lod_dta = find(~isnan(cell2mat(lod_dta(3:end,28))))+2;
non_nan_srf_dta = find(~isnan(srf_pst.srf_dta(:,ttt(1))))+2;
non_nan = intersect( non_nan_lod_dta, non_nan_srf_dta);
non_nan = intersect( non_nan, crr_sbj);

scatter( cell2mat(lod_dta(non_nan,28)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )
corrcoef( cell2mat(lod_dta(non_nan,28)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )

%% Post-op CVLT Left/3T only
lod_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    lod_dta( cellfun(@isempty,lod_dta) ) = {NaN}; 
grp_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena' '/' 'Memory_FA.csv' ] );
    grp_dta( cellfun(@isempty,grp_dta) ) = {NaN};
srf_pst = load('/home/ekaestner/Downloads/surface_correlations/JohnnyData/surf_aMRI_thickness_lhs_sm176.mat');

crr_sbj = find(strcmpi(grp_dta(3:end,end-1),'3_left'))+2;

non_nan_lod_dta = find(~isnan(cell2mat(lod_dta(3:end,5))))+2;
non_nan_srf_dta = find(~isnan(srf_pst.srf_dta(:,ttt(1))))+2;
non_nan = intersect( non_nan_lod_dta, non_nan_srf_dta);
non_nan = intersect( non_nan, crr_sbj);

scatter( cell2mat(lod_dta(non_nan,5)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )
corrcoef( cell2mat(lod_dta(non_nan,5)) , srf_pst.srf_dta(non_nan-2,ttt(1)) )    





































