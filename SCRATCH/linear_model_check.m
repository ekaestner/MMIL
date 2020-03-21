clear; clc;

%%
sbj_nme = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/sbj_nme_no_nan.csv');

sbj_cov = mmil_readtext('/home/jxrao/Desktop/Lab/surface_map/files/epd_age_covariates_no_nan.csv');

lhs_dta = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_lhs_no_nan.mat');
rhs_dta = load('/home/jxrao/Desktop/Lab/surface_map/data/mat/surf_aMRI_thickness_rhs_no_nan.mat');

% Total
avg_srf = SurfStatReadSurf({ '/home/mmilmcd/data/FSRECONS/fsaverage/surf/lh.pial' , ...
                             '/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial' });
msk_lhs = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/lh.cortex.label',' '); msk_lhs = cell2mat(msk_lhs(:,1));
    msk_lhs_hld = zeros(1,size(avg_srf.coord,2)/2);
    msk_lhs_hld(msk_lhs+1) = 1;
msk_rhs = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/rh.cortex.label',' '); msk_rhs = cell2mat(msk_rhs(:,1));
    msk_rhs_hld = zeros(1,size(avg_srf.coord,2)/2);
    msk_rhs_hld(msk_rhs+1) = 1; 
avg_msk = [ msk_lhs_hld msk_rhs_hld ];

SurfStatView( avg_msk , avg_srf , 'Masked Surface Average');

% right hemisphere
rhs_avg_srf = SurfStatReadSurf({ '/home/mmilmcd/data/FSRECONS/fsaverage/surf/rh.pial' });
msk_rhs_hld;

% Setup
sex = {'male' 'female'};
str = {'1.5T' '3T'};
                         
% ttt = SurfStatReadSurf({ '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd058_130130_20130130.102308_1/surf/lh.pial' , ...
%                          '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd058_130130_20130130.102308_1/surf/rh.pial' });
% ttt = SurfStatReadData({ '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd058_130130_20130130.102308_1/surf/lh.thickness' , ...
%                          '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd058_130130_20130130.102308_1/surf/rh.thickness' });



%%
grp_ind_one = find( strcmpi( sbj_nme(:,3) , 'EPD_Old') );
grp_ind_two = find( strcmpi( sbj_nme(:,3) , 'HC' ) );
grp_ind = [grp_ind_one ; grp_ind_two]; clear grp_ind_one grp_ind_two

for iS = 1:numel(grp_ind)
    
    dta_ind = find(strcmpi( lhs_dta.srf_dta_sbj , sbj_nme{grp_ind(iS),1} ));
    cov_ind = find(strcmpi( sbj_cov(:,1) , sbj_nme{grp_ind(iS),1} ));
        
    lhs_thk(iS,:) = lhs_dta.srf_dta(dta_ind,:);
    rhs_thk(iS,:) = rhs_dta.srf_dta(dta_ind,:);
    
    grp_cov{iS,1} = sbj_cov{ cov_ind , 6 };
    age_cov(iS,1) = sbj_cov{ cov_ind , 2 };
    str_cov{iS,1} = sex{sbj_cov{ cov_ind , 3 }};
    sex_cov{iS,1} = str{sbj_cov{ cov_ind , 5 }};
    
end

ovr_thk = [ lhs_thk rhs_thk ];

grp_hld = term(grp_cov);
age_hld = term(age_cov);
str_hld = term(str_cov);
sex_hld = term(sex_cov);

%%
mdl_hld = 1 + age_hld + str_hld + sex_hld + grp_hld;
ctr_hld = grp_hld.EPD_Old - grp_hld.HC;

mdl_slm = SurfStatLinMod( ovr_thk , mdl_hld , avg_srf );
grp_stt = SurfStatT( mdl_slm , ctr_hld );

figure(); SurfStatView( grp_stt.t , avg_srf , 't-test?' ); SurfStatColLim([-4.5 4.5]);

rsl_hld = SurfStatResels( grp_stt , avg_msk );
stt_thr = stat_threshold( rsl_hld , length(grp_stt.t) , 1 , grp_stt.df );

[ pvl_hld ~ ~ ] = SurfStatP( grp_stt , logical(avg_msk));
SurfStatView( pvl_hld , avg_srf , 'p-value?')

qvl_hld = SurfStatQ( grp_stt , logical(avg_msk) );
SurfStatView( qvl_hld , avg_srf );

%%
mdl_hld = 1 + age_hld + str_hld + sex_hld + grp_hld;
ctr_hld = grp_hld.EPD_Old - grp_hld.HC;

mdl_slm = SurfStatLinMod( rhs_thk , mdl_hld , rhs_avg_srf );
grp_stt = SurfStatT( mdl_slm , ctr_hld );

figure(); SurfStatView( grp_stt.t , rhs_avg_srf , 't-test?' ); SurfStatColLim([-4.5 4.5]);

rsl_hld = SurfStatResels( grp_stt , logical(msk_rhs_hld));
stt_thr = stat_threshold( rsl_hld , length(grp_stt.t) , 1 , grp_stt.df );

[ pvl_hld , ~ , cls_hld ] = SurfStatP( grp_stt , logical(msk_rhs_hld) , stt_thr );
figure(); SurfStatView( pvl_hld , rhs_avg_srf , 'p-value?');

qvl_hld = SurfStatQ( grp_stt , logical(msk_rhs_hld) );
SurfStatView( qvl_hld , rhs_avg_srf );

%%
mdl_hld = 1 + grp_hld;
ctr_hld = grp_hld.EPD_Old - grp_hld.HC;

mdl_slm = SurfStatLinMod( ovr_thk , mdl_hld , avg_srf );
grp_stt = SurfStatT( mdl_slm , ctr_hld );

figure(); SurfStatView( grp_stt.t , avg_srf , 't-test?' ); SurfStatColLim([-4.5 4.5]);






