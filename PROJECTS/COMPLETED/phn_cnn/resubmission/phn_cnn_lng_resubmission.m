clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/';
red_fle = 'mmil_redcap_2019_10_09.csv';

sbj_nme = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/language_connectome_final_sample.csv');
sbj_nme = sbj_nme(2:end,[1 13]);

sbj_fsr = mmil_readtext('/home/ekaestne/PROJECTS/SUBJECTS/mmilmcdRSI_freesurfer_names.csv');
for iS = 1:size(sbj_nme,1)
    sbj_ind       = find(strcmpi(sbj_fsr(:,2),sbj_nme{iS,1}));
    sbj_nme{iS,3} = sbj_fsr{sbj_ind,3};
end

sbj_nme = sbj_nme(:,[1 3 2]);

sbj_nme(1:33,4)   = strcat('UCSF - ',sbj_nme(1:33,3));
sbj_nme(34:end,4) = strcat('UCSD - ',sbj_nme(34:end,3));

%% FA/MD
fib_typ = { 'FA' 'MD' };
scl_fct = [ 1    1000 ];

tic;
for iFT = 2:numel(fib_typ)
    for iS = 1:size(sbj_nme,1)
        
        fcfg = [];
        
        fcfg.ovr_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';
        
        fcfg.sbj_nme     = sbj_nme{iS,1};
        fcfg.sbj_fsr_dir = sbj_nme{iS,2};
        
        fcfg.atl_dir = 'AtlasTrack';
        fcfg.map_dir = 'fiber_maps';
        fcfg.fib_typ = fib_typ{iFT};
        fcfg.scl_fct = scl_fct(iFT);
        
        fcfg.thr_prb = 0.08;
        
        fcfg.fib_cde = [101 102 103 104 105 106 107 108 109 110 115 116 117 118 119 120 121 122 123 133 134 135 136 137 138 141 142 143 144 145 146 147 148 149 150 1014 1024 2000 2001 2002 2003 2004];
        
        fcfg.min_val = 1e-06;
        
        [ fib_dta_hld , tot_lbl] = ejk_extract_fibers(fcfg);
        if ~isempty(fib_dta_hld); fib_dta(iS,:) = fib_dta_hld; else fib_dta(iS,:) = nan(1,size(fib_dta,2)); end
        
    end
    
    dta_out{iFT} = fib_dta;
    clear fib_dta 
    
end
toc;

cell2csv('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/FA_out.csv' , [ [{'SubjID'};sbj_nme(:,1)] [ tot_lbl ; num2cell(dta_out{1})] ])
cell2csv('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/MD_out.csv' , [ [{'SubjID'};sbj_nme(:,1)] [ tot_lbl ; num2cell(dta_out{2})] ])

%% BAI/BDI
% Load Redcap
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
fcfg.sbj_nme = sbj_nme(:,1);
[~ , ~ , ~ , ~ , sbj_emo] = mmil_load_redcap(fcfg);

% Run ttest between LI/NLI
imp_ind = find(strcmpi(sbj_nme(:,3),'Language Impaired'));
non_ind = find(strcmpi(sbj_nme(:,3),'Not Impaired'));

nanmean(sbj_emo.bai(imp_ind))
nanmean(sbj_emo.bai(non_ind))

nanmean(sbj_emo.bdi(imp_ind))
nanmean(sbj_emo.bdi(non_ind))

[~,pvl_bai] = ttest2( sbj_emo.bai(imp_ind) , sbj_emo.bai(non_ind) );
[~,pvl_bdi] = ttest2( sbj_emo.bdi(imp_ind) , sbj_emo.bdi(non_ind) );

cell2csv('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/BAI_BDI.csv', [ [{'SubjID'};sbj_nme(:,1)] [ {'BAI' 'BDI'} ; num2cell(sbj_emo.bai) num2cell(sbj_emo.bdi)] ] )

% Run ANOVA between groups
usf_imp_ind = find(strcmpi(sbj_nme(:,4),'UCSF -Language Impaired'));
usf_non_ind = find(strcmpi(sbj_nme(:,4),'UCSF -Not Impaired'));
usd_imp_ind = find(strcmpi(sbj_nme(:,4),'UCSD -Language Impaired'));
usd_non_ind = find(strcmpi(sbj_nme(:,4),'UCSD -Not Impaired'));

ttt(1,1) = nanmean(sbj_emo.bai(usd_imp_ind)); ttt(1,2) = nanstd(sbj_emo.bai(usd_imp_ind));
ttt(2,1) = nanmean(sbj_emo.bai(usd_non_ind)); ttt(2,2) = nanstd(sbj_emo.bai(usd_non_ind));
ttt(3,1) = nanmean(sbj_emo.bai(usf_imp_ind)); ttt(3,2) = nanstd(sbj_emo.bai(usf_imp_ind));
ttt(4,1) = nanmean(sbj_emo.bai(usf_non_ind)); ttt(4,2) = nanstd(sbj_emo.bai(usf_non_ind));
ttt
[ pnv(1) , tbl_hld , anv_hld ] = anova1( sbj_emo.bai , sbj_nme(:,4) , 'off' );

ttt(1,1) = nanmean(sbj_emo.bdi(usd_imp_ind)); ttt(1,2) = nanstd(sbj_emo.bdi(usd_imp_ind));
ttt(2,1) = nanmean(sbj_emo.bdi(usd_non_ind)); ttt(2,2) = nanstd(sbj_emo.bdi(usd_non_ind));
ttt(3,1) = nanmean(sbj_emo.bdi(usf_imp_ind)); ttt(3,2) = nanstd(sbj_emo.bdi(usf_imp_ind));
ttt(4,1) = nanmean(sbj_emo.bdi(usf_non_ind)); ttt(4,2) = nanstd(sbj_emo.bdi(usf_non_ind));
ttt
[ pnv(1) , tbl_hld , anv_hld ] = anova1( sbj_emo.bdi , sbj_nme(:,4) , 'off' );

%% Laterality Index
lat_out = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/language_connectome_final_sample.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['/home/ekaestne/PROJECTS/' '/' 'DATA' '/' 'ROIHOLD' '/' 'LI' '/' 'fMRI_aparc_xAlicia_N_FF_xnzvoxels.mat' ]);

lin_dta_hld = squeeze(lin_out(1,:,:)); lin_dta_hld{1,1} = '';
rhs_dta_hld = squeeze(lin_out(2,:,:)); rhs_dta_hld{1,1} = '';
lhs_dta_hld = squeeze(lin_out(3,:,:)); lhs_dta_hld{1,1} = '';

% Parcel Data
for iS = 2:size(lat_out,1)
    
    sbj_ind = find(strcmpi(lin_dta_hld(:,1),lat_out{iS,1}));
    
    sbj_lin_nme{iS-1,1} = lin_dta_hld{sbj_ind,1};
        
    sbj_lin_dta_hld(iS-1,:) = cell2mat(lin_dta_hld(sbj_ind,2:end)); 
    sbj_rhs_dta_hld(iS-1,:) = cell2mat(rhs_dta_hld(sbj_ind,2:end));
    sbj_lhs_dta_hld(iS-1,:) = cell2mat(lhs_dta_hld(sbj_ind,2:end));
    
    if (sbj_lhs_dta_hld(iS-1,:) + sbj_rhs_dta_hld(iS-1,:)) > 0
        sbj_lii(iS-1,:) = (sbj_lhs_dta_hld(iS-1,:) - sbj_rhs_dta_hld(iS-1,:)) ./ (sbj_lhs_dta_hld(iS-1,:) + sbj_rhs_dta_hld(iS-1,:));
    else
        sbj_lii(iS-1,:) = nan(1,4);
    end

    if sbj_lii(iS-1,3) >= 0.2
        sbj_lin_cat{iS-1,1} = 'L';
    elseif sbj_lii(iS-1,3) <= -0.2
        sbj_lin_cat{iS-1,1} = 'R';
    elseif sbj_lii(iS-1,3) < 0.2 && sbj_lii(iS-1,3) > -0.2
        sbj_lin_cat{iS-1,1} = 'B';
    elseif isnan(sbj_lii(iS-1,3))
        sbj_lin_cat{iS-1,1} = 'N/A';
    else
        error('blah blah')
    end
    
end

cell2csv('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/Resubmission/language_connectome_final_sample_with_laterality.csv' , [ lat_out(:,1:end-1) [ 'fMRI Laterality' ; sbj_lin_cat] ])    
    
[ sbj_nme(:,1) num2cell(sbj_lin_dta_hld(:,3)) num2cell(sbj_rhs_dta_hld(:,3)) num2cell(sbj_lhs_dta_hld(:,3)) ]

