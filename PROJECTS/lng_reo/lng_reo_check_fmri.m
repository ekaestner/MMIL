
%%
% %%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%
fmr_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Data' '/' 'Fibers' '/' 'fMRI_aparc_xa2009s_N_FF_xnzvoxels.csv' ] );
fmr_lbl = fmr_dta(1,2:end);
fmr_dta = cell2mat(fmr_dta(2:end,2:end));

fmr_dta_tot = sum( fmr_dta , 2 );

% %%%%%%%%%%%%%%%
lft_ind = find(strcmpi( grp_fle(:,3) , 'L' ));
rgh_ind = find(strcmpi( grp_fle(:,3) , 'R' ));
% lft_ind = 1:5; % find(strcmpi( , ))
con_ind = find(strcmpi( grp_fle(:,3) , 'HC' ));

% Z-Score %%%%%%%%%%%%%%%
fcfg  = [];
fcfg.con_ind = con_ind;
fmr_dta_zsc = ejk_create_zscore(fcfg,fmr_dta);

fcfg  = [];
fcfg.con_ind = con_ind;
fmr_dta_tot_zsc = ejk_create_zscore(fcfg,fmr_dta_tot);

%% Total Activation
% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
lft_tot = fmr_dta_tot( lft_ind , : );
rgh_tot = fmr_dta_tot( rgh_ind , : );
con_tot = fmr_dta_tot( con_ind , : ) ;

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot      rgh_tot      con_tot };

% fcfg.ylm     = [ 0 3000 ];

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total fMRI Activation';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/'];
fcfg.out_nme = 'Destrieaux_OverallActivation';

ejk_scatter(fcfg)

% Z-Score
lft_tot_zsc = fmr_dta_tot_zsc( lft_ind , : );
rgh_tot_zsc = fmr_dta_tot_zsc( rgh_ind , : );
con_tot_zsc = fmr_dta_tot_zsc( con_ind , : );

fcfg = [];

fcfg.xdt     = { 1            2            3 };
fcfg.ydt     = { lft_tot_zsc  rgh_tot_zsc  con_tot_zsc };

fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };

fcfg.xlb = '';
fcfg.ylb = 'Total fMRI Zscore';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/'];
fcfg.out_nme = 'Destrieaux_OverallActivation_ZScore';

ejk_scatter(fcfg)

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = {'Total Activation'};

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/'];
fcfg.out_nme = 'Destrieaux_OverallActivation';

ejk_grp_stt( fcfg , fmr_dta_tot )

%% ROIs
% Raw Data
for iRO = 1:size(fmr_lbl,2)
    
    %
    lft_roi = fmr_dta( lft_ind , iRO );
    rgh_roi = fmr_dta( rgh_ind , iRO );
    con_roi = fmr_dta( con_ind , iRO );
    
    %
    fcfg = [];
    
    fcfg.xdt     = { 1            2            3 };
    fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
    
    fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
    fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
    
    fcfg.xlb = '';
    fcfg.ylb = 'fMRI Activation';
    
    fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/' 'Activation'];
    fcfg.out_nme = mmil_spec_char( fmr_lbl{iRO} , {'-'} );
    
    ejk_scatter(fcfg)
    
end

% Z-Score
for iRO = 1:size(fmr_lbl,2)
    
    %
    lft_roi = fmr_dta_zsc( lft_ind , iRO );
    rgh_roi = fmr_dta_zsc( rgh_ind , iRO );
    con_roi = fmr_dta_zsc( con_ind , iRO );
    
    %
    if ~ (all(isnan(lft_roi)) && all(isnan(rgh_roi)) && all(isnan(con_roi)))
        
        fcfg = [];
        
        fcfg.xdt     = { 1            2            3 };
        fcfg.ydt     = { lft_roi      rgh_roi      con_roi };
        
        fcfg.fce_col = { rgb('blue')  rgb('red')   rgb('greenish grey') };
        fcfg.edg_col = { rgb('black') rgb('black') rgb('black') };
        
        fcfg.xlb = '';
        fcfg.ylb = 'fMRI Zscore';
        
        fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/' 'Activation_zscore'];
        fcfg.out_nme = [ mmil_spec_char( fmr_lbl{iRO} , {'-'} ) '_' 'zscore'];
        
        ejk_scatter(fcfg)
        
    end
    
end

% Stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw Data 
fcfg = [];

fcfg.ind     = { lft_ind rgh_ind con_ind };
fcfg.ind_nme = { 'left'  'right' 'HC' };

fcfg.dta_lbl = fmr_lbl';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Destrieaux' '/'];
fcfg.out_nme = 'Destrieaux_ROIActivation';

ejk_grp_stt( fcfg , fmr_dta )

%% Surface Plots
% prj_dir     = '/home/ekaestne/PROJECTS/';
% prc_dir     = '/home/ekaestne/PROJECTS/DATA/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
% fmr_dir     = '/home/mmilmcd/data/MCD_BOLD/subjects';
% fmr_fsr_dir = '/home/mmilmcd/data/FSRECONS';
% 
% prj_nme = 'LanguageReorganization';
% sbj_grp_nme = 'Lng_Reo_subject_list.csv';
% 
% sbj_grp = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' sbj_grp_nme ]);
% sbj_frs = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv' ]);
% sbj_bld = mmil_readtext([ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv' ]);
% 
% ejk_chk_dir([ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Surface' '/' ])
% 
% % Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f__mri_stt = 'N-FF_GLT#0_Tstat';
% f__mri_nme = 'N-FF';
% f__mri_srf = {'smoothwm' 'pial'};
% 
% fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
% 
% fmr_rng_num = [-7 7];
% cls_sze = 20;
% sig_thr = 2.59;
% 
% f_mri_typ = {'-nzvoxels'};
% 
% % Data Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iS = 23 %1:size( sbj_nme , 1 )
% 
%     sbj_fmr_ind = find(strcmpi( sbj_bld(:,2) , sbj_grp{iS,1} ));
%     sbj_fsr_ind = find(strcmpi( sbj_frs(:,2) , sbj_grp{iS,1} ));
%     
%     % fMRI location
%     if     exist([fmr_dir '/' sbj_bld{sbj_fmr_ind,3} '/' 'orig' '/' f__mri_nme '_' 'cs20+orig.BRIK'])==2
%         org_loc = 'orig';
%         dta_loc = 'orig';
%     elseif exist([fmr_dir '/' sbj_bld{sbj_fmr_ind,3} '/' 'orig.BLOCK' '/' f__mri_nme '_' 'cs20+orig.BRIK'])==2
%         org_loc = 'orig';
%         dta_loc = 'orig.BLOCK';
%     elseif exist([fmr_dir '/' sbj_bld{sbj_fmr_ind,3} '/' 'orig2' '/' f__mri_nme '_' 'cs20+orig.BRIK'])==2
%         org_loc = 'orig';
%         dta_loc = 'orig.BLOCK';    
%     else
%         dta_loc = [];
%     end
% 
%     % Make the Stat Mask
%     fcfg = [];
%     
%     fcfg.prj_dir = prj_dir;
%     fcfg.sbj_nme = sbj_bld{sbj_fmr_ind,2};
%         
%     fcfg.bld_dir = fmr_dir;
%     
%     fcfg.fsr_sbj_nme = sbj_frs{sbj_fsr_ind,3};
%     fcfg.bld_sbj_nme = sbj_bld{sbj_fmr_ind,3};
%     
%     fcfg.f__mri_stt = f__mri_stt;
%     fcfg.f__mri_nme = f__mri_nme;
%     
%     fcfg.cls_sze = cls_sze;
%     fcfg.sig_thr = sig_thr;
%     
%     ejk_fmri_stat_create(fcfg)
%     
%     for iSU = 1:numel(f__mri_srf)
%         
%         %% Pull The Surface Out
%         fcfg = [];
%         
%         fcfg.prj_dir = prj_dir;
%         fcfg.sbj_nme = sbj_bld{sbj_fmr_ind,2};
%         
%         fcfg.fsr_dir = fmr_fsr_dir;
%         fcfg.bld_dir = fmr_dir;
%         
%         fcfg.fsr_sbj_nme = sbj_frs{sbj_fsr_ind,3};
%         fcfg.bld_sbj_nme = sbj_bld{sbj_fmr_ind,3};
%               
%         fcfg.f__mri_stt = f__mri_stt;
%         fcfg.f__mri_nme = f__mri_nme;
%         fcfg.f__mri_srf = f__mri_srf{iSU};
%         
%         fcfg.cls_sze = cls_sze;
%         fcfg.sig_thr = sig_thr;
%         
%         ejk_fmri_surf_extract(fcfg)
%         
%         %% Make Plot
%         fcfg = [];
%         
%         fcfg.prj_dir = prj_dir;
%         fcfg.prj_nme = prj_nme;
%                 
%         fcfg.sbj_nme = sbj_bld{sbj_fmr_ind,2};
%         
%         fcfg.fsr_dir = fmr_fsr_dir;
%         fcfg.bld_dir = fmr_dir;
%         
%         fcfg.fsr_sbj_nme = sbj_frs{sbj_fsr_ind,3};
%         fcfg.bld_sbj_nme = sbj_bld{sbj_fmr_ind,3};
%         
%         fcfg.f__mri_stt = f__mri_stt;
%         fcfg.f__mri_nme = f__mri_nme;
%         fcfg.f__mri_srf = f__mri_srf{iSU};
%         
%         fcfg.fmr_col_map = fmr_col_map;
%         fcfg.fmr_rng_num = fmr_rng_num;
%         
%         fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'DataCheck' '/' 'Fibers' '/' 'Surface' '/' f__mri_srf{iSU} '/' ];
%         fcfg.out_nme = sbj_bld{sbj_fmr_ind,2};
%         
%         mmil_fmri_surf_plot(fcfg)
%         
%     end
%     
% end

























