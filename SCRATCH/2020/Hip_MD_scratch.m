clear; clc;

%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'scratch';

%
red_fle = 'sbj000_total_2019_03_27.csv';

%
dta_dir_dti = '/home/mmilmcdRSI/MetaData/MCD_RSI/ROI_Summaries';
dta_nme_dti = 'DTI_all.csv';

dta_dir_t1d = '/home/mmilmcd/MetaData/MRI/ROI_Summaries';
dta_nme_t1d = 'MRI_all.csv';

dfa_mes_nme = { 'aseg_FA-Left-Hippocampus' ...
                'aseg_FA-Right-Hippocampus' ...
                ...
                'aseg_FA-Left-Amygdala' ...
                'aseg_FA-Right-Amygdala' ...                
                ...
                'wmparc_FA-wm-lh-entorhinal' ...
                'wmparc_FA-wm-rh-entorhinal' ...
                ...
                'wmparc_FA-wm-lh-parahippocampal' ...
                'wmparc_FA-wm-rh-parahippocampal' ...
                };

dmd_mes_nme = { 'aseg_MD-Left-Hippocampus' ...
                'aseg_MD-Right-Hippocampus' ...
                ...
                'aseg_MD-Left-Amygdala' ...
                'aseg_MD-Right-Amygdala' ...
                ...
                'wmparc_MD-wm-lh-entorhinal' ...
                'wmparc_MD-wm-rh-entorhinal' ...
                ...
                'wmparc_MD-wm-lh-parahippocampal' ...
                'wmparc_MD-wm-rh-parahippocampal' ...
                };

gry_mes_nme = { 'subcort_vol-Left-Hippocampus' ...
                'subcort_vol-Right-Hippocampus' ...
                ...
               'subcort_vol-Left-Amygdala' ...
                'subcort_vol-Right-Amygdala' ...
                ...
                'cort_thick-ctx-lh-entorhinal' ...
                'cort_thick-ctx-rh-entorhinal' ...
                ...
                'cort_thick-ctx-lh-parahippocampal' ...
                'cort_thick-ctx-rh-parahippocampal' ...
                };
            
pre_tst_nme = { 'log_mem_nor_scr_one' ...
                'log_mem_nor_scr_two' ...
                'cvl_lfr_nor_scr' ...
                'vp1_nor_scr' ...
                'cvl_tot_nor_scr' ...
                'bnt_nor_scr' ...
                'cat_flu_nor_scr' ...
                };
            
pst_tst_nme = { 'log_mem_scr_one' ...
                'log_mem_scr_two' ...
                'cvl_lfr_scr' ...
                'vp1_scr' ...
                'cvl_tot_scr' ...
                'bnt_scr' ...
                'cat_flu_scr' ...
                };
        
%% Load Data
% DTI Data
dti_dta = mmil_readtext( [ dta_dir_dti '/' dta_nme_dti ] );

dfa_mtx = [ dti_dta(:,1:2) cell( size(dti_dta,1) , numel(dfa_mes_nme)) ];
for iC = 1:numel(dfa_mes_nme)
    dfa_mtx(:,iC+2) = dti_dta( : , strcmpi( dti_dta(1,:) , dfa_mes_nme{iC} ) );
end

dmd_mtx = [ dti_dta(:,1:2) cell( size(dti_dta,1) , numel(dmd_mes_nme)) ];
for iC = 1:numel(dmd_mes_nme)
    dmd_mtx(:,iC+2) = dti_dta( : , strcmpi( dti_dta(1,:) , dmd_mes_nme{iC} ) );
end

% DTI Data
mri_dta = mmil_readtext( [ dta_dir_t1d '/' dta_nme_t1d ] );

gry_mtx = [ mri_dta(:,1:2) cell( size(mri_dta,1) , numel(gry_mes_nme)) ];
for iC = 1:numel(gry_mes_nme)
    gry_mtx(:,iC+2) = mri_dta( : , strcmpi( mri_dta(1,:) , gry_mes_nme{iC} ) );
end

% Cognitive Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);

cfg = [];
[ pst_cog_dta_hld , ~ ] = ejk_post_cognitive(cfg,sbj_cog,sbj_scn);

% Check DTI rows
fsr_fld = cell( numel( sbj_cog.sbj_nme ) , 3 );
for iS = 1:numel( sbj_cog.sbj_nme )
    fsr_fld{ iS , 1 } = sbj_cog.sbj_nme{iS};
    fsr_fld{ iS , 2 } = [  dfa_mtx( strcmpi( dfa_mtx(:,1) , sbj_cog.sbj_nme{iS} ) , 2 ) ];
    fsr_fld{ iS , 3 } = 1;    
end

fsr_mri_fld = cell( numel( sbj_cog.sbj_nme ) , 3 );
for iS = 1:numel( sbj_cog.sbj_nme )
    fsr_mri_fld{ iS , 1 } = sbj_cog.sbj_nme{iS};
    fsr_mri_fld{ iS , 2 } = [  gry_mtx( strcmpi( gry_mtx(:,1) , sbj_cog.sbj_nme{iS} ) , 2 ) ];
    fsr_mri_fld{ iS , 3 } = 1;    
end

fsr_mri_fld( cellfun( @isempty , fsr_fld(:,2) ) , : ) = [];
fsr_fld( cellfun( @isempty , fsr_fld(:,2) ) , : ) = [];
ttt = find( cellfun( @numel , fsr_fld(:,2) ) > 1 );

fsr_fld{ttt(1),3}  = 1; 
fsr_fld{ttt(2),3}  = 1;
fsr_fld{ttt(3),3}  = 1;
fsr_fld{ttt(4),3}  = 3;
fsr_fld{ttt(5),3}  = 2;
fsr_fld{ttt(6),3}  = 3;
fsr_fld{ttt(7),3}  = 3;
fsr_fld{ttt(8),3}  = 1;
fsr_fld{ttt(9),3}  = 1;
fsr_fld{ttt(10),3} = 1;
fsr_fld{ttt(11),3} = 1;
fsr_fld{ttt(12),3} = 1;
fsr_fld{ttt(13),3} = 2;
fsr_fld{ttt(14),3} = 2;
fsr_fld{ttt(15),3} = 1;
fsr_fld{ttt(16),3} = 1;
fsr_fld{ttt(17),3} = 1;
fsr_fld{ttt(18),3} = 1;
fsr_fld{ttt(19),3} = 2;
fsr_fld{ttt(20),3} = 2;
fsr_fld{ttt(21),3} = 2;
fsr_fld{ttt(22),3} = 2;
fsr_fld{ttt(23),3} = 2;
fsr_fld{ttt(24),3} = 2;
fsr_fld{ttt(25),3} = 2;
fsr_fld{ttt(26),3} = 1;
fsr_fld{ttt(27),3} = 1;
fsr_fld{ttt(28),3} = 1;
fsr_fld{ttt(29),3} = 1;
fsr_fld{ttt(30),3} = 1;
fsr_fld{ttt(31),3} = 2;
fsr_fld{ttt(32),3} = 1;
fsr_fld{ttt(33),3} = 2;
fsr_fld{ttt(34),3} = 2;
fsr_fld{ttt(35),3} = 3;
fsr_fld{ttt(36),3} = 1;
fsr_fld{ttt(37),3} = 1;
fsr_fld{ttt(38),3} = 2;
fsr_fld{ttt(39),3} = 1;
fsr_fld{ttt(40),3} = 2;
fsr_fld{ttt(41),3} = 1;
fsr_fld{ttt(42),3} = 1;
fsr_fld{ttt(43),3} = 1;
fsr_fld{ttt(44),3} = 1;
fsr_fld{ttt(45),3} = 1;
fsr_fld{ttt(46),3} = 1;
fsr_fld{ttt(47),3} = 2;
fsr_fld{ttt(48),3} = 2;
fsr_fld{ttt(49),3} = 2;
fsr_fld{ttt(50),3} = 1;
fsr_fld{ttt(51),3} = 2;
fsr_fld{ttt(52),3} = 2;

%% Make Data
img_dfa_dta     = cell( size(fsr_fld,1) + 1 , numel(dfa_mes_nme) + 1 );
img_dmd_dta     = cell( size(fsr_fld,1) + 1 , numel(dmd_mes_nme) + 1 );
img_vol_dta     = cell( size(fsr_fld,1) + 1 , numel(gry_mes_nme) + 1 );
pre_neu_psy_dta = cell( size(fsr_fld,1) + 1 , numel(pre_tst_nme) + 1 );
pst_neu_psy_dta = cell( size(fsr_fld,1) + 1 , numel(pst_tst_nme) + 1 );
lat_dta         = cell( size(fsr_fld,1) + 1 , 3 );

for iS = [ 1:106 108:size(fsr_fld,1)-2 ]
    
    % DTI FA Imaging Data
    img_dfa_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    for iFA = 1:numel(dfa_mes_nme)
        img_dfa_dta{ iS+1 , iFA+1 } = dfa_mtx{ strcmpi( dfa_mtx(:,2) , fsr_fld{iS,2}{fsr_fld{iS,3}} ) , strcmpi( dfa_mtx(1,:) , dfa_mes_nme{iFA} ) };
    end    
    
    % DTI MD Imaging Data
    img_dmd_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    for iMD = 1:numel(dmd_mes_nme)
        img_dmd_dta{ iS+1 , iMD+1 } = dmd_mtx{ strcmpi( dmd_mtx(:,2) , fsr_fld{iS,2}{fsr_fld{iS,3}} ) , strcmpi( dmd_mtx(1,:) , dmd_mes_nme{iMD} ) };
    end   
    
    % Grey Matter Data
    img_vol_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    for iMD = 1:numel(gry_mes_nme)
        if sum(strcmpi( gry_mtx(:,2) , fsr_fld{iS,2}{fsr_fld{iS,3}} ))~=0
            img_vol_dta{ iS+1 , iMD+1 } = gry_mtx{ strcmpi( gry_mtx(:,2) , fsr_fld{iS,2}{fsr_fld{iS,3}} ) , strcmpi( gry_mtx(1,:) , gry_mes_nme{iMD} ) };
        else
            if numel(fsr_mri_fld{iS,2})==1
                img_vol_dta{ iS+1 , iMD+1 } = gry_mtx{ strcmpi( gry_mtx(:,2) , fsr_mri_fld{iS,2}{fsr_mri_fld{iS,3}} ) , strcmpi( gry_mtx(1,:) , gry_mes_nme{iMD} ) };
            else
               error('problem')
            end
        end
    end
    
    % Neuropsych Data
    pre_neu_psy_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    for iNP = 1:numel(pre_tst_nme)
        pre_neu_psy_dta{ iS+1 , iNP+1 } = sbj_cog.( pre_tst_nme{iNP} )( strcmpi( sbj_cog.sbj_nme , fsr_fld{iS,1} ) );
    end
          
    % Post Neuropsych Data
    pst_neu_psy_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    for iNP = 1:numel(pst_tst_nme)
        pst_neu_psy_dta{ iS+1 , iNP+1 } = pst_cog_dta_hld.( pst_tst_nme{iNP} )( strcmpi( pst_cog_dta_hld.sbj_nme , fsr_fld{iS,1} ) );
    end

    % Laterality
    lat_dta{ iS+1 , 1 } = fsr_fld{iS,1};
    hld_chk = strfind( fsr_fld{iS,1} , 'fc' );
    if isempty( hld_chk )
        lat_dta{ iS+1 , 2 } = 'EPD';
        lat_dta{ iS+1 , 3 } = sbj_sze.( 'sbj_sde_ons' ){ strcmpi( sbj_sze.sbj_nme , fsr_fld{iS,1} ) };
    else
        lat_dta{ iS+1 , 2 } = 'HC';
        lat_dta{ iS+1 , 3 } = 'HC';
    end
    
end

img_dfa_dta(1,:)     = { 'SubjID' dfa_mes_nme{:} };
img_dmd_dta(1,:)     = { 'SubjID' dmd_mes_nme{:} };
img_vol_dta(1,:)     = { 'SubjID' gry_mes_nme{:} };
pre_neu_psy_dta(1,:) = { 'SubjID' pre_tst_nme{:} };
pst_neu_psy_dta(1,:) = { 'SubjID' pst_tst_nme{:}};
lat_dta(1,:)         = { 'SubjID' 'Laterality' 'Group' };

%% Format Data
grp_dta = lat_dta([ 2:107 109:size(fsr_fld,1)-1 ] , :);

% DTI FA Imaging Data
dfa_dta.sbj_nme = grp_dta(:,1);
for iC = 1:numel(dfa_mes_nme)
    dfa_dta.(['x' mmil_spec_char(dfa_mes_nme{iC},{'-'})]) = cell2mat(img_dfa_dta( 2:end ,iC + 1));
end

% DTI MD Imaging Data
dmd_dta.sbj_nme = grp_dta(:,1);
for iC = 1:numel(dmd_mes_nme)
    dmd_dta.(['x' mmil_spec_char(dmd_mes_nme{iC},{'-'})]) = cell2mat(img_dmd_dta( 2:end ,iC + 1));
end

% Grey Matter Imaging Data
gry_dta.sbj_nme = grp_dta(:,1);
for iC = 1:numel(gry_mes_nme)
    gry_dta.(['x' mmil_spec_char(gry_mes_nme{iC},{'-'})]) = cell2mat(img_vol_dta( 2:end ,iC + 1));
end

% Pre-Neuropsych Data
pre_cog_dta.sbj_nme = grp_dta(:,1);
for iC = 1:numel(pre_tst_nme)
    pre_cog_dta.(['x' mmil_spec_char(pre_tst_nme{iC},{'-'})]) = cell2mat(pre_neu_psy_dta( 2:end ,iC + 1));
end

% Post-Neuropsych Data
pst_cog_dta.sbj_nme = grp_dta(:,1);
for iC = 1:numel(pst_tst_nme)
    pst_cog_dta.(['x' mmil_spec_char(pst_tst_nme{iC},{'-'})]) = cell2mat(pst_neu_psy_dta( 2:end ,iC + 1));
end

%% Make Plots
% Pre-Operative Scores, MD
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'gry_post'  'gry_pre'   'wmd_post'  'wmd_pre'   'wfa_post'  'wfa_pre' };
fcfg.xdt     = { pst_cog_dta pre_cog_dta pst_cog_dta pre_cog_dta pst_cog_dta pre_cog_dta }; % Secondary Folders
fcfg.ydt     = { gry_dta     gry_dta     dmd_dta     dmd_dta     dfa_dta     dfa_dta }; % Master Folders

fcfg.grp     = grp_dta;
fcfg.grp_col = [ 1                    2 ];
fcfg.grp_nme = { { 'EPD'    'HC' }    { 'L'    'R' } };
fcfg.grp_clr = { { 'orange' 'green' } { 'blue' 'red' } };

ejk_roi_scatter(fcfg)














