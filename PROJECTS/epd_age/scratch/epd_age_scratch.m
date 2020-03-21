clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/SUBJECTS/projects/epd_age/';
prj_nme = 'Epilepsy_and_Aging';

%% Category of Subjects
% Orig Dataset
grp_adn = mmil_readtext([ prj_dir '/' 'ADNI_VisitInfo_cleaned.csv' ]); % mmil_readtext([ prj_dir '/' 'ADNI_VisitInfo.csv' ]);
grp_adn = grp_adn(:,[ 1 2 5 ]); % grp_adn(:,[ 2 3 7 ]);
grp_adn{1,3} = 'Diagnosis';
grp_adn(:,3) = strrep( grp_adn(:,3) , 'Cluster-Derived Normal' , 'CDN' );
grp_adn(:,3) = strrep( grp_adn(:,3) , 'Normal Control' , 'NormalControl' );
mci_ind = 1:size(grp_adn,1); mci_ind(strcmpi(grp_adn(:,3),'NormalControl')) = []; mci_ind(1) = [];
con_ind = 1:size(grp_adn,1); con_ind(mci_ind) = []; con_ind(1) = [];
grp_adn(:,4) = grp_adn(:,3);
grp_adn{1,3} = 'Diagnosis'; grp_adn{1,4} = 'Phenotype';
grp_adn(mci_ind,3) = {'MCI_OldAge'};
grp_adn(con_ind,3) = {'HC_OldAge'};
grp_adn = strrep(grp_adn,'NormalControl','HC_OldAge');

grp_age = mmil_readtext([ prj_dir '/' 'Matched_dataset.csv' ]);
grp_age = grp_age(:,[ 3 4 2 ]);
grp_age = grp_age(1:159,:);

grp_aus = mmil_readtext([ prj_dir '/' 'ucsd_data_082019_cleaned.csv' ]);
grp_aus = strrep( grp_aus , 'EPD' , 'EPD_OldAge' );
grp_aus(:,4) = grp_aus(:,3);

unique(grp_adn(:,3))
unique(grp_age(:,3))

grp_epd = mmil_readtext([ prj_dir '/' 'subjects_phenotype_update2.csv' ]);
grp_epd( cellfun(@isempty,grp_epd(:,3)) , :) = [];
grp_epd = grp_epd( : , [1 2 6 3] );
grp_epd = strrep( grp_epd , 'EPD' , 'EPD_MiddleAge' );
grp_epd = strrep( grp_epd , 'HC' , 'HC_MiddleAge' );
grp_epd = strrep( grp_epd , 'Language & memory' , 'Language_memory' );
grp_epd = strrep( grp_epd , 'No impairment'     , 'Noimpairment' );

%
con_inc = intersect( grp_age( strcmpi(grp_age(:,3),'HC') ,1) , grp_adn(2:end,1) );
con_mss = setxor( con_inc , grp_age( strcmpi(grp_age(:,3),'HC') ,1) );

mci_inc = intersect( grp_age( strcmpi(grp_age(:,3),'MCI') ,1) , grp_adn(2:end,1) );
mci_mss = setxor( mci_inc , grp_age( strcmpi(grp_age(:,3),'MCI') ,1) );

[ ~ , ~ , adn_ind ] = intersect( mci_inc , grp_adn(2:end,1) );
tabulate( grp_adn(adn_ind,4) )

grp_nme = cell(0);
for iS = 2:size( grp_age , 1 )
    grp_nme = { grp_nme{:} grp_adn{ strcmpi( grp_adn(:,1) , grp_age{iS,1} ) , 4} };
end
tabulate(grp_nme)

%% Category of Subjects - Updated to 80
% Older EPD - (80 EPD)
grp_aus = mmil_readtext([ prj_dir '/' 'ucsd_data_090319_cleaned.csv' ]);
grp_aus = strrep( grp_aus , 'EPD' , 'EPD_OldAge' );
grp_aus(:,4) = grp_aus(:,3);

epd_non_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
epd_usd_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf';

loc_aus = repmat( {epd_non_dir} , size(grp_aus(2:end,:),1) , 1);
loc_aus(string_find( grp_aus(:,1) , 'epd0' )) = repmat( {epd_usd_dir} , numel(string_find( grp_aus(:,1) , 'epd0' )) , 1);

% Middle Age EPD - ( 69 EPD , 46 HC )
grp_epd = mmil_readtext([ prj_dir '/' 'subjects_phenotype_update2_young.csv' ]);
grp_epd = grp_epd(2:end,:);
grp_epd( cellfun(@isempty,grp_epd(:,3)) , :) = [];
grp_epd = grp_epd( : , [1 2 6 3] );
grp_epd = strrep( grp_epd , 'EPD' , 'EPD_MiddleAge' );
grp_epd = strrep( grp_epd , 'HC' , 'HC_MiddleAge' );
grp_epd = strrep( grp_epd , 'Language & memory' , 'Language_memory' );
grp_epd = strrep( grp_epd , 'No impairment'     , 'Noimpairment' );

epd_usd_dir = '/home/mmilmcd/data/FSRECONS'; % '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf';
loc_epd = repmat( {epd_usd_dir} , size(grp_epd(2:end,:),1) , 1);

% Load MCI data
grp_adn = mmil_readtext([ prj_dir '/' 'ADNI_VisitInfo_cleaned.csv' ]); % mmil_readtext([ prj_dir '/' 'ADNI_VisitInfo.csv' ]);
grp_adn = grp_adn(:,[ 1 2 5 3]); % grp_adn(:,[ 2 3 7 ]);
grp_adn{1,3} = 'Diagnosis';
grp_adn(:,3) = strrep( grp_adn(:,3) , 'Cluster-Derived Normal' , 'CDN' );
grp_adn(:,3) = strrep( grp_adn(:,3) , 'Normal Control' , 'NormalControl' );
mci_ind = 1:size(grp_adn,1); mci_ind(strcmpi(grp_adn(:,3),'NormalControl')) = []; mci_ind(1) = [];
con_ind = 1:size(grp_adn,1); con_ind(mci_ind) = []; con_ind(1) = [];
grp_adn(:,5) = grp_adn(:,3);
grp_adn{1,3} = 'Diagnosis'; grp_adn{1,5} = 'Phenotype';
grp_adn(mci_ind,3) = {'MCI_OldAge'};
grp_adn(con_ind,3) = {'HC_OldAge'};
grp_adn(:,[3 5]) = strrep(grp_adn(:,[3 5]),'NormalControl','HC_OldAge');
grp_adn = grp_adn( : , [ 1 2 3 5 4]);

% MCI-Total - ( 80 MCI )
pos_ind = [ string_find( grp_adn(:,4) , 'Amnestic' ) ; string_find( grp_adn(:,4) , 'Dysexecutive' ) ; string_find( grp_adn(:,4) , 'Dysnomic' ) ];
[ ~ , age_srt ] = sort(cell2mat(grp_adn(pos_ind,5)));

grp_mci_tot = grp_adn(pos_ind(age_srt(1:80)),1:4);
grp_mci_tot(:,4) = strrep(grp_mci_tot(:,4) , 'Amnestic' , '' );
grp_mci_tot(:,4) = strrep(grp_mci_tot(:,4) , 'Dysexecutive' , '' );
grp_mci_tot(:,4) = strrep(grp_mci_tot(:,4) , 'Dysnomic' , '' );

adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
loc_mci_tot = repmat( {adn_dir} , size(grp_mci_tot,1) , 1);

mci_tot_age = mean(cell2mat(grp_adn(pos_ind(age_srt(1:80)),5)));

% MCI-Amnestic - ( 80 MCI )
pos_ind = string_find( grp_adn(:,4) , 'Amnestic' );
[ ~ , age_srt ] = sort(cell2mat(grp_adn(pos_ind,5)));

grp_mci_amn = grp_adn(pos_ind(age_srt(1:80)),1:4);
grp_mci_amn(:,3) = strrep(grp_mci_amn(:,3) , 'MCI_OldAge' , '' );

adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
loc_mci_amn = repmat( {adn_dir} , size(grp_mci_amn,1) , 1);

mci_amn_age = mean(cell2mat(grp_adn(pos_ind(age_srt(1:80)),5)));

% MCI-Dysnomic - ( XX MCI )
pos_ind = string_find( grp_adn(:,4) , 'Dysnomic' );
[ ~ , age_srt ] = sort(cell2mat(grp_adn(pos_ind,5)));

grp_mci_dys = grp_adn(pos_ind(age_srt(1:80)),1:4);
grp_mci_dys(:,3) = strrep(grp_mci_dys(:,3) , 'MCI_OldAge' , '' );

adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
loc_mci_dys = repmat( {adn_dir} , size(grp_mci_dys,1) , 1);

mci_dys_age = mean(cell2mat(grp_adn(pos_ind(age_srt(1:80)),5)));

% MCI-CDN - ( XX MCI )
pos_ind = string_find( grp_adn(:,4) , 'CDN' );
[ ~ , age_srt ] = sort(cell2mat(grp_adn(pos_ind,5)));

grp_mci_cdn = grp_adn(pos_ind(age_srt(1:80)),1:4);
grp_mci_cdn(:,3) = strrep(grp_mci_cdn(:,3) , 'MCI_OldAge' , '' );

adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
loc_mci_cdn = repmat( {adn_dir} , size(grp_mci_cdn,1) , 1);

mci_cdn_age = mean(cell2mat(grp_adn(pos_ind(age_srt(1:80)),5)));

% Older HC - ( 80 HC )
pos_ind = string_find( grp_adn(:,4) , 'HC_OldAge' );
[ ~ , age_srt ] = sort(cell2mat(grp_adn(pos_ind,5)));

grp_mci_con = grp_adn(pos_ind(age_srt(1:80)),1:4);

adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
loc_mci_con = repmat( {adn_dir} , size(grp_mci_con,1) , 1);

mci_con_age = mean(cell2mat(grp_adn(pos_ind(age_srt(1:80)),5)));

% Put together
dta_use = [ grp_aus ; ...
            grp_epd ; ...
            grp_mci_tot ; ...
            grp_mci_amn ; ...
            grp_mci_dys ; ...
            grp_mci_cdn ; ...
            grp_mci_con ];

dta_loc = [ loc_aus ; ...
            loc_epd ; ...
            loc_mci_tot ; ...
            loc_mci_amn ; ...
            loc_mci_dys ; ...
            loc_mci_cdn ; ...
            loc_mci_con ];

%% Load Surface
prj_dir = '/home/ekaestne/PROJECTS/';

hms     = {'lhs' 'rhs'};

%
sbj_nme = dta_use; %[ grp_adn [ {'Phenotype'} ; repmat({''},size(grp_adn,1)-1,1)] ; grp_aus(2:end,:) repmat({''},size(grp_aus(2:end,:),1),1) ; grp_epd(2:end,:) ];

%
% epd_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
% adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
% mid_dir = '/home/mmilmcdRSI/data/fsurf';
mri_fsr_dir = dta_loc; % [ repmat( {adn_dir} , size(grp_adn,1) , 1)  ; repmat( {epd_dir} , size(grp_aus(2:end,:),1) , 1) ; repmat( {mid_dir} , size(grp_epd(2:end,:),1) , 1) ];

%
old_dta.lhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs.mat'  ]);
old_dta.rhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs.mat'  ]);

tic;
fcfg = [];

fcfg.prj_dir = prj_dir;

fcfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 2819; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

fcfg.anl_dir = 'analysis';

for iH = 1:numel(hms)
    
    srf_dta  = nan(size(sbj_nme,1)-1,163842);
    fcfg.hms = hms{iH}(1:2);
    
    for iS = 2:size(sbj_nme,1)
        
%         if any(strcmpi( old_dta.([hms{iH} '_' 'dta']).srf_dta_sbj , sbj_nme{iS,1} )) && ~strcmpi(mri_fsr_dir{iS},'/home/mmilmcd/data/FSRECONS')
%             
%             sbj_ind = find(strcmpi( old_dta.([hms{iH} '_' 'dta']).srf_dta_sbj , sbj_nme{iS,1} ));
%             
%             srf_dta(iS-1,:) = old_dta.([hms{iH} '_' 'dta']).srf_dta( sbj_ind , : );
%             
%         else
            
            fcfg.prc_dir     = mri_fsr_dir{iS-1};
            fcfg.sbj_fsr_dir = sbj_nme{iS,2};
            
            srf_dta_hld = ejk_extract_vertices(fcfg);
            if ~isempty(srf_dta_hld); srf_dta(iS-1,:) = srf_dta_hld; end
            
%         end
%         
    end
    srf_dta_sbj = sbj_nme(2:end,1);
    save([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% Load ROIs
prj_dir = '/home/ekaestne/PROJECTS/';

hms     = {'lhs' 'rhs'};

%
sbj_nme = dta_use; %[ grp_adn [ {'Phenotype'} ; repmat({''},size(grp_adn,1)-1,1)] ; grp_aus(2:end,:) repmat({''},size(grp_aus(2:end,:),1),1) ; grp_epd(2:end,:) ];

%
% epd_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
% adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
% mid_dir = '/home/mmilmcdRSI/data/fsurf';
mri_fsr_dir = dta_loc; % [ repmat( {adn_dir} , size(grp_adn,1) , 1)  ; repmat( {epd_dir} , size(grp_aus(2:end,:),1) , 1) ; repmat( {mid_dir} , size(grp_epd(2:end,:),1) , 1) ];

prc_nme = { '' '.a2009s' }; % '' '.a2009s'

tic;

fcfg = [];

fcfg.prj_dir = prj_dir;

fcfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 2819; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

fcfg.anl_dir = 'analysis';

for iPR = 1:numel(prc_nme)
    
    fcfg.prc_nme = prc_nme{iPR};
    
    for iS = 2:size(sbj_nme,1)
        
        fcfg.ovr_dir     = mri_fsr_dir{iS-1};
        fcfg.sbj_fsr_dir = sbj_nme{iS,2};
        
        [ gry_thk_dta_hld , tot_lbl ]= ejk_extract_grey_thickness(fcfg);
        try
            if ~isempty(gry_thk_dta_hld)
                gry_thk_dta(iS-1,:) = gry_thk_dta_hld;
                fprintf( [ sbj_nme{iS,2} ' : ' 'Data Loaded\n' ] )
            else
                gry_thk_dta(iS-1,:) = nan(1,size(gry_thk_dta,2));
                fprintf( [ sbj_nme{iS,2} ' : ' 'Missing\n' ] )
            end
        catch
            gry_thk_dta(iS-1,:) = nan(1,size(gry_thk_dta,2));
            fprintf( [ sbj_nme{iS,2} ' : ' 'Missing\n' ] )
        end
    end
    
    sve_sbj_nme = sbj_nme(2:end,1);
    cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI' '_' 'thickness' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(gry_thk_dta,1))] [tot_lbl ; num2cell(gry_thk_dta)] ]);
    clear gry_thk_dta tot_lbl
    
end
toc

%% Load Volumes
prj_dir = '/home/ekaestne/PROJECTS/';

%
sbj_nme = dta_use; %[ grp_adn [ {'Phenotype'} ; repmat({''},size(grp_adn,1)-1,1)] ; grp_aus(2:end,:) repmat({''},size(grp_aus(2:end,:),1),1) ; grp_epd(2:end,:) ];

%
% epd_dir = '/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf';
% adn_dir = '/home/mmilmcdRSI/data_ADNI/fsurf/';
% mid_dir = '/home/mmilmcdRSI/data/fsurf';
mri_fsr_dir = dta_loc; % [ repmat( {adn_dir} , size(grp_adn,1) , 1)  ; repmat( {epd_dir} , size(grp_aus(2:end,:),1) , 1) ; repmat( {mid_dir} , size(grp_epd(2:end,:),1) , 1) ];

%
tic;
fcfg = [];

fcfg.prj_dir = prj_dir;

for iS = 2:size(sbj_nme,1)
    
    fcfg.ovr_dir     = mri_fsr_dir{iS-1};
    fcfg.sbj_fsr_dir = sbj_nme{iS,2};
    
    [ vol_dta_hld , tot_lbl] = ejk_extract_volumes(fcfg);
    if ~isempty(vol_dta_hld); vol_dta(iS-1,:) = vol_dta_hld; else vol_dta(iS-1,:) = nan(1,size(vol_dta,2)); end
    
end

sve_sbj_nme = sbj_nme(2:end,1);
cell2csv([prj_dir '/' 'OUTPUT' '/' prj_nme '/'  'Volumes' '_' 'aparc' '.csv'],[ ['SubjId' ; sve_sbj_nme] [tot_lbl ; num2cell(vol_dta)] ]);
clear vol_dta

toc

%% ROI Plots
prj_dir = '/home/ekaestne/PROJECTS/';

sbj_nme = dta_use; %[ grp_adn ; grp_aus(2:end,:) ; grp_epd(2:end,:) ];

%
dkn_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI_thickness_aparc_.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = dkn_dta;
dkn_dta = ejk_load_mcd_data(fcfg);

dst_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI_thickness_aparc_xa2009s.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = dst_dta;
dst_dta = ejk_load_mcd_data(fcfg);

vol_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Volumes_aparc.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = vol_dta;
vol_dta = ejk_load_mcd_data(fcfg);

% BAR
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'volumes' 'desikan' 'destrieaux' }; % 'sbj_cln_cnt'
fcfg.ydt     = {  vol_dta   dkn_dta   dst_dta     }; % sbj_cln_cnt

fcfg.grp     = sbj_nme;
fcfg.grp_col = [ 3 ...
                 4 ...
                 nan ...
                 3 ...
                 3 ...
                 4 ...
                 3 ...
                 4 ...
                 nan ];
fcfg.nme_col = [ sbj_nme(1,3) ...
                 sbj_nme(1,4) ...
                 {''} ...
                 sbj_nme(1,3) ...
                 sbj_nme(1,4) ...
                 {''} ...
                 sbj_nme(1,3) ...
                 sbj_nme(1,4) ...
                 {''} ];

fcfg.grp_nme = { { 'HC_OldAge' 'MCI_OldAge' } ...
                 { 'HC_OldAge' 'CDN' 'Amnestic' 'Dysnomic' }  ...
                 {} ...
                 { 'HC_OldAge' 'EPD_OldAge' } ...
                 { 'EPD_OldAge' 'MCI_OldAge' } ...
                 { 'EPD_OldAge' 'CDN' 'Amnestic'   'Dysnomic' } ...
                 { 'HC_MiddleAge' 'EPD_MiddleAge' } ...
                 { 'HC_MiddleAge' 'Noimpairment' 'Language' 'Memory' 'Language_memory' } ...
                 {} };
fcfg.grp_clr = { { 'black' 'blue' } ...
                 { 'black' 'bluish grey' 'light blue' 'dark blue' } ...
                 {} ...
                 { 'black' 'orange' } ...
                 { 'orange' 'blue' } ...
                 { 'orange' 'bluish grey' 'light blue' 'dark blue' } ...
                 { 'black' 'green' } ...
                 { 'black' 'greenish grey' 'light green' 'dark green' 'greenish yellow' } ...
                 { } };
fcfg.xdt     = { [ 1 2 ] ...
                 [ 1 3 4 5  ] ...
                 [] ...
                 [ 1 2 ] ...
                 [ 1 2 ] ...
                 [ 1 3 4 5  ] ...
                 [ 1 2 ] ...
                 [ 1 3 4 5 6  ] ...
                 [] };

fcfg.plt_cmp = { { [1 2] } ...
                 { [1 2] [1 3] [1 4] } ...
                 {} ...
                 { [1 2] } ...
                 { [1 2] } ...
                 { [1 2] [1 3] [1 4] } ...
                 { [1 2] } ...
                 { [1 2] [1 3] [1 4] [1 5]} ...
                 {} };
fcfg.plt_anv = { [] ...
                 [ 1 2 3 4 ] ...
                 [] ...
                 [] ...
                 [] ...
                 [ 1 2 3 4 ] ...
                 [] ...
                 [ 1 2 3 4 5 ] ...
                 [] };

ejk_roi_bar(fcfg)

% SCATTER 
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'volumes' 'desikan' 'destrieaux' };
fcfg.ydt     = {  vol_dta   dkn_dta   dst_dta     };

fcfg.grp     = sbj_nme;
fcfg.grp_col = [ 3 ...
                 4 ...
                 nan ...
                 3 ...
                 3 ...
                 4 ...
                 3 ...
                 4 ...
                 nan ];
fcfg.grp_nme = { { 'HC_OldAge' 'MCI_OldAge' } ...
                 { 'HC_OldAge' 'CDN' 'Amnestic' 'Dysnomic' }  ...
                 {} ...
                 { 'HC_OldAge' 'EPD_OldAge' } ...
                 { 'EPD_OldAge' 'MCI_OldAge' } ...
                 { 'EPD_OldAge' 'CDN' 'Amnestic'   'Dysnomic' } ...
                 { 'HC_MiddleAge' 'EPD_MiddleAge' } ...
                 { 'HC_MiddleAge' 'Noimpairment' 'Language' 'Memory' 'Language_memory' } ...
                 {} };
fcfg.xdt     = { [ 1 2 ] ...
                 [ 1 3 4 5  ] ...
                 [] ...
                 [ 1 2 ] ...
                 [ 1 2 ] ...
                 [ 1 3 4 5  ] ...
                 [ 1 2 ] ...
                 [ 1 3 4 5 6  ] ...
                 [] };
fcfg.grp_clr = { { 'black' 'blue' } ...
                 { 'black' 'bluish grey' 'light blue' 'dark blue' } ...
                 {} ...
                 { 'black' 'orange' } ...
                 { 'orange' 'bluish grey' 'light blue' 'dark blue' } ...
                 {} ...
                 { 'black' 'green' } ...
                 { 'black' 'greenish grey' 'light green' 'dark green' 'greenish yellow' } ...
                 { } };
ejk_roi_scatter(fcfg)

% %% Volume Plots
% prj_dir = '/home/ekaestne/PROJECTS/';
% 
% sbj_nme = [ grp_adn ; grp_aus(2:end,:) ];
% 
% %
% vol_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Volumes_aparc.csv' ];
% fcfg = [];
% fcfg.sbj_nme = sbj_nme;
% fcfg.fle_nme = vol_dta;
% vol_dta = ejk_load_mcd_data(fcfg);
% 
% % BAR
% fcfg = [];
% 
% fcfg.prj_dir = prj_dir;
% fcfg.prj_nme = prj_nme;
% 
% fcfg.dta_lbl = { 'volumes' }; % 'sbj_cln_cnt'
% fcfg.ydt     = {  vol_dta  }; % sbj_cln_cnt
% 
% fcfg.grp     = sbj_nme;
% fcfg.grp_col = [ 3                                3                                                                                      3 ];
% fcfg.nme_col = sbj_nme(1,fcfg.grp_col);
% 
% fcfg.grp_nme = { { 'NormalControl'    'EPD'}     { 'NormalControl' 'Cluster-Derived Normal' 'Amnestic'   'Dysnomic'  'Dysexecutive' }  { 'EPD'    'Cluster-Derived Normal' 'Amnestic'   'Dysnomic'  'Dysexecutive' } };
% fcfg.grp_clr = { { 'black'             'orange' } { 'black'          'bluish grey'            'light blue' 'dark blue' 'bluish purple' } { 'orange' 'bluish grey'            'light blue' 'dark blue' 'bluish purple' }    };
% fcfg.xdt     = { [ 1                   2 ]        [ 1                3                        4            5           6  ]              [ 1                3                        4            5           6  ] };
% 
% fcfg.plt_cmp = { { [1 2] }                        { [1 2] [1 3] [1 4] [1 5]}                                                             { [1 2] [1 3] [1 4] [1 5]} };
% fcfg.plt_anv = { []                               [ 1 2 3 4 5 ]                                                                          [ 1 2 3 4 5 ] };
% 
% ejk_roi_bar(fcfg)
% 
% % SCATTER 
% fcfg = [];
% 
% fcfg.prj_dir = prj_dir;
% fcfg.prj_nme = prj_nme;
% 
% fcfg.dta_lbl = { 'volumes' }; % 'sbj_cln_cnt'
% fcfg.ydt     = {  vol_dta  }; % sbj_cln_cnt
% 
% fcfg.grp     = sbj_nme;
% fcfg.grp_col = [ 3                               3                                                                                     3 ];
% fcfg.grp_nme = { { 'NormalControl'    'EPD'}     { 'NormalControl' 'Cluster-Derived Normal' 'Amnestic'   'Dysnomic'  'Dysexecutive' }  { 'EPD'    'Cluster-Derived Normal' 'Amnestic'   'Dysnomic'  'Dysexecutive' } };
% fcfg.xdt     = { [ 1                   2 ]       [ 1                3                        4            5           6  ]             [ 1         3                        4            5           6  ]  };
% fcfg.grp_clr = { { 'black'             'orange' } { 'black'          'bluish grey'            'light blue' 'dark blue' 'bluish purple' } { 'orange' 'bluish grey'            'light blue' 'dark blue' 'bluish purple' }    };
% 
% ejk_roi_scatter(fcfg)

%% Surface Map Plots
%
sbj_nme = dta_use;
sbj_nme{1,3} = 'Diagnosis';
sbj_nme{1,4} = 'Phenotype';

%
lhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs.mat'  ]);
rhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs.mat'  ]);

%
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.sbj_grp     = sbj_nme;
fcfg.sbj_grp_col = { 'Diagnosis' ...
                     'Diagnosis' ...
                     'Diagnosis' ...
                     'Phenotype' ...
                     'Phenotype' };
fcfg.sbj_grp_nme = { { 'HC_OldAge'    'HC_MiddleAge' 'EPD_MiddleAge'} ...
                     { 'HC_OldAge'    'MCI_OldAge' 'EPD_OldAge' } ...
                     { 'HC_MiddleAge' 'EPD_MiddleAge+' } ...
                     { 'HC_OldAge' 'CDN' 'Amnestic'   'Dysnomic' } ...
                     { 'HC_MiddleAge' 'Noimpairment' 'Language' 'Memory' 'Language_memory' } };
fcfg.sbj_grp_cmp = { { [2 1] [3 1] } ...
                     { [2 1] [3 1] } ...
                     { [2 1] } ...
                     { [2 1] [3 1] [4 1] } ...
                     { [2 1] [3 1] [4 1] [5 1] } };

fcfg.dta_lbl = { 'corticalthickness' ...
                 'corticalthickness' ...
                 'corticalthickness' ...
                 'corticalthickness' ...
                 'corticalthickness' };
fcfg.dta     = { lhs_dta rhs_dta }; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

fcfg.hms     = {'lhs' 'rhs'};

ejk_surf_group_avg(fcfg)

%% Category of Subjects - Anny
% Orig Dataset
grp_any = mmil_readtext([ prj_dir '/' 'Matched_dataset.csv' ]); % mmil_readtext([ prj_dir '/' 'ADNI_VisitInfo.csv' ]);
grp_any = grp_any(1:159,[4 3 2]);
grp_any{1,3} = 'Diagnosis';

%
lhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_lhs.mat'  ]);
rhs_dta = load([ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging' '/' 'surf_aMRI_thickness_rhs.mat'  ]);

%
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.sbj_grp     = grp_any;
fcfg.sbj_grp_col = { 'Diagnosis' };
fcfg.sbj_grp_nme = { {'HC' 'MCI' } };
fcfg.sbj_grp_cmp = { { [2 1] } };

fcfg.dta_lbl = { 'corticalthickness' };
fcfg.dta     = { lhs_dta rhs_dta }; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

fcfg.hms     = {'lhs' 'rhs'};

ejk_surf_group_avg(fcfg)

%% 

%% ROI Plots for Presentation
prj_dir = '/home/ekaestne/PROJECTS/';

sbj_nme = dta_use; %[ grp_adn ; grp_aus(2:end,:) ; grp_epd(2:end,:) ];

%
dkn_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'MRI_thickness_aparc_.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = dkn_dta;
dkn_dta = ejk_load_mcd_data(fcfg);

vol_dta = [prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'Volumes_aparc.csv' ];
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.fle_nme = vol_dta;
vol_dta = ejk_load_mcd_data(fcfg);

% BAR
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'volumes_presentation' 'desikan_presentation' }; % 'sbj_cln_cnt'
fcfg.ydt     = {  vol_dta                dkn_dta  }; % sbj_cln_cnt

fcfg.grp     = sbj_nme;
fcfg.grp_col = [ 3 ...
                 3 ...
                 3 ...
                 nan nan nan nan nan nan ];
fcfg.nme_col = [ sbj_nme(1,3) ...
                 sbj_nme(1,3) ...
                 sbj_nme(1,3) ...
                 {''} {''} {''} {''} {''} {''} ];

fcfg.grp_nme = { { 'HC_OldAge' 'EPD_OldAge' 'MCI_OldAge' } { 'HC_OldAge' 'MCI_OldAge' } { 'HC_MiddleAge' 'EPD_MiddleAge' } ...
                 {} {} {} {} {} {} };
fcfg.grp_clr = { { 'black' 'orange' 'blue' }               { 'black' 'dark blue' }      { 'light grey' 'light orange' } ...
                 {''} {''} {''} {''} {''} {''} };
fcfg.xdt     = { [ 1 2 3 ] [ 1 2 ] [ 1 2 ] ...
                 [] [] [] [] [] [] };

fcfg.plt_cmp = { { [1 2] [2 3] } { [1 2] } { [1 2] } ...
                 {} {} {} {} {} {} };
fcfg.plt_anv = { [1 2 3] [] [] ...
                 [] [] [] [] [] [] };

ejk_roi_bar(fcfg)

% SCATTER 
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;

fcfg.dta_lbl = { 'volumes_presentation' 'desikan_presentation' }; % 'sbj_cln_cnt'
fcfg.ydt     = {  vol_dta                dkn_dta  }; % sbj_cln_cnt

fcfg.grp     = sbj_nme;
fcfg.grp_col = [ 3 ...
                 3 ...
                 3 ...
                 nan nan nan nan nan nan ];
fcfg.grp_nme = { { 'HC_OldAge' 'EPD_OldAge' 'MCI_OldAge' } { 'HC_OldAge' 'MCI_OldAge' } { 'HC_MiddleAge' 'EPD_MiddleAge' } ...
                 {} {} {} {} {} {} };
fcfg.xdt     = { [ 1 2 3 ] [ 1 2 ] [ 1 2 ] ...
                 [] [] [] [] [] [] };
fcfg.grp_clr = { { 'black' 'orange' 'blue' }               { 'black' 'dark blue' }      { 'light grey' 'light orange' } ...
                 {''} {''} {''} {''} {''} {''} };
             
ejk_roi_scatter(fcfg)




