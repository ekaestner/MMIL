% load groups
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

% Load Cognitive
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

% Load Clinical
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

% Load Redcap
% Scan-Date | Race | Ethnicity | 
fcfg = [];
fcfg.sbj_nme = cln_dta_sbj;
fcfg.red_fle = red_cap_fle;
fcfg.sep     = '|';
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

% Load Emory
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Emory' '/' 'Emory_SLAH_Memory_11.22.21_cleanNONAMES.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );



%% Fix Groups
grp_nme = { 'ltle_slah' 'ltle_atl' 'rtle_slah' 'rtle_atl' };
for iG = 1:numel(grp_nme)
    grp.surgery.pst_cog_dti_slm.(grp_nme{iG}) = grp.surgery.pst_cog_dti.(grp_nme{iG})(~cellfun(@isnan,cog_dta(grp.surgery.pst_cog_dti.(grp_nme{iG}),strcmpi(cog_dta_col,'lm2_pct'))));
end

%% Combine Emory
emy_dta_col_ind = [numel(emy_dta_col)+1 numel(emy_dta_col)+numel(emy_add_dta_col) ];
emy_dta_col = [ emy_dta_col emy_add_dta_col ];

for iS = 1:numel(emy_add_dta_sbj)
    emy_add_ind = find(strcmpi(emy_dta_sbj,emy_add_dta_sbj{iS}));
    if ~isempty(emy_add_ind)
        emy_dta(emy_add_ind,emy_dta_col_ind(1):emy_dta_col_ind(2)) = emy_add_dta(iS,:);
    end    
end

%% Get UCSD/UCSF data

cln_nui_dta_sbj = cln_dta_sbj;

cln_nui_dta_col = { 'neu_psy_tst_dte' 'neu_psy_tst_dte_pst' 'srg_dte'     'scn_dte' ...
                    'wtr_int_nor_scr' 'wtr_int_nor_scr_pst' 'sbj_sze_frq' 'sbj_mts' ...
                    'eng_out'         'lng_lat'             'lng_lat_mth' 'sbj_eth'     'sbj_rce' 'scn_typ' };

cln_nui_dta = [ sbj_cog.neu_psy_tst_dte sbj_cog.neu_psy_tst_dte_pst sbj_srg.srg_dte sbj_dem.sbj_scn_dte ...
                num2cell(sbj_cog.wtr_int_nor_scr) num2cell(sbj_cog.wtr_int_nor_scr_pst) num2cell(sbj_sze.sbj_sze_frq) sbj_sze.sbj_mts ...
                sbj_srg.eng_out repmat({''},numel(sbj_dem.sbj_nme),1)  repmat({''},numel(sbj_dem.sbj_nme),1) repmat({''},numel(sbj_dem.sbj_nme),1) repmat({''},numel(sbj_dem.sbj_nme),1) repmat({'MMIL'},numel(sbj_dem.sbj_nme),1)];
cln_nui_dta = [ cln_nui_dta ; cell(numel(emy_dta_sbj),size(cln_nui_dta,2))];

% 


%% Add on Emory data
usd_emy_cph  = {  'neu_psy_tst_dte'     'Date of Neuropsych_preop' ; ...
                  'neu_psy_tst_dte_pst' 'Post-Operative NP test date' ; ...
                  'scn_dte'             'Neuroimaging Scan Date' ; ...
                  'srg_dte'             'Surgery Date' ; ...
                  'wtr_int_nor_scr'     '' ; ...
                  'wtr_int_nor_scr_pst' '' ; ...
                  'sbj_sze_frq'         'Seizure frequency' ; ...
                  'sbj_mts'             'MRI Results' ; ...
                  'eng_out'             'Engel Outcome' ; ...
                  'lng_lat'             'Language' ; ...
                  'lng_lat_mth'         'Language' ; ...
                  'sbj_eth'             'Ethnicity (Non-Hispanic/Hispanic/Unknown' ; ...
                  'sbj_rce'             'Race (White/Black/Asian/NH/Unknown' ; ...
                  'scn_typ'             'BITC/CMRR sequence' };
 
for iC = 1:numel(cln_nui_dta_col)
    org_col_ind = strcmpi(usd_emy_cph(:,1),cln_nui_dta_col{iC});
    emy_col_ind = strcmpi(emy_dta_col, usd_emy_cph{org_col_ind,2});
    num_val     = isnumeric(cln_nui_dta{1,iC});
    for iS = 1:size(emy_dta_sbj,1)
        sbj_ind = strcmpi(cln_nui_dta_sbj, emy_dta_sbj{iS});
                
        if isempty(find(emy_col_ind)) || isempty(emy_dta{iS,emy_col_ind})
            if num_val
                cln_nui_dta{sbj_ind,iC} = NaN;
            else
                cln_nui_dta{sbj_ind,iC} = '';
            end
        else
            cln_nui_dta{sbj_ind,iC} = emy_dta{iS,emy_col_ind};
        end               
        
    end
end

%% Fix Category Labels
% MTS
col_tgh{find(strcmpi(cln_nui_dta_col,'sbj_mts'))} = { 'L'        'yes' ; ...
                                                  'Left MTS' 'yes' ; ...
                                                  'left MTS' 'yes' ; ...
                                                  'left MTS and glioma involving optic chaism' 'yes' ; ...
                                                  'MTS' 'yes' ; ...
                                                  'R'   'yes' ; ...
                                                  'right MTS'    'yes' ; ...
                                                  'Right MTS'    'yes' ; ...
                                                  'right MTS and diffuse volume loss'    'yes' ; ...
                                                  'right MTS and small right TL'         'yes' ; ...
                                                  'Right MTS & Broad region of Ill-defined gray-white matter junction involving the right temporal lobe as described representative of secondary acquired temporal sclerosis.'    'yes' ;  ...
                                                  'normal MRI' 'no' ; ...
                                                  'other' 'no' ; ...
                                                  'N/A' 'no' ; ...
                                                  '' 'no'};
                                                 
% Language Laterality
col_tgh{find(strcmpi(cln_nui_dta_col,'lng_lat'))} = { 'bilateral (fMRI)'               'B' ; ...
                                                  'bilateral (left>right by Wada)' 'B' ; ...
                                                  'bilateral (left>right) (Wada)'  'B' ; ...
                                                  'bilateral (Wada)'               'B' ; ...
                                                  'left (fmri/wada)' 'L' ; ...
                                                  'Left (Wada & fMRI)' 'L' ; ...
                                                  'left (fMRI and Wada)' 'L' ; ...
                                                  'Left (fMRI and Wada)' 'L' ; ...
                                                  'Left (fMRI and Wada); failed ICA Wada ; but passed PCA' 'L' ; ...
                                                  'left (fMRI)' 'L' ; ...
                                                  'Left (fMRI)' 'L' ; ...
                                                  'right (Wada)' 'R' ; ...
                                                  'unknown no fMRI/Wada' '' ; ...
                                                  'right (fmri)' 'R' };
                                              
% Language Laterality Method
col_tgh{find(strcmpi(cln_nui_dta_col,'lng_lat_mth'))} = { 'bilateral (fMRI)'               'fMRI' ; ...
                                                      'bilateral (left>right by Wada)' 'Wada' ; ...
                                                      'bilateral (left>right) (Wada)'  'Wada' ; ...
                                                      'bilateral (Wada)'               'Wada' ; ...
                                                      'left (fMRI and Wada)' 'Wada' ; ...
                                                      'Left (fMRI and Wada)' 'Wada' ; ...
                                                      'Left (fMRI and Wada); failed ICA Wada ; but passed PCA' 'Wada' ; ...
                                                      'left (fmri/wada)' 'Wada' ; ...
                                                      'Left (Wada & fMRI)' 'Wada' ; ...
                                                      'left (fMRI)' 'fMRI' ; ...
                                                      'Left (fMRI)' 'fMRI' ; ...
                                                      'right (Wada)' 'fMRI' ; ...
                                                      'unknown no fMRI/Wada' '' ; ...
                                                      'right (fmri)' 'fMRI' };

% Engel outcome
col_tgh{find(strcmpi(cln_nui_dta_col,'eng_out'))} = { '1A'               'I' ; ...
                                                      '1B'               'I' ; ...
                                                      '1D'               'I' ; ...
                                                      '2A'               'II' ; ...
                                                      '2B'               'II' ; ...
                                                      '2D'               'II' ; ...
                                                      '3A'               'III' ; ...
                                                      '4B'               'IV' ; ...
                                                      '4C'               'IV' ; ...
                                                      }; 


% Run
for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            cln_nui_dta(strcmpi(cln_nui_dta(:,iC),col_tgh{iC}{iT,1}),iC) = {col_tgh{iC}{iT,2}};
        end
    end
end

clear col_tgh

% Simplified Engle Outcome
cln_nui_dta_col = [ cln_nui_dta_col 'eng_out_smp' ];
cln_nui_dta = [ cln_nui_dta cln_nui_dta(:,strcmpi(cln_nui_dta_col,'eng_out'))];    

col_tgh{find(strcmpi(cln_nui_dta_col,'eng_out_smp'))} = { 'I'               'I'   ; ...
                                                          'II'              'II+' ; ...
                                                          'III'             'II+' ; ...
                                                          'IV'              'II+' }; 

% Run
for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            cln_nui_dta(strcmpi(cln_nui_dta(:,iC),col_tgh{iC}{iT,1}),iC) = {col_tgh{iC}{iT,2}};
        end
    end
end


%% Fix Emory Date
dte_fld = { 'neu_psy_tst_dte' 'neu_psy_tst_dte_pst' 'scn_dte' 'srg_dte' };

for iF = 1:numel(dte_fld)
    dte_col = strcmpi(cln_nui_dta_col,dte_fld{iF});
    for iS = 1:size(cln_nui_dta_sbj,1)
        str_hld = regexpi(cln_nui_dta{iS,dte_col},'/','split');
        if numel(str_hld)==3
            cln_nui_dta{iS,dte_col} = [ str_hld{3} '-' str_hld{1} '-' str_hld{2} ];
        end
    end
end

%% Save out for Dan
% sbj_ind = [ grp.surgery.pst_cog.ltle_atl ; grp.surgery.pst_cog.ltle_slah ; grp.surgery.pst_cog.rtle_atl ; grp.surgery.pst_cog.rtle_slah];
% dan_dta = [ {'sbj_nme'}  cln_nui_dta_col ; cln_nui_dta_sbj(sbj_ind,:) cln_nui_dta(sbj_ind,:) ];
% cell2csv(['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Analysis/output/Data/Emory/additional_variables.csv'],dan_dta)

%% Calculate
% Distance-to-duration calculations
dur_col_one = { 'srg_dte'             'neu_psy_tst_dte' 'scn_dte'    'neu_psy_tst_dte' };
dur_col_two = { 'neu_psy_tst_dte_pst' 'srg_dte'         'srg_dte'    'scn_dte' };
dur_dta_col = { 'srg_to_pst'          'pre_to_srg'      'img_to_srg' 'pre_to_img' };
dur_dta     = cell( numel(cln_nui_dta_sbj), numel(dur_dta_col) );

for iC = 1:numel(dur_dta_col)
    
    fst_col = strcmpi(cln_nui_dta_col,dur_col_one{iC});
    scd_col = strcmpi(cln_nui_dta_col,dur_col_two{iC});
    
    for iS = 1:numel(cln_nui_dta_sbj)
               
        % Duration: Surgery-to-PostOp %%%%%%%%%%%%%%
        if ( ~isempty(cln_nui_dta{iS,fst_col}) && ~strcmpi(cln_nui_dta{iS,fst_col},'') ) && ...
           ( ~isempty(cln_nui_dta{iS,scd_col}) && ~strcmpi(cln_nui_dta{iS,scd_col},'') )
            
            fst_dte = cellfun(@(x) str2num(x),regexp(cln_nui_dta{iS,fst_col},'-','split'));
            if numel(fst_dte)==1
                fst_dte = cellfun(@(x) str2num(x),regexp(cln_nui_dta{iS,fst_col},'/','split'));
                fst_dte = (fst_dte(3)*365) + (fst_dte(1)*30) + fst_dte(2);
            elseif numel(fst_dte)==3
                fst_dte = (fst_dte(1)*365) + (fst_dte(2)*30) + fst_dte(3);
            end
            
            
            scd_dte = cellfun(@(x) str2num(x),regexp(cln_nui_dta{iS,scd_col},'-','split'));
            if numel(scd_dte)==1;
                scd_dte = cellfun(@(x) str2num(x),regexp(cln_nui_dta{iS,scd_col},'/','split'));
                scd_dte = (scd_dte(3)*365) + (scd_dte(1)*30) + scd_dte(2);
            elseif numel(scd_dte)==3
                scd_dte = (scd_dte(1)*365) + (scd_dte(2)*30) + scd_dte(3);
            end
            
            dur_dta{iS,iC} = roundsd((fst_dte - scd_dte) / 365,3);
        else
            dur_dta{iS,iC} = nan;
        end
        
    end
end

% Output usuable Numbers
clc

all_ind = [ grp.surgery.pst_cog_dti_slm.ltle_atl  ; ...
            grp.surgery.pst_cog_dti_slm.ltle_slah ; ...
            grp.surgery.pst_cog_dti_slm.rtle_atl  ; ...
            grp.surgery.pst_cog_dti_slm.rtle_slah ];

for iD = 1:numel(dur_dta_col)
    men_dur = roundsd(nanmean( cell2mat(dur_dta(all_ind,iD)))*12,3);
    std_dur = roundsd(nanstd(  cell2mat(dur_dta(all_ind,iD)))*12,3);
    fprintf('%s : median=%f; std: %f\n',dur_dta_col{iD},men_dur,std_dur)
end

%% Subject #'s
% Per Site & Sequence
tbl_hld = tabulate(cln_nui_dta(all_ind,strcmpi(cln_nui_dta_col,'scn_typ')));


usd = numel(string_find(cln_nui_dta_sbj(all_ind),'epd'))-numel(string_find(cln_nui_dta_sbj(all_ind),'epd_ucsf'));
usf = numel(string_find(cln_nui_dta_sbj(all_ind),'epd_ucsf'));
bit = sum(strcmpi(cln_nui_dta(all_ind,strcmpi(cln_nui_dta_col,'scn_typ')),'bitc'));
cmr = sum(strcmpi(cln_nui_dta(all_ind,strcmpi(cln_nui_dta_col,'scn_typ')),'cmrr'));

fprintf('\n\nUCSD: %i\n',usd)
fprintf('UCSF: %i\n',usf)
fprintf('Emory: %i\n',bit+cmr)
fprintf('BITC: %i\n',bit)
fprintf('CMRR: %i\n',cmr)

%% Categorical Variables
% Engel Outcomes / # Decliner's
tbl_typ = { 'count' 'lng_lat'      'L/R/B' 1; ...
            'count' 'eng_out_smp'  'I/II+' 1  };

% Stats: Fishers
fld_nme = fieldnames(grp.surgery.pst_cog_dti_slm);
fsh_var = [find(strcmpi(cln_nui_dta_col,'lng_lat')) find(strcmpi(cln_nui_dta_col,'eng_out_smp'))];

[~, lvl_ind] = ismember(tbl_typ(:,2),cln_nui_dta_col(fsh_var));
        clear cln_lvl; for iLI = 1:numel(lvl_ind); if lvl_ind(iLI); cln_lvl(lvl_ind(iLI)) = tbl_typ(iLI,3); end; end

fcfg = [];
fcfg.grp     = grp.surgery.pst_cog_dti_slm;
fcfg.grp_inc = {fld_nme};
fcfg.grp_nme = {fld_nme};
fcfg.dta = cln_nui_dta(:,fsh_var);
fcfg.sbj = cln_nui_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

grp_dta{1}(cellfun(@isempty,grp_dta{1})) = {'N/A'};

fcfg = [];
fcfg.sbj = grp_sbj{1};
fcfg.grp_nme = grp_nme{iG};
fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cln_nui_dta_col(fsh_var);
fcfg.lvl     = cln_lvl;
fcfg.dta_two = repmat(grp_typ{1},1,numel(fsh_var));
fcfg.lbl_two = strcat( 'group_', cln_nui_dta_col(fsh_var));
fcfg.out_dir = [ prj_dir '/' prj_nme  '/' 'FinalAnalysis_restricted' '/'  'stats' '/' 'Clinical_Fishers' '_' 'NuisanceTest' '/'  ];
ejk_fisher_test( fcfg );

% Stats: Fishers L-TLE only
fld_nme = fieldnames(grp.surgery.pst_cog_dti_slm);
fsh_var = [find(strcmpi(cln_nui_dta_col,'lng_lat')) find(strcmpi(cln_nui_dta_col,'eng_out_smp'))];

[~, lvl_ind] = ismember(tbl_typ(:,2),cln_nui_dta_col(fsh_var));
        clear cln_lvl; for iLI = 1:numel(lvl_ind); if lvl_ind(iLI); cln_lvl(lvl_ind(iLI)) = tbl_typ(iLI,3); end; end

fcfg = [];
fcfg.grp     = grp.surgery.pst_cog_dti_slm;
fcfg.grp_inc = {fld_nme(1:2)};
fcfg.grp_nme = {fld_nme(1:2)};
fcfg.dta = cln_nui_dta(:,fsh_var);
fcfg.sbj = cln_nui_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

grp_dta{1}(cellfun(@isempty,grp_dta{1})) = {'N/A'};

fcfg = [];
fcfg.sbj = grp_sbj{1};
fcfg.grp_nme = grp_nme{iG};
fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = cln_nui_dta_col(fsh_var);
fcfg.lvl     = cln_lvl;
fcfg.dta_two = repmat(grp_typ{1},1,numel(fsh_var));
fcfg.lbl_two = strcat( 'group_', cln_nui_dta_col(fsh_var));
fcfg.out_dir = [ prj_dir '/' prj_nme  '/' 'FinalAnalysis_restricted' '/'  'stats' '/' 'Clinical_Fishers' '_' 'NuisanceTest' '_' 'Left '/'  ];
ejk_fisher_test( fcfg );

% Table



%% Save out
% out_dta = [ {'sbj_nme'}  cln_nui_dta_col ; cln_nui_dta_sbj cln_nui_dta ];
% cell2csv(['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Analysis/output/Data/additional_variables_complete.csv'],out_dta)
% 
%% Calculate partial correlations
grp_nme     = { 'surgery' };
grp_nme_sub = { {'pst_cog_dti'} };

% [ cln_nui_dta, cln_nui_dta_sbj, cln_nui_dta_col ]

% Load Final Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv'];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'additional_variables.csv'];
fcfg.dta_col = 2;
[ cln_nui_dta, cln_nui_dta_sbj, cln_nui_dta_col] = ejk_dta_frm( fcfg );

% Combine
fnl_dta     = [ fnl_dta     cln_nui_dta];
fnl_dta_col = [ fnl_dta_col cln_nui_dta_col ];

% Partial Correlations
cog_col = { 'lm2_chg' 'lm2_rci' };
neu_col = { 'Hippocampus' 'ILF' 'Unc' 'fusiform' 'lateralorbitofrontal' };
cov_col = { 'lm2_pre' 'NumAEDs' 'AgeOfSeizureOnset' 'Educ' 'AgeAtSurgery' };

cog_col = find(ismember(fnl_dta_col,cog_col));
neu_col = find(ismember(fnl_dta_col,neu_col));
cov_col = find(ismember(fnl_dta_col,cov_col));

% Run Correlations
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        for iFN = 1:numel(fld_nme)
            grp_ind = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iFN});
            
            fcfg = [];
            
            fcfg.sbj_nme = fnl_dta_sbj( grp_ind, 1);
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_one = cell2mat(fnl_dta( grp_ind, cog_col));
            fcfg.lbl_one = fnl_dta_col(cog_col);
            
            fcfg.dta_two = cell2mat(fnl_dta( grp_ind, neu_col));
            fcfg.lbl_two = fnl_dta_col(neu_col);
            
            fcfg.dta_cov = cell2mat(fnl_dta( grp_ind, cov_col));
            fcfg.lbl_cov = fnl_dta_col(cov_col);
            
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'partial_correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '/' fld_nme{iFN}  '/' ];
            
            ejk_partial_correlation( fcfg );
            
        end
        
    end
end








