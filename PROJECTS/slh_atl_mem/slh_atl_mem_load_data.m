%% Get participants
% Load UCSD/UCSF Redcap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.red_fle = red_cap_fle;
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

% Check who has post operative scores
sbj_nme_hld = cell(0);
for iT = 1:numel(cog_tst_nme)
    sbj_nme_hld = [ sbj_nme_hld ; sbj_cog.sbj_nme(~isnan(sbj_cog.(cog_tst_nme{iT}))) ];
end
sbj_nme = unique(sbj_nme_hld);

clear sbj_cog

% Re-Load Redcap
fcfg = [];
fcfg.sbj_nme = sbj_nme;
fcfg.red_fle = red_cap_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

% Load Emory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Emory' '/' 'Emory_SLAH_Memory_11.22.21_cleanNONAMES.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );

%% Add on Emory
nme_add = { 'sbj_cog' 'sbj_dem' 'sbj_sze' 'sbj_srg' };
hld_add.sbj_cog = sbj_cog;
hld_add.sbj_dem = sbj_dem;
hld_add.sbj_sze = sbj_sze;
hld_add.sbj_srg = sbj_srg;

usd_emy_cph{1}  = { 'neu_psy_tst_dte'     'Date of Neuropsych_preop' ; ...
    'neu_psy_tst_dte_pst' '' ; ...
    'wms_ver_log_mem'     'WMS version' ; ...
    'wms_ver_vpa'         'WMS version' ; ...
    'wms_ver_log_mem_pst' 'WMS version post' ; ...
    'wms_ver_vpa_pst'     'WMS version post' ; ...
    'log_mem_raw_scr_one' 'LM I raw_pre' ; ...
    'log_mem_nor_scr_one' 'LM I norm_ss_pre' ; ...
    'log_mem_raw_scr_two' 'LM II raw_pre' ; ...
    'log_mem_nor_scr_two' 'LM II norm_ss_pre' ; ...
    'vp1_raw_scr' 'VPA I raw_pre' ; ...
    'vp1_nor_scr' 'VPA I norm_ss_pre' ; ...
    'vp2_raw_scr' 'VPA II raw_pre' ; ...
    'vp2_nor_scr' 'VPA II norm_ss_pre'; ...
    'cvl_tot_raw_scr' '' ; ...
    'cvl_tot_nor_scr' '' ; ...
    'cvl_lfr_raw_scr' '' ; ...
    'cvl_lfr_nor_scr' '' ; ...
    ...
    'log_mem_raw_scr_one_pst' 'LM I raw_6m.post' ; ...
    'log_mem_nor_scr_one_pst' 'LM I norm_ss_6m.post' ; ...
    'log_mem_raw_scr_two_pst' 'LM II raw_6m.post' ; ...
    'log_mem_nor_scr_two_pst' 'LM II norm_ss_6m.post' ; ...
    'vp1_raw_scr_pst' 'VPA I raw_1y.post' ; ...
    'vp1_nor_scr_pst' 'VPA I norm_ss__1y.post' ; ...
    'vp2_raw_scr_pst' 'VPA II raw_1y.post' ; ...
    'vp2_nor_scr_pst' 'VPA II norm_ss__1y.post' ; ...
    'cvl_tot_raw_scr_pst' '' ; ...
    'cvl_tot_nor_scr_pst' '' ; ...
    'cvl_lfr_raw_scr_pst' '' ; ...
    'cvl_lfr_nor_scr_pst' '' };

usd_emy_cph{2}  = { 'sbj_edu'     'Education' ; ...
    'sbj_sex'     'Sex (F/M)' ; ...
    'sbj_hnd'     'Handedness (R/L)' };

usd_emy_cph{3}  = { 'sbj_sde_ons' 'Side Of Seizure Onset' ; ...
    'sbj_age_ons' 'Age of Seizure Onset' ; ...
    'sbj_aed_num' 'Number of ASMs' ; ...
    'sbj_sze_frq' 'Seizure frequency' ; ...
    'sbj_mts'     '' };

usd_emy_cph{4}  = { 'srg_age' 'Age' ; ...
    'srg_typ' 'Surgery' ; ...
    'srg_sde' 'Side_Surgery' ; ...
    'eng_out' '' };

for iFN = 1:numel(nme_add)
    
    sbj_num = numel(hld_add.(nme_add{iFN}).sbj_nme);
    
    for iS = 1:size(emy_dta_sbj,1)
        hld_add.(nme_add{iFN}).sbj_nme{sbj_num+iS,1} = emy_dta_sbj{iS,1};
        
        for iAD = 1:size(usd_emy_cph{iFN},1)
            
            if ~strcmpi(class(hld_add.(nme_add{iFN}).(usd_emy_cph{iFN}{iAD,1})),'double') && any(strcmpi(unique(cellfun(@class,hld_add.(nme_add{iFN}).(usd_emy_cph{iFN}{iAD,1}),'uni',0)),'char'))
                if isempty(usd_emy_cph{iFN}{iAD,2})
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} ){sbj_num+iS} = '';
                elseif isempty(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})}) || strcmpi(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})},'n/a')
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} ){sbj_num+iS} = '';
                else
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} ){sbj_num+iS} = emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})};
                end
            elseif strcmpi(class(hld_add.(nme_add{iFN}).(usd_emy_cph{iFN}{iAD,1})),'double')
                if isempty(usd_emy_cph{iFN}{iAD,2})
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
                elseif isempty(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})})
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
                elseif strcmpi(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})},'n/a')
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
                else
                    hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})};
                end
            else
                error('Check class agreement')
            end
        end
    end
    
end

sbj_cog = hld_add.sbj_cog;
sbj_dem = hld_add.sbj_dem;
sbj_sze = hld_add.sbj_sze;
sbj_srg = hld_add.sbj_srg;

%% Save Cognitive Spreadsheet
% Clean up WMS versions
for iS = 1:numel(sbj_cog.sbj_nme)
    % Pre-op
    if isnumeric(sbj_cog.wms_ver_log_mem{iS})
        if sbj_cog.wms_ver_log_mem{iS}==3
            sbj_cog.wms_ver_log_mem{iS} = 'III';
            sbj_cog.wms_ver_vpa{iS}     = 'III';
        elseif sbj_cog.wms_ver_log_mem{iS}==4
            sbj_cog.wms_ver_log_mem{iS} = 'IV';
            sbj_cog.wms_ver_vpa{iS}     = 'IV';
        elseif isnan(sbj_cog.wms_ver_log_mem{iS})
            if ~isnan(sbj_cog.log_mem_nor_scr_two(iS))
                sbj_cog.wms_ver_log_mem{iS} = 'III';
            else
                sbj_cog.wms_ver_log_mem{iS} = '';
            end
            if ~isnan(sbj_cog.vp2_nor_scr(iS))
                sbj_cog.wms_ver_vpa{iS}     = 'III';
            else
                sbj_cog.wms_ver_vpa{iS}     = '';
            end
        elseif ~isnan(sbj_cog.wms_ver_log_mem{iS})
            error('Other WMS version numeric')
        else
            sbj_cog.wms_ver_log_mem{iS} = '';
            sbj_cog.wms_ver_vpa{iS}     = '';
        end
    elseif ~isempty(strfind(sbj_cog.wms_ver_log_mem{iS},';'))
        pre_wms_hld = sbj_cog.wms_ver_log_mem{iS};
        if contains(pre_wms_hld,'WMS-IV LM') || contains(pre_wms_hld,'WSM-IV LM') || contains(pre_wms_hld,'WSM -IV LM')
            sbj_cog.wms_ver_log_mem{iS} = 'IV';
        elseif contains(pre_wms_hld,'WMS-III LM')
            sbj_cog.wms_ver_log_mem{iS} = 'III';
        elseif isnan(sbj_cog.log_mem_nor_scr_two(iS))
            sbj_cog.wms_ver_log_mem{iS} = '';
        else
            error('Problem with LM WMS')
        end
        pre_wms_hld = sbj_cog.wms_ver_vpa{iS};
        if contains(pre_wms_hld,'WMS-IV VPA') || contains(pre_wms_hld,'WSM-IV VPA') || contains(pre_wms_hld,'WSM -IV VPA') || contains(pre_wms_hld,'WMS-IV LM and VPA') || contains(pre_wms_hld,'WMS-IV LM; VR and VPA')  || contains(pre_wms_hld,'WMS-IV LM; WMS- VR and VPA')
            sbj_cog.wms_ver_vpa{iS} = 'IV';
        elseif contains(pre_wms_hld,'WMS-III VPA') || contains(pre_wms_hld,'WMS-III LM & VPA')
            sbj_cog.wms_ver_vpa{iS} = 'III';
        elseif contains(pre_wms_hld,'WMS-R VPA')
            sbj_cog.wms_ver_vpa{iS} = 'II';
        elseif isnan(sbj_cog.vp2_nor_scr(iS))
            sbj_cog.wms_ver_vpa{iS} = '';
        else
            error('Problem with VPA WMS')
        end
    end
    % Post-op
    if isnumeric(sbj_cog.wms_ver_log_mem_pst{iS})
        if sbj_cog.wms_ver_log_mem_pst{iS}==3
            sbj_cog.wms_ver_log_mem_pst{iS} = 'III';
            sbj_cog.wms_ver_vpa_pst{iS}     = 'III';
        elseif sbj_cog.wms_ver_log_mem_pst{iS}==4
            sbj_cog.wms_ver_log_mem_pst{iS} = 'IV';
            sbj_cog.wms_ver_vpa_pst{iS}     = 'IV';
        elseif isnan(sbj_cog.wms_ver_log_mem{iS})
            if ~isnan(sbj_cog.log_mem_nor_scr_two_pst(iS))
                sbj_cog.wms_ver_log_mem_pst{iS} = 'III';
            else
                sbj_cog.wms_ver_log_mem_pst{iS} = '';
            end
            if ~isnan(sbj_cog.vp2_nor_scr_pst(iS))
                sbj_cog.wms_ver_vpa_pst{iS}     = 'III';
            else
               sbj_cog.wms_ver_vpa_pst{iS}     = ''; 
            end
        elseif ~isnan(sbj_cog.wms_ver_log_mem_pst{iS})
            error('Other WMS version numeric')
        else
            sbj_cog.wms_ver_log_mem_pst{iS} = '';
            sbj_cog.wms_ver_vpa_pst{iS}     = '';
        end
    elseif ~isempty(strfind(sbj_cog.wms_ver_log_mem_pst{iS},';'))
        pre_wms_hld = sbj_cog.wms_ver_log_mem_pst{iS};
        if contains(pre_wms_hld,'WMS-IV LM') || contains(pre_wms_hld,'WSM-IV LM') || contains(pre_wms_hld,'WSM -IV LM')
            sbj_cog.wms_ver_log_mem_pst{iS} = 'IV';
        elseif contains(pre_wms_hld,'WMS-III LM')
            sbj_cog.wms_ver_log_mem_pst{iS} = 'III';
        elseif isnan(sbj_cog.log_mem_nor_scr_two_pst(iS))
            sbj_cog.wms_ver_log_mem_pst{iS} = '';
        else
            error('Problem with LM WMS post')
        end
        pre_wms_hld = sbj_cog.wms_ver_vpa_pst{iS};
        if contains(pre_wms_hld,'WMS-IV VPA') || contains(pre_wms_hld,'WSM-IV VPA') || contains(pre_wms_hld,'WSM -IV VPA') || contains(pre_wms_hld,'WMS-IV LM and VPA') || contains(pre_wms_hld,'WMS-IV LM; VR and VPA') || contains(pre_wms_hld,'WSM-IV LM & VPA') || contains(pre_wms_hld,'WMS-4 VPA')  || contains(pre_wms_hld,'WMS-IV LM; WMS- VR and VPA')
            sbj_cog.wms_ver_vpa_pst{iS} = 'IV';
        elseif contains(pre_wms_hld,'WMS-III VPA') || contains(pre_wms_hld,'WMS-III LM & VPA')
            sbj_cog.wms_ver_vpa_pst{iS} = 'III';
        elseif contains(pre_wms_hld,'WMS-R VPA') || contains(pre_wms_hld,'WMS-R - VPA')
            sbj_cog.wms_ver_vpa_pst{iS} = 'II';
        else
            error('Problem with VPA WMS post')
        end
    end
end

% Convert CVLT Z to SS
fcfg = [];
fcfg.typ = 'z_to_ss';
sbj_cog.cvl_lfr_nor_scr     = ejk_convert_neuropsych(fcfg,sbj_cog.cvl_lfr_nor_scr);
sbj_cog.cvl_lfr_nor_scr_pst = ejk_convert_neuropsych(fcfg,sbj_cog.cvl_lfr_nor_scr_pst);

% Calculate Post-op Cognitive
fcfg = [];
pst_cog_dta = ejk_post_cognitive(fcfg,sbj_cog);

% Initial put together
cog_out_col = { 'lm2_pre' 'lm2_pst' 'lm2_chg' 'lm2_rci' 'lm2_pct' ...
                'vp2_pre' 'vp2_pst' 'vp2_chg' 'vp2_rci' 'vp2_pct' ...
                'cv2_pre' 'cv2_pst' 'cv2_chg' 'cv2_rci' 'cv2_pct' };
cog_out_dta = [ sbj_cog.log_mem_nor_scr_two sbj_cog.log_mem_nor_scr_two_pst pst_cog_dta.raw.log_mem_nor_scr_two_pst pst_cog_dta.rci.log_mem_nor_scr_two_pst pst_cog_dta.pct.log_mem_nor_scr_two_pst ...
                sbj_cog.vp2_nor_scr         sbj_cog.vp2_nor_scr_pst         pst_cog_dta.raw.vp2_nor_scr_pst         pst_cog_dta.rci.vp2_nor_scr_pst         pst_cog_dta.pct.vp2_nor_scr_pst ...
                sbj_cog.cvl_lfr_nor_scr     sbj_cog.cvl_lfr_nor_scr_pst     pst_cog_dta.raw.cvl_lfr_nor_scr_pst     pst_cog_dta.rci.cvl_lfr_raw_scr_pst     pst_cog_dta.pct.cvl_lfr_nor_scr_pst ];

% Save out
cell2csv( [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'] , [ [{'sbj_nme'} cog_out_col] ; [pst_cog_dta.sbj_nme num2cell(cog_out_dta)] ])

%% Save out Clinical
% Initial put together
cln_dta_nme = { 'SubjID'      'VisitID'           'SideOfSeizureFocus' 'FieldStrength' 'AgeAtSurgery' ...
    'Educ'        'AgeOfSeizureOnset' 'NumAEDs'            'SeizureFreq'   'SurgicalSide' ...
    'SurgeryType' 'Sex'               'Handedness'         'MTS'           'EngelOutcome' ...
    'wms_lm2' 'wms_vp2' };
cln_dta_out = [ sbj_dem.sbj_nme           repmat({''},numel(sbj_dem.sbj_nme),1)          sbj_sze.sbj_sde_ons           repmat({NaN},numel(sbj_dem.sbj_nme),1)          num2cell(sbj_srg.srg_age) ...
    num2cell(sbj_dem.sbj_edu) num2cell(sbj_sze.sbj_age_ons) num2cell(sbj_sze.sbj_aed_num) num2cell(sbj_sze.sbj_sze_frq) sbj_srg.srg_sde ...
    sbj_srg.srg_typ           sbj_dem.sbj_sex               sbj_dem.sbj_hnd               sbj_sze.sbj_mts               sbj_srg.eng_out ...
    sbj_cog.wms_ver_log_mem  sbj_cog.wms_ver_vpa ];

% Code location
cln_dta_nme{end+1} = 'Location';
cln_dta_out(string_find(cln_dta_out(:,1),{'epd\d'}),end+1)      = {'UCSD'};
cln_dta_out(string_find(cln_dta_out(:,1),{'epd_ucsf\d'}),end) = {'UCSF'};
cln_dta_out(string_find(cln_dta_out(:,1),{'fc\d'}),end)       = {'HC'};
cln_dta_out(string_find(cln_dta_out(:,1),{'EPK\d'}),end)      = {'Emory'};

% Tighten up clinical coding
col_tgh{find(strcmpi(cln_dta_nme,'SideOfSeizureFocus'))} = { 'L' 'left' ;
    'R' 'right' };

col_tgh{find(strcmpi(cln_dta_nme,'SurgicalSide'))} = { 'L'   'left'  ;
    'R'   'right' ;
    ''    'n/a'};

col_tgh{find(strcmpi(cln_dta_nme,'SurgeryType'))} = { 'SLAH' 'left SLAH'   ; ...
    'SLAH' 'Right SLAH'  ; ...
    'ATL'  'Right ATL'    ; ...
    'SLAH' 'neuroablation (hippocampus & amygdala)'};

col_tgh{find(strcmpi(cln_dta_nme,'Sex'))} = { 'F'   'female' ;
    'M'   'male'   };

col_tgh{find(strcmpi(cln_dta_nme,'Handedness'))} = { 'L' 'left' ;
    'R' 'right' };

for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            cln_dta_out(strcmpi(cln_dta_out(:,iC),col_tgh{iC}{iT,2}),iC) = {col_tgh{iC}{iT,1}};
        end
    end
end

% Update WMS version

% Save out
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical' '.csv'], [ cln_dta_nme ; cln_dta_out ]);

%% Neurobio
% Reconfile update


%% Misc
% Emory - What's left?
emy_col_nme = [ usd_emy_cph{1}(:,2) ; usd_emy_cph{2}(:,2) ; usd_emy_cph{3}(:,2) ; usd_emy_cph{4}(:,2) ];
emy_col_nme(cellfun(@isempty,emy_col_nme)) = [];
emy_col_nme_rmn = setxor(emy_dta_col,emy_col_nme);
cell2csv( [prj_dir '/' prj_nme '/' 'Data' '/' 'EmoryUnusedVariables' '.csv'], emy_col_nme_rmn);

% Complicated Surgeries What's left?


