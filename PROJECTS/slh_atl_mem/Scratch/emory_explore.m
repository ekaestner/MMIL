out_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/EmoryExplore';

%% Load Clinical & Cognitive
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Emory' '/' 'Emory_SLAH_Memory_11.22.21_cleanNONAMES.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = '/home/ekaestne/PROJECTS/DATA/csv/redcap/Redcap_2022_01_10.csv';
fcfg.dta_col = 2;
[ red_cap_dta, red_cap_dta_sbj, red_cap_dta_col] = ejk_dta_frm( fcfg );

%% Concatenate Data
usd_emy_cln_cph = { 'VisitID'            'FreesurferRecon' ; ...
                    'SideOfSeizureFocus' 'Side Of Seizure Onset' ; ...
                    'FieldStrength'      '' ; ...
                    'AgeAtSurgery'       'Age' ; ...
                    'Educ'               'Education' ; ...
                    'AgeOfSeizureOnset'  'Age of Seizure Onset' ; ...
                    'NumAEDs'            'Number of ASMs' ; ...
                    'SeizureFreq'        'Seizure frequency' ; ...
                    'SurgicalSide'       'Side_Surgery' ; ...
                    'SurgeryType'        'Surgery' ; ...
                    'Sex'                'Sex (F/M)' ; ...
                    'Handedness'         'Handedness (R/L)' ; ...
                    'MTS'                '' ; ...
                    'EngelOutcome'       '' };

usd_emy_cog_cph = { 'log_mem_nor_scr_two'         'LM II norm_ss_pre' ; ...
                    'vp2_nor_scr'                 'VPA II norm_ss_pre'; ...
                    'log_mem_nor_scr_two_zsc'     '' ; ...
                    'vp2_nor_scr_zsc'             '' ; ...
                    'log_mem_nor_scr_two_cat'     '' ; ...
                    'vp2_nor_scr_cat'             '' ; ...
                    'log_mem_nor_scr_two_pst'     'LM II norm_ss_6m.post' ; ...
                    'vp2_nor_scr_pst'             'VPA II norm_ss__1y.post' ; ...
                    'log_mem_nor_scr_two_pst_cat' '' ; ...
                    'vp2_nor_scr_pst_cat'         '' };           
               
% concatenate clinical
cmb_cln_dta_sbj = [ cln_dta_sbj ; emy_dta_sbj];  
cmb_cln_dta_col = cln_dta_col;
for iC = 1:size(usd_emy_cln_cph,1)
    non_emp_ind = find(cellfun(@isempty,cln_dta(:,iC))); if isempty(non_emp_ind); non_emp_ind=1; end
    if ~isempty(usd_emy_cln_cph{iC,2})
        cmb_cln_dta(:,iC) = [ cln_dta( :, strcmpi(cln_dta_col,usd_emy_cln_cph{iC,1}) ) ; ...
                              emy_dta( :, strcmpi(emy_dta_col,usd_emy_cln_cph{iC,2}) ) ];
        if strcmpi( class(cln_dta{non_emp_ind(1),iC}), 'char')
            pos_usd_cde{iC} = unique(cln_dta( :, strcmpi(cln_dta_col,usd_emy_cln_cph{iC,1}) ));
            pos_emy_cde{iC} = unique(emy_dta( :, strcmpi(emy_dta_col,usd_emy_cln_cph{iC,2}) ));
        end
    elseif isempty(usd_emy_cln_cph{iC,2}) && strcmpi( class(cln_dta{non_emp_ind(1),iC}), 'double')
        cmb_cln_dta(:,iC) = [ cln_dta( :, strcmpi(cln_dta_col,usd_emy_cln_cph{iC,1}) ) ; ...
                              repmat({NaN},numel(emy_dta_sbj),1) ];
    elseif isempty(usd_emy_cln_cph{iC,2}) && strcmpi( class(cln_dta{non_emp_ind(1),iC}), 'char')
        cmb_cln_dta(:,iC) = [ cln_dta( :, strcmpi(cln_dta_col,usd_emy_cln_cph{iC,1}) ) ; ...
                              repmat({''},numel(emy_dta_sbj),1) ];
    end
        
    non_emp_ind = find(cellfun(@isempty,cln_dta(:,iC))); if isempty(non_emp_ind); non_emp_ind=1; end
    if strcmpi( class(cln_dta{non_emp_ind(1),iC}), 'char')
        cmb_cln_dta(cellfun(@isempty,cmb_cln_dta(:,iC)),iC) = {''};       
    elseif strcmpi( class(cln_dta{non_emp_ind(1),iC}), 'double')
        cmb_cln_dta(cellfun(@isempty,cmb_cln_dta(:,iC)),iC) = {NaN};
    end
end

% concatenate cognitive
cmb_cog_dta_sbj = [ cog_dta_sbj ; emy_dta_sbj];  
cmb_cog_dta_col = cog_dta_col;
for iC = 1:size(usd_emy_cog_cph,1)
    if ~isempty(usd_emy_cog_cph{iC,2})
        cmb_cog_dta(:,iC) = [ cog_dta( :, strcmpi(cog_dta_col,usd_emy_cog_cph{iC,1}) ) ; ...
                              emy_dta( :, strcmpi(emy_dta_col,usd_emy_cog_cph{iC,2}) ) ];
    elseif isempty(usd_emy_cog_cph{iC,2}) && strcmpi( class(cog_dta{1,iC}), 'double')
        cmb_cog_dta(:,iC) = [ cog_dta( :, strcmpi(cog_dta_col,usd_emy_cog_cph{iC,1}) ) ; ...
                              repmat({NaN},numel(emy_dta_sbj),1) ];
    elseif isempty(usd_emy_cog_cph{iC,2}) && strcmpi( class(cog_dta{1,iC}), 'char')
        cmb_cog_dta(:,iC) = [ cog_dta( :, strcmpi(cog_dta_col,usd_emy_cog_cph{iC,1}) ) ; ...
                              repmat({''},numel(emy_dta_sbj),1) ];
    end 
    
    non_emp_ind = find(cellfun(@isempty,cog_dta(:,iC))); if isempty(non_emp_ind); non_emp_ind=1; end
    if strcmpi( class(cog_dta{non_emp_ind(1),iC}), 'char')
        cmb_cog_dta(cellfun(@isempty,cmb_cog_dta(:,iC)),iC) = {''};
    elseif strcmpi( class(cog_dta{non_emp_ind(1),iC}), 'double')
        cmb_cog_dta(cellfun(@isempty,cmb_cog_dta(:,iC)),iC) = {NaN};
    end
end

cmb_cog_dta(cellfun(@(x) strcmpi(x,'n/a'),cmb_cog_dta)) = {NaN};

% label location
sbj_loc = cell(numel(cmb_cog_dta_sbj),1);
sbj_loc(string_find(cmb_cog_dta_sbj,{'epd\d'}))      = {'UCSD'};
sbj_loc(string_find(cmb_cog_dta_sbj,{'epd_ucsf\d'})) = {'UCSF'};
sbj_loc(string_find(cmb_cog_dta_sbj,{'fc\d'}))       = {'HC'};
sbj_loc(string_find(cmb_cog_dta_sbj,{'EPK\d'}))      = {'Emory'};

cmb_cln_dta_col = [cmb_cln_dta_col {'Location'}];
cmb_cln_dta     = [cmb_cln_dta sbj_loc];

%% WMS version - obtain
wms_hld = cell(numel(cmb_cln_dta_sbj),5); % sbj_nme pre_lm pre_vpa pst_lm pst_vpa
for iS = 1:size(wms_hld,1)
   wms_hld{iS,1} =  cmb_cln_dta_sbj{iS};
   if any(strcmpi(cln_dta_sbj,wms_hld{iS,1}))
       sbj_ind = find(strcmpi(red_cap_dta_sbj,wms_hld{iS,1}));
       wms_ind(1) = find(strcmpi(red_cap_dta_col,'wms_version'));
       wms_ind(2) = find(strcmpi(red_cap_dta_col,'post_wms_version'));
       
       pre_wms_hld = red_cap_dta{sbj_ind,wms_ind(1)};
       pst_wms_hld = red_cap_dta{sbj_ind,wms_ind(2)};
       
       if (ischar(pre_wms_hld) && contains(pre_wms_hld,'IV')) || (~isempty(pre_wms_hld) && isnumeric(pre_wms_hld) && pre_wms_hld==4)
            wms_hld{iS,2} = 'IV';     
            wms_hld{iS,3} = 'IV';     
       elseif (ischar(pre_wms_hld) && contains(pre_wms_hld,'III')) || (~isempty(pre_wms_hld) && isnumeric(pre_wms_hld) && pre_wms_hld==3)
           wms_hld{iS,2} = 'III'; 
           wms_hld{iS,3} = 'III'; 
       elseif isempty(pre_wms_hld)
           wms_hld{iS,2} = '';
           wms_hld{iS,3} = ''; 
       else
           wms_hld{iS,2} = error(pre_wms_hld);
       end
             
       if (ischar(pst_wms_hld) && contains(pst_wms_hld,'IV')) || (~isempty(pst_wms_hld) && isnumeric(pst_wms_hld) && pst_wms_hld==4)
            wms_hld{iS,4} = 'IV';     
            wms_hld{iS,5} = 'IV';      
       elseif (ischar(pst_wms_hld) && contains(pst_wms_hld,'III')) || (~isempty(pst_wms_hld) && isnumeric(pst_wms_hld) && pst_wms_hld==3)
          wms_hld{iS,4} = 'III';     
            wms_hld{iS,5} = 'III';  
       elseif isempty(pst_wms_hld)
           wms_hld{iS,4} = '';
           wms_hld{iS,5} = ''; 
       else
            error(pst_wms_hld)
       end       
       
   elseif any(strcmpi(emy_dta_sbj,wms_hld{iS,1}))
       sbj_ind = find(strcmpi(emy_dta_sbj,wms_hld{iS,1}));
       wms_ind = find(strcmpi(emy_dta_col,'WMS version'));
       
       pre_wms_hld = emy_dta{sbj_ind,wms_ind(1)};
       pst_wms_hld = emy_dta{sbj_ind,wms_ind(2)};
       
       if contains(pre_wms_hld,'WMS-IV LM') || contains(pre_wms_hld,'WSM-IV LM') || contains(pre_wms_hld,'WSM -IV LM')
            wms_hld{iS,2} = 'IV';     
       elseif contains(pre_wms_hld,'WMS-III LM')
           wms_hld{iS,2} = 'III'; 
       else
           wms_hld{iS,2} = pre_wms_hld;
       end
       
       if contains(pre_wms_hld,'WMS-IV VPA') || contains(pre_wms_hld,'WSM-IV VPA') || contains(pre_wms_hld,'WSM -IV VPA') || contains(pre_wms_hld,'WMS-IV LM and VPA') || contains(pre_wms_hld,'WMS-IV LM; VR and VPA')
            wms_hld{iS,3} = 'IV';     
       elseif contains(pre_wms_hld,'WMS-III VPA') || contains(pre_wms_hld,'WMS-III LM & VPA')
           wms_hld{iS,3} = 'III';        
       elseif contains(pre_wms_hld,'WMS-R VPA')
           wms_hld{iS,3} = 'II'; 
       else
           wms_hld{iS,3} = pre_wms_hld;
       end
      
       if contains(pst_wms_hld,'WMS-IV LM') || contains(pst_wms_hld,'WSM-IV LM') || contains(pst_wms_hld,'WSM -IV LM')
            wms_hld{iS,4} = 'IV';     
       elseif contains(pst_wms_hld,'WMS-III LM')
           wms_hld{iS,4} = 'III'; 
       else
            error(pst_wms_hld)
           wms_hld{iS,4} = pst_wms_hld;
       end
       
       if contains(pst_wms_hld,'WMS-IV VPA') || contains(pst_wms_hld,'WSM-IV VPA') || contains(pst_wms_hld,'WSM -IV VPA') || contains(pst_wms_hld,'WMS-IV LM and VPA') || contains(pst_wms_hld,'WMS-IV LM & VPA') || contains(pst_wms_hld,'WSM-IV LM & VPA') || contains(pst_wms_hld,'WMS-IV LM; VR and VPA') || contains(pst_wms_hld,'WMS-4 VPA')
            wms_hld{iS,5} = 'IV';     
       elseif contains(pst_wms_hld,'WMS-III VPA') || contains(pst_wms_hld,'WMS-III LM & VPA')
           wms_hld{iS,5} = 'III';        
       elseif contains(pst_wms_hld,'WMS-R VPA') || contains(pst_wms_hld,'WMS-R - VPA')
           wms_hld{iS,5} = 'II'; 
       else
           wms_hld{iS,5} = pst_wms_hld;
       end
       
   end
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/output/Data/wms_version.csv',wms_hld);

%% Tighten up clinical coding
col_tgh{2} = { 'L' 'left' ;
               'R' 'right' };
           
col_tgh{9} = { 'L'   'left'  ;
               'R'   'right' ;
               ''    'n/a'};

col_tgh{10} = { 'SLAH' 'left SLAH'   ; ...
                'SLAH' 'Right SLAH'  ; ...
                'ATL'  'Right ATL'    ; ...
                'SLAH' 'neuroablation (hippocampus & amygdala)'};           

col_tgh{11} = { 'F'   'female' ;
                'M'   'male'   };

col_tgh{12} = { 'L' 'left' ;
                'R' 'right' };

for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            cmb_cln_dta(strcmpi(cmb_cln_dta(:,iC),col_tgh{iC}{iT,2}),iC) = {col_tgh{iC}{iT,1}};
        end
    end
end

%% Reconfile update
mri_dev = mmil_readtext(mri_dev_fle);

rcn_out = cell(numel(emy_dta_sbj),3);
for iS = 1:numel(emy_dta_sbj)
    tst_str = emy_dta{iS,strcmpi(emy_dta_col,'Date of Neuropsych_preop')};
    tst_num = cellfun(@(x) str2num(x),regexp(tst_str,'/','split'));
    tst_num = (tst_num(3)*365) + (tst_num(1)*30) + tst_num(2);
    
    pos_ind = string_find(mri_dev(:,1),emy_dta_sbj{iS});
 
    
    fprintf('\n')
    fprintf('%s - %s\n',emy_dta_sbj{iS},tst_str)
    for iPS = 1:numel(pos_ind)
        scn_str = regexp(mri_dev{pos_ind(iPS),1},'20[\d{6}](?![_])');
        scn_str = mri_dev{pos_ind(iPS),1}(scn_str:scn_str+7);
        scn_num = (str2num(scn_str(1:4))*365) + (str2num(scn_str(5:6))*30) + str2num(scn_str(7:8));
        fprintf('%s : %d days\n',mri_dev{pos_ind(iPS),1},scn_num-tst_num)
        day_hld(iPS) = scn_num-tst_num;
    end
    fprintf('\n')
    
    if ~isempty(pos_ind)
        [ ~, use_ind] = min(abs(day_hld));
        rcn_out(iS,:) = [ emy_dta_sbj(iS) mri_dev{pos_ind(use_ind),2} {day_hld(use_ind)} ];
    else
        rcn_out(iS,:) = [ emy_dta_sbj(iS) {''} {NaN} ];
    end
    clear day_hld
    
end

cell2csv( [ out_dir '/' 'emory_recon_files.csv' ], rcn_out );

%% Load Neuroimaging Data
inp_fle     = { mri_dev_fle dti_dev_fle };
inp_pre     = { mri_dev_pre dti_dev_pre };
inp_suf     = { mri_dev_suf dti_dev_suf };
inp_mse     = { mri_dev_mse dti_dev_mse };
inp_roi     = { mri_dev_roi dti_dev_roi };
inp_lat     = { mri_dev_lat dti_dev_lat };
inp_icv     = { mri_dev_icv dti_dev_icv };

for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [ out_dir '/' 'Data' '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [ out_dir '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '/' ];
            else
                roi_fle_nme = [ out_dir '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [ out_dir '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM}  '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '/' ];
            end
            
            if inp_roi{iF}{iM}(iR)==0
                roi_nme_hld = [];
                roi_hms_hld = [];
            else
                roi_nme_hld = { [ roi_loc '/' 'lh.' roi_nme{inp_roi{iF}{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{inp_roi{iF}{iM}(iR)}] };
                roi_hms_hld = { 'lh' 'rh' };
            end
            
            % Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fcfg = [];
            fcfg.sbj_nme = cmb_cln_dta_sbj;
            fcfg.fle_nme = inp_fle{iF};
            fcfg.mes_nme = inp_mse{iF}{iM};
            fcfg.rcn_nme = rcn_fle;
            fcfg.roi_nme = roi_nme_hld;
            fcfg.roi_hms = roi_hms_hld;
            
            [ neu_bio_dta, neu_bio_wrn ] = ejk_extract_mmps_roi( fcfg );
            neu_bio_dta(1,:) = cellfun(@(x) strrep(x,[inp_mse{iF}{iM} '_'],''),neu_bio_dta(1,:),'uni',0);
            
            % Save
            cell2csv( roi_fle_nme, neu_bio_dta)
            if ~isempty(neu_bio_wrn); cell2csv( [roi_fle_nme(1:end-4) '_warnings.txt' ], neu_bio_wrn); end
            
            % QC
            fcfg = [];
            fcfg.sbj_nme = neu_bio_dta(2:end,1);
            fcfg.dta     = cell2mat(neu_bio_dta(2:end,5:end));
            fcfg.dta_lbl = neu_bio_dta(1,5:end);
            fcfg.out_dir     = qal_dir;
            fcfg.out_pre_fix = inp_pre{iF};
            ejk_qc_roi(fcfg)
            
            % ICV normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmpi(inp_icv{iF}{iM},'')
                
                fcfg = [];
                
                fcfg.dta     = neu_bio_dta;
                fcfg.cor_col = inp_icv{iF}{iM};
                
                neu_bio_dta_icv = ejk_cor_roi( fcfg );
                
                % Save
            cell2csv( [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '.csv'], neu_bio_dta_icv)
                
            end
            
            % Laterality Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if inp_lat{iF}(iM)==1
                
                [ dta_out , dta_lbl ] = ejk_create_laterality_index( neu_bio_dta(2:end,:) , neu_bio_dta(1,:));
                neu_bio_dta_lat = [ neu_bio_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
                
                % Save
                cell2csv( [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '.csv'], neu_bio_dta_lat)
                
                % QC
                fcfg = [];
                fcfg.sbj_nme = neu_bio_dta(2:end,1);
                fcfg.dta     = cell2mat(neu_bio_dta(2:end,5:end));
                fcfg.dta_lbl = neu_bio_dta(1,5:end);
                fcfg.out_dir     = [ qal_dir(1:end-1) '_' 'LateralityIndex' ];
                fcfg.out_pre_fix = inp_pre{iF};
                ejk_qc_roi(fcfg)
                
            end
            
            % Clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear neu_bio_dta neu_bio_dta_icv neu_bio_dta_lat neu_bio_wrn
            
        end
    end
end

%% Usable Neuroimaging Data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'subcort_vol_dev_norm_IntracranialVolume.csv'];
fcfg.dta_col = 2;
[ mri_dev_dta, mri_dev_dta_sbj, mri_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'fiber_FA_dev.csv'];
fcfg.dta_col = 2;
[ dti_fib_dev_dta, dti_fib_dev_dta_sbj, dti_fib_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_dev.csv'];
fcfg.dta_col = 2;
[ dti_wmp_dev_dta, dti_wmp_dev_dta_sbj, dti_wmp_dev_dta_col] = ejk_dta_frm( fcfg );

use_mri_dta     = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(mri_dev_dta_col,'Left_Hippocampus')))) );
use_dti_fib_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_fib_dev_dta_col,'L_IFO')))) );

use_mri_lbl = repmat({'no'},numel(cmb_cln_dta_sbj),1);
    use_mri_lbl(use_mri_dta) = {'yes'};
use_dti_lbl = repmat({'no'},numel(cmb_cln_dta_sbj),1);
    use_dti_lbl(use_dti_fib_dta) = {'yes'};

cmb_cln_dta_col = [cmb_cln_dta_col {'MRI'}     {'DTI'}];
cmb_cln_dta     = [cmb_cln_dta     use_mri_lbl use_dti_lbl ];

%% Groups
% Grouping 1: Site
grp.site.control = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'HC'));
grp.site.emory   = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'Emory'));
grp.site.ucsf    = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'UCSF'));
grp.site.ucsd    = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'UCSD'));

% Grouping 2: presurgical-x-condition
has_pre_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr')))) );

ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));
con_trl     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'HC'));

grp.pre_cog.control = intersect(has_pre_cog,con_trl);
grp.pre_cog.ltle   = intersect(has_pre_cog,ltl_tle);
grp.pre_cog.rtle    = intersect(has_pre_cog,rgh_tle);

% Grouping 3: postsurgical-x-condition
has_pst_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two_pst')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr_pst')))) );

ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));

grp.pst_cog.ltle    = intersect(has_pst_cog,ltl_tle);
grp.pst_cog.rtle    = intersect(has_pst_cog,rgh_tle);

% Grouping 4: presurgical-x-condition-x-imaging
use_mri_dta     = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(mri_dev_dta_col,'Left_Hippocampus')))) );
use_dti_fib_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_fib_dev_dta_col,'L_IFO')))) );
use_dti_wmp_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_wmp_dev_dta_col,'lh_parsorbitalis')))) );

has_pre_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr')))) );

ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));
con_trl     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'HC'));

grp.pre_cog_img.control = intersect(intersect(has_pre_cog,con_trl),use_dti_fib_dta);
grp.pre_cog_img.ltle    = intersect(intersect(has_pre_cog,ltl_tle),use_dti_fib_dta);
grp.pre_cog_img.rtle    = intersect(intersect(has_pre_cog,rgh_tle),use_dti_fib_dta);

% Grouping 5: postsurgical-x-condition-x-imaging
use_mri_dta     = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(mri_dev_dta_col,'Left_Hippocampus')))) );
use_dti_fib_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_fib_dev_dta_col,'L_IFO')))) );
use_dti_wmp_dta = find( ~isnan(cell2mat(mri_dev_dta(:,strcmpi(dti_wmp_dev_dta_col,'lh_parsorbitalis')))) );

has_pst_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two_pst')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr_pst')))) );

ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));

grp.pst_cog_img.ltle    = intersect(intersect(has_pst_cog,ltl_tle),use_dti_fib_dta);
grp.pst_cog_img.rtle    = intersect(intersect(has_pst_cog,rgh_tle),use_dti_fib_dta);

% Grouping 6: SurgeryType
ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));

grp.surgery.ltle_slah = intersect(ltl_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'SLAH')));
grp.surgery.ltle_atl  = intersect(ltl_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'ATL')));
grp.surgery.rtle_slah = intersect(rgh_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'SLAH')));
grp.surgery.rtle_atl  = intersect(rgh_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'ATL')));

% Grouping 7: WMS version - explore
grp.wms_ver.lm2_two   = intersect(find(strcmpi(wms_hld(:,2),'II')),grp.site.emory);
grp.wms_ver.lm2_three = intersect(find(strcmpi(wms_hld(:,2),'III')),grp.site.emory);
grp.wms_ver.lm2_four  = intersect(find(strcmpi(wms_hld(:,2),'IV')),grp.site.emory);
grp.wms_ver.vp2_two   = intersect(find(strcmpi(wms_hld(:,3),'II')),grp.site.emory);
grp.wms_ver.vp2_three = intersect(find(strcmpi(wms_hld(:,3),'III')),grp.site.emory);
grp.wms_ver.vp2_four  = intersect(find(strcmpi(wms_hld(:,3),'IV')),grp.site.emory);

% Grouping 8: L-TLE only
ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));

grp.surgery_ltle.ltle_slah = intersect(ltl_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'SLAH')));
grp.surgery_ltle.ltle_atl  = intersect(ltl_tle,find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'ATL')));

% Grouping 9: Surgery-Only
grp.surgery_only.ltle_slah = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'SLAH'));
grp.surgery_only.ltle_atl  = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SurgeryType')),'ATL'));

%% New cognitive data
spr_cmb_cog_dta     = cmb_cog_dta(:,ismember(cmb_cog_dta_col,{'log_mem_nor_scr_two' 'vp2_nor_scr' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst'}));
spr_cmb_cog_dta_col = cmb_cog_dta_col(:,ismember(cmb_cog_dta_col,{'log_mem_nor_scr_two' 'vp2_nor_scr' 'log_mem_nor_scr_two_pst' 'vp2_nor_scr_pst'}));

emy_sbj = find(strcmpi( cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location'),:),'Emory'));
for iS = 1:numel(emy_sbj)
    spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'log_mem_nor_scr_two_pst')} = ...
        ( ( spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'log_mem_nor_scr_two')} - ...
            spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'log_mem_nor_scr_two_pst')}) - ...
            (12.5-10.2) ) / 2.07;
        
        
    spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'vp2_nor_scr_pst')} = ...
        ( ( spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'vp2_nor_scr')} - ...
            spr_cmb_cog_dta{emy_sbj(iS),strcmpi(spr_cmb_cog_dta_col,'vp2_nor_scr_pst')}) - ...
            (11.1-10.5) ) / 1.88;
end


has_pre_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr')))) );
has_pst_cog = find(~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'log_mem_nor_scr_two_pst')))) | ...
                   ~isnan(cell2mat(cmb_cog_dta(:,strcmpi(cmb_cog_dta_col,'vp2_nor_scr_pst')))) );

use_pre_cog_lbl = repmat({'no'},numel(cmb_cln_dta_sbj),1);
    use_pre_cog_lbl(has_pre_cog) = {'yes'};
use_pst_cog_lbl = repmat({'no'},numel(cmb_cln_dta_sbj),1);
    use_pst_cog_lbl(has_pst_cog) = {'yes'};

cmb_cln_dta_col = [cmb_cln_dta_col {'PreCog'}      {'PstCog'}];
cmb_cln_dta     = [cmb_cln_dta     use_pre_cog_lbl use_pst_cog_lbl ];

%% New cognitive data - NO RCI's
% Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lgm_pre_raw = { 'lm2_raw'       'LM II raw_pre' };
lgm_pre_scl = { 'lm2_ss'        'LM II norm_ss_pre' };
vpa_pre_raw = { 'vpa2_raw'      'VPA II raw_pre' };
vpa_pre_scl = { 'vpa2_ss'       'VPA II norm_ss_pre' };
lgm_pst_raw = { 'post_lm2_raw'  'LM II raw_6m.post' };
lgm_pst_scl = { 'post_lm2_ss'   'LM II norm_ss_6m.post' };
vpa_pst_raw = { 'post_vpa2_raw' 'VPA II raw_1y.post' };
vpa_pst_scl = { 'post_vpa2_ss'  'VPA II norm_ss__1y.post' };

new_cog_dta_sbj = cmb_cln_dta_sbj;
new_cog_dta_col = { 'lm2_pre_raw' 'lm2_pre_scl' 'vp2_pre_raw' 'vp2_pre_scl' ...
                    'lm2_pst_raw' 'lm2_pst_scl' 'vp2_pst_raw' 'vp2_pst_scl' ...
                    'lm2_chg_raw' 'lm2_chg_scl' 'vp2_chg_raw' 'vp2_chg_scl' ...
                    'lm2_pct_raw' 'lm2_pct_scl' 'vp2_pct_raw' 'vp2_pct_scl' };
new_cog_dta     = cell(numel(new_cog_dta_sbj),numel(new_cog_dta_col));

for iS = 1:numel(cmb_cln_dta_sbj)
    if any(ismember(grp.site.ucsd,iS)) || any(ismember(grp.site.ucsf,iS)) || any(ismember(grp.site.control,iS))
        new_cog_dta{iS,1} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,lgm_pre_raw{1})};
        new_cog_dta{iS,2} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,lgm_pre_scl{1})};
        new_cog_dta{iS,3} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,vpa_pre_raw{1})};
        new_cog_dta{iS,4} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,vpa_pre_scl{1})};
        new_cog_dta{iS,5} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,lgm_pst_raw{1})};
        new_cog_dta{iS,6} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,lgm_pst_scl{1})};
        new_cog_dta{iS,7} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,vpa_pst_raw{1})};
        new_cog_dta{iS,8} = red_cap_dta{strcmpi(red_cap_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(red_cap_dta_col,vpa_pst_scl{1})};
    elseif any(ismember(grp.site.emory,iS))
        new_cog_dta{iS,1} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,lgm_pre_raw{2})};
        new_cog_dta{iS,2} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,lgm_pre_scl{2})};
        new_cog_dta{iS,3} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,vpa_pre_raw{2})};
        new_cog_dta{iS,4} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,vpa_pre_scl{2})};
        new_cog_dta{iS,5} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,lgm_pst_raw{2})};
        new_cog_dta{iS,6} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,lgm_pst_scl{2})};
        new_cog_dta{iS,7} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,vpa_pst_raw{2})};
        new_cog_dta{iS,8} = emy_dta{strcmpi(emy_dta_sbj,cmb_cln_dta_sbj{iS}),strcmpi(emy_dta_col,vpa_pst_scl{2})};        
    end   
end

new_cog_dta(cellfun(@isempty,new_cog_dta)) = {NaN};
new_cog_dta(cellfun(@ischar,new_cog_dta)) = {NaN};
new_cog_dta = cell2mat(new_cog_dta);

for iS = 1:numel(cmb_cln_dta_sbj)
    for iC = 9:12
        new_cog_dta(iS,iC) = new_cog_dta(iS,iC-4) - new_cog_dta(iS,iC-8);
    end
    for iC = 13:16
        new_cog_dta(iS,iC) = round(((new_cog_dta(iS,iC-8) - new_cog_dta(iS,iC-12)) / new_cog_dta(iS,iC-12))*100);        
    end
end

new_cog_dta(isinf(new_cog_dta)) = NaN;

%% Surgical Explore
srg_tab = tabulate(cmb_cln_dta(:,10));
srg_tab(strcmpi(srg_tab(:,1),'SLAH'),:) = [];
srg_tab(strcmpi(srg_tab(:,1),'ATL'),:) = [];

tbl_out = cell(sum(cell2mat(srg_tab(:,2))),2);
tbl_cnt = 1;
for iSU = 1:size(srg_tab,1)
    sbj_ind = find(strcmpi(cmb_cln_dta(:,10),srg_tab{iSU,1}));
    for iSB = 1:numel(sbj_ind)
        tbl_out{tbl_cnt,1} = srg_tab{iSU,1};
        tbl_out{tbl_cnt,2} = cmb_cln_dta_sbj{sbj_ind(iSB)};
        tbl_cnt = tbl_cnt+1;
    end    
end

cell2csv([ out_dir '/' 'complicated_surgeries.csv'],tbl_out); clear tbl_out

%% Cognitive Stats
tbl_typ = { 'mean/std' new_cog_dta_col{1}  ''         1 ; ...
            'mean/std' new_cog_dta_col{2}  ''         1 ; ...
            'mean/std' new_cog_dta_col{3}  ''         1 ; ...
            'mean/std' new_cog_dta_col{4}  ''         1 ;
            'mean/std' new_cog_dta_col{5}  ''         1 ; ...
            'mean/std' new_cog_dta_col{6}  ''         1 ; ...
            'mean/std' new_cog_dta_col{7}  ''         1 ; ...
            'mean/std' new_cog_dta_col{8}  ''         1 ;
            'mean/std' new_cog_dta_col{9}  ''         1 ; ...
            'mean/std' new_cog_dta_col{10} ''         1 ; ...
            'mean/std' new_cog_dta_col{11} ''         1 ; ...
            'mean/std' new_cog_dta_col{12} ''         1 ;
            'mean/std' new_cog_dta_col{13} ''         1 ; ...
            'mean/std' new_cog_dta_col{14} ''         1 ; ...
            'mean/std' new_cog_dta_col{15} ''         1 ; ...
            'mean/std' new_cog_dta_col{16} ''         1 };

grp_nme = { 'site' 'pre_cog' 'pst_cog' 'pre_cog_img' 'pst_cog_img' 'surgery' 'wms_ver' 'surgery_ltle' 'surgery_only' };

new_cog_dta_use = num2cell(new_cog_dta);

for iG = 1:numel(grp_nme)
       
    % get group data
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    fcfg = [];
    fcfg.grp     = grp.(grp_nme{iG});
    fcfg.grp_inc = {fld_nme};
    fcfg.grp_nme = {fld_nme};
    fcfg.dta = new_cog_dta_use(:,cellfun(@isnumeric,new_cog_dta_use(1,:)));
    fcfg.sbj = new_cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    if numel(fld_nme)>2
        % ANOVA
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = new_cog_dta_col(cellfun(@isnumeric,new_cog_dta_use(1,:)));
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive' '/' 'ANOVA' '/' ];
        ejk_1way_anova( fcfg )
    elseif numel(fld_nme)==2
        % t-test
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = new_cog_dta_col(cellfun(@isnumeric,new_cog_dta_use(1,:)));
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive' '/' 'ttest2' '/' ];
        ejk_ttest2_independent( fcfg );
        
    end
    
end

%% Cognitive Table
% Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dta_inp{1} = [ new_cog_dta_col ; num2cell(new_cog_dta) ];
      
grp_nme = { 'site' 'pre_cog' 'pst_cog' 'pre_cog_img' 'pst_cog_img' 'surgery' 'wms_ver' 'surgery_ltle' 'surgery_only' };
      
for iG = 1:numel(grp_nme)
       
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    if numel(fld_nme)==2
        dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive' '/' 'ttest2' '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
    elseif numel(fld_nme)>2
        dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive' '/' 'ANOVA' '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
    end
    
    fcfg = [];
    
    fcfg.tbl = cell(size(tbl_typ,1),numel(fld_nme));
    for iR = 1:size(tbl_typ,1)
        for iN = 1:numel(fld_nme)            
            fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
        end
        fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
    end
    
    fcfg.dta = dta_inp;
    fcfg.grp = grp.(grp_nme{iG});
    tbl_out  = ejk_create_table( fcfg );
    
    num_sbj{1} = 'N';
    tbl_lbl{1} = '';
    for iN = 1:numel(fld_nme)
        num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(fld_nme{iN}));
        tbl_lbl{iN+1} = fld_nme{iN};
    end
    num_sbj{iN+2} = '-'; 
    tbl_lbl{iN+2} = 'Stats';
    
    tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
    
    cell2csv([ out_dir '/' 'wms_explore' '/' 'cognitive_table_' grp_nme{iG} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
end
15 14

% Scatterplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    ejk_chk_dir([ out_dir '/' 'scatters' '/' grp_nme{iG} '/']);
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    for iR = 1:size(tbl_typ,1)
        
        fcfg = [];
        
        fcfg.xdt     = num2cell(1:numel(fld_nme));
        for iN = 1:numel(fld_nme)
            fcfg.ydt{iN} = cell2mat(dta_inp{tbl_typ{iR,4}}(grp.(grp_nme{iG}).(fld_nme{iN})+1 ,strcmpi(dta_inp{tbl_typ{iR,4}}(1,:),tbl_typ{iR,2})));      
            fcfg.ydt{iN}(fcfg.ydt{iN}>100)    = 100;
            fcfg.ydt{iN}(fcfg.ydt{iN}<-100)   = -100;
        end
        
        fcfg.fce_col = { rgb('black') rgb('slate blue') rgb('brick red') rgb('orange') rgb('light teal') rgb('mauve')};
        fcfg.edg_col = { rgb('black') rgb('black')      rgb('black')     rgb('black')  rgb('black')      rgb('black') };
        
        if contains(tbl_typ{iR,2},'_pct_')
           fcfg.ylm = [-100 100]; 
        end
        
        fcfg.xlb = cellfun(@(x) mmil_spec_char(x,{'_'},{' '}),fld_nme,'UniformOutput',false)';
        fcfg.ylb = {mmil_spec_char(tbl_typ{iR,2},{'_'},{' '})};
        
        fcfg.out_dir = [ out_dir '/' 'wms_explore' '/' 'scatters' '/' grp_nme{iG} '/'];
        fcfg.out_nme = tbl_typ{iR,2};
        
        ejk_scatter(fcfg)
    end
    
end

clear dta_inp

%% Clinical Stats
grp_nme = { 'site' 'pre_cog' 'pst_cog' 'pre_cog_img' 'pst_cog_img' 'surgery' 'wms_ver' 'surgery_ltle' 'surgery_only' };

tbl_typ = { 'count'    'Location'                'HC/UCSD/UCSF/Emory' 1; ...
            'mean/std' 'AgeAtSurgery'            ''         1 ; ...
            'mean/std' 'Educ'                    ''         1 ; ...
            'count'    'Sex'                     'M/F'      1 ; ...
            'count'    'Handedness'              'R/L'      1 ; ...
            'mean/std' 'AgeOfSeizureOnset'       ''         1 ; ...
            'mean/std' 'NumAEDs'                 ''         1 ; ...
            'count'    'SideOfSeizureFocus'      'R/L'      1 ; ...
            'mean/std' 'SeizureFreq'             ''         1 ; ...
            'count'    'SurgeryType'             'ATL/SLAH' 1 ; ...
            'count'    'MRI'                     'yes/no'   1 ; ...
            'count'    'DTI'                     'yes/no'   1 ;
            'count'    'PstCog'                  'yes/no'   1 };
        
[~, use_dta_col, use_tbl_col ] = intersect( cmb_cln_dta_col, tbl_typ(:,2) );
use_dta = cmb_cln_dta(:,use_dta_col);
use_col = cmb_cln_dta_col(use_dta_col);
use_lvl = tbl_typ(use_tbl_col,3);

for iG = 1:numel(grp_nme)
       
    % get numeric group data
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    fcfg = [];
    fcfg.grp     = grp.(grp_nme{iG});
    fcfg.grp_inc = {fld_nme};
    fcfg.grp_nme = {fld_nme};
    fcfg.dta = use_dta(:,cellfun(@isnumeric,use_dta(1,:)));
    fcfg.sbj = cmb_cln_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    if numel(fld_nme)>2
        % ANOVA
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = use_col(cellfun(@isnumeric,use_dta(1,:)));
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical' '/' 'ANOVA' '/' ];
        ejk_1way_anova( fcfg )
    elseif numel(fld_nme)==2
        % t-test
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = use_col(cellfun(@isnumeric,use_dta(1,:)));
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical' '/' 'ttest2' '/' ];
        ejk_ttest2_independent( fcfg );        
    end
    
    % get factor group data
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    fcfg = [];
    fcfg.grp     = grp.(grp_nme{iG});
    fcfg.grp_inc = {fld_nme};
    fcfg.grp_nme = {fld_nme};
    fcfg.dta = use_dta(:,~cellfun(@isnumeric,use_dta(1,:)));
    fcfg.sbj = cmb_cln_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    % Fishers
    fcfg = [];    
    fcfg.sbj = grp_sbj{1};
    fcfg.grp_nme = grp_nme{iG};
    fcfg.dta_one = grp_dta{1};
    fcfg.lbl_one = use_col(~cellfun(@isnumeric,use_dta(1,:))); 
    fcfg.lvl     = use_lvl(~cellfun(@isnumeric,use_dta(1,:)));
    fcfg.dta_two = repmat(grp_typ{1},1,sum(~cellfun(@isnumeric,use_dta(1,:))));
    fcfg.lbl_two = strcat( 'group_', use_col(~cellfun(@isnumeric,use_dta(1,:))));
    fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Clinical' '/' 'Fishers' '/' ];    
    ejk_fisher_test( fcfg );
    
end

%% Clinical Table
% Table
dta_inp{1} = [ cmb_cln_dta_col     ; cmb_cln_dta ];

grp_nme = { 'site' 'pre_cog' 'pst_cog' 'pre_cog_img' 'pst_cog_img' 'surgery' 'surgery_ltle' 'surgery_only' };
      
for iG = 1:numel(grp_nme)
       
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    if numel(fld_nme)==2
        dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical' '/' 'ttest2' '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
    elseif numel(fld_nme)>2
        dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical' '/' 'ANOVA' '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv']);
    end
    dta_inp{3} = mmil_readtext([ out_dir '/' 'stats' '/' 'Clinical' '/' 'Fishers' '/' grp_nme{iG} '/' 'output_table.csv']);
    
    fcfg = [];
    
    for iR = 1:size(tbl_typ)
        if strcmpi(tbl_typ{iR,1},'mean/std')
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
            end
            fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
        elseif  strcmpi(tbl_typ{iR,1},'count')
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2} ',' tbl_typ{iR,3}];
            end
            fcfg.tbl{iR,iN+1} = ['copy,3,' tbl_typ{iR,2} ',report'];
        end
    end
    
    fcfg.dta = dta_inp;
    fcfg.grp = grp.(grp_nme{iG});
    tbl_out  = ejk_create_table( fcfg );
    
    num_sbj{1} = 'N';
    tbl_lbl{1} = '';
    for iN = 1:numel(fld_nme)
        num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(fld_nme{iN}));
        tbl_lbl{iN+1} = fld_nme{iN};
    end
    num_sbj{iN+2} = '-'; 
    tbl_lbl{iN+2} = 'Stats';
    
    tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
    
    cell2csv([ out_dir '/' 'wms_explore' '/' 'clinical_table_' grp_nme{iG} '.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
end

% Scatterplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    ejk_chk_dir([ out_dir '/' 'scatters' '/' grp_nme{iG} '/']);
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    for iR = 1:size(tbl_typ,1)
        if strcmpi(tbl_typ{iR,1},'mean/std')
        fcfg = [];
        
        fcfg.xdt     = num2cell(1:numel(fld_nme));
        for iN = 1:numel(fld_nme)
            fcfg.ydt{iN} = cell2mat(dta_inp{tbl_typ{iR,4}}(grp.(grp_nme{iG}).(fld_nme{iN})+1 ,strcmpi(dta_inp{tbl_typ{iR,4}}(1,:),tbl_typ{iR,2})));      
        end
        
        fcfg.fce_col = { rgb('black') rgb('slate blue') rgb('brick red') rgb('orange') rgb('light teal') rgb('mauve')};
        fcfg.edg_col = { rgb('black') rgb('black')      rgb('black')     rgb('black')  rgb('black')      rgb('black') };
                
        fcfg.xlb = cellfun(@(x) mmil_spec_char(x,{'_'},{' '}),fld_nme,'UniformOutput',false)';
        fcfg.ylb = {mmil_spec_char(tbl_typ{iR,2},{'_'},{' '})};
        
        fcfg.out_dir = [ out_dir '/' 'wms_explore' '/' 'scatters' '/' grp_nme{iG} '/'];
        fcfg.out_nme = tbl_typ{iR,2};
        
        ejk_scatter(fcfg)
        end
    end
    
end

clear dta_inp

%% Harmonization
% Load Laterality ROI data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_dev_LateralityIndex.csv'];
fcfg.dta_col = 2;
[ dti_wmp_lat_dev_dta, dti_wmp_lat_dev_dta_sbj, dti_wmp_lat_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'subcort_vol_dev_LateralityIndex.csv'];
fcfg.dta_col = 2;
[ vol_lat_dev_dta, vol_lat_dev_dta_sbj, vol_lat_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'fiber_FA_dev_LateralityIndex.csv'];
fcfg.dta_col = 2;
[ dti_fib_lat_dev_dta, dti_fib_lat_dev_dta_sbj, dti_fib_lat_dev_dta_col] = ejk_dta_frm( fcfg );

% Find relevant LI ROIs
neu_use = { {'entorhinal' 'parahippocampal' 'lateralorbitofrontal' 'medialorbitofrontal' } ...
            {'Hippocampus' 'Thalamus_Proper'} ...
            {'Unc' 'IFO' 'ILF'} };

neu_use_lat_dta     = [ dti_wmp_lat_dev_dta(:,ismember(dti_wmp_lat_dev_dta_col,neu_use{1})) ...
                        vol_lat_dev_dta(:,ismember(vol_lat_dev_dta_col,neu_use{2})) ...
                        dti_fib_lat_dev_dta(:,ismember(dti_fib_lat_dev_dta_col,neu_use{3})) ];
neu_use_lat_dta_col = [ dti_wmp_lat_dev_dta_col(:,ismember(dti_wmp_lat_dev_dta_col,neu_use{1})) ...
                        vol_lat_dev_dta_col(:,ismember(vol_lat_dev_dta_col,neu_use{2})) ...
                        dti_fib_lat_dev_dta_col(:,ismember(dti_fib_lat_dev_dta_col,neu_use{3})) ];
neu_use_lat_dta_sbj = dti_wmp_lat_dev_dta_sbj;

% Load Raw ROI data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'wmparc_FA_wm_aparc_annot_dev.csv'];
fcfg.dta_col = 2;
[ dti_wmp_dev_dta, dti_wmp_dev_dta_sbj, dti_wmp_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'subcort_vol_dev.csv'];
fcfg.dta_col = 2;
[ vol_dev_dta, vol_dev_dta_sbj, vol_dev_dta_col] = ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data' '/' 'fiber_FA_dev.csv'];
fcfg.dta_col = 2;
[ dti_fib_dev_dta, dti_fib_dev_dta_sbj, dti_fib_dev_dta_col] = ejk_dta_frm( fcfg );
     
% Find relevant raw ROIs
neu_use = { {'entorhinal' 'parahippocampal' 'lateralorbitofrontal' 'medialorbitofrontal' } ...
            {'Hippocampus' 'Thalamus_Proper'} ...
            {'Unc' 'IFO' 'ILF'} };

neu_use_raw_dta     = [ dti_wmp_dev_dta(:,string_find(dti_wmp_dev_dta_col,neu_use{1})) ...
                        vol_dev_dta(:,string_find(vol_dev_dta_col,neu_use{2})) ...
                        dti_fib_dev_dta(:,string_find(dti_fib_dev_dta_col,neu_use{3})) ];
neu_use_raw_dta_col = [ dti_wmp_dev_dta_col(:,string_find(dti_wmp_dev_dta_col,neu_use{1})) ...
                        vol_dev_dta_col(:,string_find(vol_dev_dta_col,neu_use{2})) ...
                        dti_fib_dev_dta_col(:,string_find(dti_fib_dev_dta_col,neu_use{3})) ];
neu_use_raw_dta_sbj = dti_wmp_dev_dta_sbj;

% Batch
grp.site.control = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'HC'));
grp.site.emory   = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'Emory'));
grp.site.ucsf    = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'UCSF'));
grp.site.ucsd    = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'UCSD'));

btc_dta = cell(size(neu_use_raw_dta,1),1);
btc_dta(grp.site.control) = {'UCSD'};
btc_dta(grp.site.emory) = {'Emory'};
btc_dta(grp.site.ucsf) = {'UCSF'};
btc_dta(grp.site.ucsd) = {'UCSD'};

% Covariates
ltl_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'L'));
rgh_tle     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'SideOfSeizureFocus')),'R'));
con_trl     = find(strcmpi(cmb_cln_dta(:,strcmpi(cmb_cln_dta_col,'Location')),'HC'));

cov_dta = cell(size(neu_use_raw_dta,1),1);
cov_dta(ltl_tle) = {'LTLE'};
cov_dta(rgh_tle) = {'RTLE'};
cov_dta(con_trl) = {'HC'};
cov_dta(cellfun(@isempty,cov_dta)) = {'BTLE'};

% Run
fcfg = [];

fcfg.sbj_nme = neu_use_lat_dta_sbj;

fcfg.dta     = cell2mat([ neu_use_raw_dta neu_use_lat_dta  ]);
fcfg.dta_nme = [ neu_use_raw_dta_col neu_use_lat_dta_col  ];

fcfg.btc     = btc_dta;
fcfg.btc_nme = {'Site'};

fcfg.cov     = cov_dta;
fcfg.cov_nme = {'Diagnosis'};

fcfg.plt = 1;
% fcfg.ylm = [-1 1];

fcfg.out_dir = [ prj_dir '/' prj_nme '/' 'EmoryExplore' '/' 'Data'];

com_bat_epd_dta = ejk_ComBat(fcfg);

% Save
use_dta     = [ com_bat_epd_dta(1:numel(neu_use_raw_dta_col)) ];
use_dta_col = [ neu_use_raw_dta_col neu_use_lat_dta_col  ];
use_dta_sbj = neu_use_lat_dta_sbj;

%% QC

%% Neurobio Correlations
[]

%%

