clear; clc;

prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/';

dta_dir = [ prj_dir '/' 'data' '/' ];

%% Load Data %%%%%%%%%%%%%%%%%%%%%
% Overall data file
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'cov.csv'];
[ cov_dta, cov_dta_sbj, cov_dta_col ] = ejk_dta_frm(fcfg);
    cov_dta_sbj( cellfun(@isnumeric,cov_dta_sbj)) = cellfun(@num2str,cov_dta_sbj( cellfun(@isnumeric,cov_dta_sbj)),'uni',0);
    
% Rush Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'rush_cov.csv'];
[ rsh_cov_dta, rsh_cov_dta_sbj, rsh_cov_dta_col ] = ejk_dta_frm(fcfg);
    rsh_cov_dta_sbj( cellfun(@isnumeric,rsh_cov_dta_sbj)) = cellfun(@num2str,rsh_cov_dta_sbj( cellfun(@isnumeric,rsh_cov_dta_sbj)),'uni',0);
   
% UCSF/UCSD/Emory Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Clinical.csv'];
[ thr_cov_dta, thr_cov_dta_sbj, thr_cov_dta_col ] = ejk_dta_frm(fcfg);

% Emory Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'emy_cov.csv'];
[ emy_cov_dta, emy_cov_dta_sbj, emy_cov_dta_col ] = ejk_dta_frm(fcfg);

% CCF Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'ccf_cov.csv'];
[ ccf_cov_dta, ccf_cov_dta_sbj, ccf_cov_dta_col ] = ejk_dta_frm(fcfg);

% MUSC/Bonn Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'msc_cov.csv'];
[ msc_cov_dta, msc_cov_dta_sbj, msc_cov_dta_col ] = ejk_dta_frm(fcfg);

% Data Column Names %%%%%%%%%%%%%%%%%%%%%
          % Original       % Rush         % UCSD/UCSF/Emory    % Emory        % CCF             % MUSC/Bonn                                      % Type
col_hld = { 'Site'           ''             ''                   ''           ''                ''                                               'str' ; ...
            'Side'           'SDx'          'SideOfSeizureFocus' 'side'       'side'            'SideL0R1'                                       'str' ; ...
            'Age'            'Age'          'AgeAtSurgery'       'age'        'age'             'Age'                                            'num' ; ...
            'Sex'            'Sex'          'Sex'                'sex'        'sex'             'Sex'                                            'str' ; ...
            'Handedness'     'Handedness'   'Handedness'         'handedness' 'handedness'      'Handedness'                                     'str' ; ...
            'Education'      'Education'    'Educ'               'education'  'education'       'Education'                                      'num' ; ...
            'AO'             'AO'           'AgeOfSeizureOnset'  'onset'      'onset'           'Age of Onset'                                   'num' ; ...
            'DURILL'         'DURILL'       ''                   'duration'   'duration'        'Duration of seizures (years)'                   'num' ; ...
            'MTS'            'mts'          'MTS'                'mts'        'mts'             'Mesial Temporal Sclerosis'                      'str' ; ...
            'ASMs'           'ASMs'         'NumAEDs'            'AEDs'       'AEDs'            'Number of Anti-Seizure Meds at time of imaging' 'num' ; ...
            'Surgery'        'Surgery Type' 'SurgeryType'        ''           'surgery'         'SurgeryLas0Res1'                                'str' ; ...
            'Engel'          'ENGEL'        'EngelOutcome'       ''           'surgery_outcome' 'Engel'                                          'str' };


% Check columns
if ~all(strcmpi( col_hld(:,1)', cov_dta_col)); error('check column names'); end

%% Tighten up UCSD/UCSF/Emory data            
clear fix_str

fix_str{find(strcmpi(thr_cov_dta_col,'MTS'))}        =  { 'NonMTS'   'no' ; ...
                                                          'N/A'   'no' ; ...
                                                          'MTS'     'yes' ; ...
                                                          'L'     'yes' ; ...
                                                          'R'     'yes' ; ...
                                                          'empty' 'no'};

fix_str{find(strcmpi(thr_cov_dta_col,'SurgeryType'))} =  { 'ATL'                     'yes' ; ...
                                                           'ATL +'                   'yes' ; ...
                                                           'amygdalohippocampectomy' 'yes' ; ...
                                                           'Left SAH (open)' 'yes' ; ...
                                                           'Left RNS' 'yes' ; ...
                                                           'Right ATL (limited open)' 'yes' ; ...
                                                           'right SLAH; right ATL; right FL resection' 'yes' ; ...
                                                           'Left SLAH with repeat SLAH for more complete ablation' 'yes' ; ...
                                                           'Left SLAH- pt. suffered a surgical complication involving a bleed' 'yes' ; ...
                                                           'Left SLA involving encephalocele in left temporal lobe' 'yes' ; ...
                                                           'right superiror temporal gyrus RF ablation' 'yes' ; ...
                                                           'right SLAH (temporal pole extending amygdla; sparring hippocampus)' 'yes' ; ...
                                                           'left Temporal lobe RNS'  'yes' ; ...
                                                           'Right TL open resection (not standard)' 'yes' ; ...
                                                           'Two right SLAH procedures (memory declined after 2nd - took broader medial TL structures' 'yes' ; ...
                                                           'lesionectomy'            'yes' ; ...
                                                           'SLAH'                    'yes' ; ...
                                                           'no surgery-just had SEEG - still being considered for surgery'   'no'  ; ...
                                                           'no surgery up to date'   'no'  ; ...
                                                           'empty'                   'no'}; 
                                                     
                                                     
for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        for iT = 1:size(fix_str{iC},1)
            if strcmpi(fix_str{iC}{iT,1},'empty')
                thr_cov_dta(cellfun(@isempty,thr_cov_dta(:,iC)),iC) = {fix_str{iC}{iT,2}};
            else
                thr_cov_dta(strcmpi(thr_cov_dta(:,iC),fix_str{iC}{iT,1}),iC) = {fix_str{iC}{iT,2}};
            end
        end
    end
end

%% Tighten up Rush data
% Add in new data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'missing_data' '/' 'additional_data' '/' 'Rush_data_frame_5_5_2022.csv'];
[ new_rsh_cov_dta, new_rsh_cov_dta_sbj, new_rsh_cov_dta_col ] = ejk_dta_frm(fcfg);
    new_rsh_cov_dta_sbj( cellfun(@isnumeric,new_rsh_cov_dta_sbj)) = cellfun(@num2str,new_rsh_cov_dta_sbj( cellfun(@isnumeric,new_rsh_cov_dta_sbj)),'uni',0);

new_fll_nme = { 'Education' 'ASMs' };
old_fll_nme = { 'Education' 'ASMs' };
for iS = 1:numel(rsh_cov_dta_sbj)
    new_ind = strcmpi(new_rsh_cov_dta_sbj,rsh_cov_dta_sbj{iS});
    if sum(new_ind)==1
        for iC = 1:numel(new_fll_nme)
            rsh_cov_dta{iS,strcmpi(rsh_cov_dta_col,old_fll_nme{iC})} = new_rsh_cov_dta{new_ind,strcmpi(new_rsh_cov_dta_col,new_fll_nme{iC})};
        end
    end
end

%
clear fix_str

fix_str{find(strcmpi(rsh_cov_dta_col,'SDx'))}      = { 3 'L' ; ...
                                                       5 'L' ; ...
                                                       4 'R' ; ...
                                                       6 'R' };

fix_str{find(strcmpi(rsh_cov_dta_col,'Sex'))}        = { 1 'M' ; ...
                                                         2 'F' };
                                                       
fix_str{find(strcmpi(rsh_cov_dta_col,'Handedness'))} = { 1 'L' ; ...
                                                         2 'R' };

fix_str{find(strcmpi(rsh_cov_dta_col,'mts'))}        = { 3 'yes' ; ...
                                                         4 'yes' ; ...
                                                         5 'no' ; ...
                                                         6 'no' };

fix_str{find(strcmpi(rsh_cov_dta_col,'Engel'))}      = { 1 'I' ; ...
                                                         2 'II' ; ...
                                                         3 'III'  ; ...
                                                         4 'IV'  ; ...
                                                         5 ''    };
                                                     
for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        hld_num = rsh_cov_dta(:,iC);
        hld_num(cellfun(@isempty,hld_num))={NaN};
        hld_num = cell2mat(hld_num);
        for iT = 1:size(fix_str{iC},1)            
            rsh_cov_dta(hld_num==fix_str{iC}{iT,1},iC) = {fix_str{iC}{iT,2}};
        end
    end
end
  
clear fix_str
fix_str{find(strcmpi(rsh_cov_dta_col,'Surgery Type'))} = { 'Resection' 'yes' ; ...
                                                           'LITT'      'yes' ; ...
                                                           'empty'     'no'  };

for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        for iT = 1:size(fix_str{iC},1)
            if strcmpi(fix_str{iC}{iT,1},'empty')
                rsh_cov_dta(cellfun(@isempty,rsh_cov_dta(:,iC)),iC) = {fix_str{iC}{iT,2}};
            else
                rsh_cov_dta(strcmpi(rsh_cov_dta(:,iC),fix_str{iC}{iT,1}),iC) = {fix_str{iC}{iT,2}};
            end
        end
    end
end

%% Tighten up CCF data            
% Add in new data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'missing_data' '/' 'additional_data' '/' 'CNN Paper Update_4-25-2022.csv'];
[ new_ccf_cov_dta, new_ccf_cov_dta_sbj, new_ccf_cov_dta_col ] = ejk_dta_frm(fcfg);
    new_ccf_cov_dta_sbj = new_ccf_cov_dta(:,1);

new_fll_nme = { 'Patient Type' 'Seizure Free at 6 Month Follow-Up?' };
old_fll_nme = { 'surgery'      'surgery_outcome' };
for iS = 1:numel(ccf_cov_dta_sbj)
    new_ind = strcmpi(new_ccf_cov_dta_sbj,ccf_cov_dta_sbj{iS});
    if sum(new_ind)==1
        for iC = 1:numel(new_fll_nme)
            ccf_cov_dta{iS,strcmpi(ccf_cov_dta_col,old_fll_nme{iC})} = new_ccf_cov_dta{new_ind,strcmpi(new_ccf_cov_dta_col,new_fll_nme{iC})};
        end
    end
end

%
clear fix_str

fix_str{find(strcmpi(ccf_cov_dta_col,'side'))}      = { 'Left'  'L' ; ...
                                                        'Right' 'R' };

fix_str{find(strcmpi(ccf_cov_dta_col,'sex'))}        = { 'Male'   'M' ; ...
                                                         'Female' 'F' };
                                                       
fix_str{find(strcmpi(ccf_cov_dta_col,'handedness'))} = { 'Left'  'L' ; ...
                                                         'Right' 'R' ; ...
                                                         'Ambi'  'L'};

fix_str{find(strcmpi(ccf_cov_dta_col,'mts'))}        = { 'MTS'    'yes' ; ...
                                                         'NonMTS' 'no'};

fix_str{find(strcmpi(ccf_cov_dta_col,'surgery'))}    = { 'Surgical' 'yes' ; ...
                                                         'Nonsurgical - EEG Lateralized Left Temporal' 'no'};

fix_str{find(strcmpi(ccf_cov_dta_col,'surgery_outcome'))}  = { 'Yes'    'I' ; ...
                                                               'Yes (1 year followup)' 'I' ; ...
                                                               'Yes (4 year followup)' 'I' ; ...
                                                               'seizure 2 years postop after discontinuing ASMs; 3 years postop questionable seizure activity; followed locally (nonCCF) note from phone call' 'I' ; ...
                                                               'Non-reoccuring seizures (e.g. post op seizures; seizure off meds)' 'other' ; ...                                                               
                                                               '1 seizure 9 months after surgery in context of illness; seizure free otherwise (2 years post op)' 'other' ; ...
                                                               '1 questionable event 6 months after surgery in context of illness; seizure free otherwise (10 months postop)' 'other' ; ...
                                                               '2 questionable events since surgery (4 months postop)' 'other' ; ...
                                                               'Seizure free 2 months after postop; pt deceased and no further followup' 'other' ; ...
                                                               'Deceased (3 days postop)' 'other'};
                                                     
for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        for iT = 1:size(fix_str{iC},1)
            ccf_cov_dta(strcmpi(ccf_cov_dta(:,iC),fix_str{iC}{iT,1}),iC) = {fix_str{iC}{iT,2}};
        end
    end
end

%% Tighten up Emory data
clear fix_str

fix_str{find(strcmpi(emy_cov_dta_col,'side'))}      = { 'Left'  'L' ; ...
                                                        'Right' 'R' };

fix_str{find(strcmpi(emy_cov_dta_col,'sex'))}        = { 'Male'   'M' ; ...
                                                         'Female' 'F' };
                                                       
fix_str{find(strcmpi(emy_cov_dta_col,'handedness'))} = { 'Left'  'L' ; ...
                                                         'Right' 'R' ; ...
                                                         'ambi' 'L' };

fix_str{find(strcmpi(emy_cov_dta_col,'mts'))}        = { 'MTS'    'yes' ; ...
                                                         'NonMTS' 'no'};

                                                     
for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        for iT = 1:size(fix_str{iC},1)
            emy_cov_dta(strcmpi(emy_cov_dta(:,iC),fix_str{iC}{iT,1}),iC) = {fix_str{iC}{iT,2}};
        end
    end
end

%% Tighten up MUSC/Bonn data
% Add in new data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'missing_data' '/' 'additional_data' '/' 'MUSC_demographics_table.csv'];
[ new_msc_cov_dta, new_msc_cov_dta_sbj, new_msc_cov_dta_col ] = ejk_dta_frm(fcfg);
    new_msc_cov_dta(cellfun(@isempty,new_msc_cov_dta(:,strcmpi(new_msc_cov_dta_col,'Mesial Temporal Sclerosis'))),strcmpi(new_msc_cov_dta_col,'Mesial Temporal Sclerosis')) = {0};
    new_msc_cov_dta(cellfun(@isnumeric,new_msc_cov_dta(:,strcmpi(new_msc_cov_dta_col,'Engel'))),strcmpi(new_msc_cov_dta_col,'Engel')) = cellfun(@num2str,new_msc_cov_dta(cellfun(@isnumeric,new_msc_cov_dta(:,strcmpi(new_msc_cov_dta_col,'Engel'))),strcmpi(new_msc_cov_dta_col,'Engel')),'uni',0);
    
new_fll_nme = { 'Mesial Temporal Sclerosis' 'Number of Anti-Seizure Meds at time of imaging' 'Engel' 'Age' 'Sex' 'Handedness' 'Education' 'Age of Onset' 'Duration of seizures (years)' };
old_fll_nme = { 'Mesial Temporal Sclerosis' 'Number of Anti-Seizure Meds at time of imaging' 'Engel' 'Age' 'Sex' 'Handedness' 'Education' 'Age of Onset' 'Duration of seizures (years)' };
for iS = 1:numel(msc_cov_dta_sbj)
    new_ind = strcmpi(new_msc_cov_dta_sbj,msc_cov_dta_sbj{iS});
    if sum(new_ind)==1
        for iC = 1:numel(new_fll_nme)
            msc_cov_dta{iS,strcmpi(msc_cov_dta_col,old_fll_nme{iC})} = new_msc_cov_dta{new_ind,strcmpi(new_msc_cov_dta_col,new_fll_nme{iC})};
        end
    end
end

%
clear fix_str

fix_str{find(strcmpi(msc_cov_dta_col,'SideL0R1'))}      = { 0 'L' ; ...
                                                            1 'R' };

fix_str{find(strcmpi(msc_cov_dta_col,'SurgeryLas0Res1'))}        = { 0 'yes' ; ...
                                                                     1 'yes' };
                                                                 
fix_str{find(strcmpi(msc_cov_dta_col,'Sex'))}        = { 0 'F' ; ...
                                                         1 'M' };  
                                                                    
fix_str{find(strcmpi(msc_cov_dta_col,'Handedness'))} = { 0 'L' ; ...
                                                         1 'R' };  
                                                     
fix_str{find(strcmpi(msc_cov_dta_col,'Mesial Temporal Sclerosis'))} = { 0 'no' ; ...
                                                                        1 'yes' };  

                                                     
for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        hld_num = msc_cov_dta(:,iC);
        hld_num(cellfun(@isempty,hld_num))={NaN};
        hld_num = cell2mat(hld_num);
        for iT = 1:size(fix_str{iC},1)            
            msc_cov_dta(hld_num==fix_str{iC}{iT,1},iC) = {fix_str{iC}{iT,2}};
        end
    end
end

%
clear fix_str
fix_str{find(strcmpi(msc_cov_dta_col,'Engel'))} = { '1'  'I' ; ...
                                                    '1a' 'I' ; ...
                                                    '1b' 'I' ; ...
                                                    '1c' 'I' ; ...
                                                    '1d' 'I' ; ...
                                                    '2'  'II' ; ...
                                                    '2a' 'II' ; ...
                                                    '2b' 'II' ; ...
                                                    '2d' 'II' ; ...
                                                    '3'  'III' ; ...
                                                    '3a' 'III' ; ...
                                                    '4'  'IV' ; ...
                                                    '4b' 'IV'};

for iC = 1:numel(fix_str)
    if ~isempty(fix_str{iC})
        for iT = 1:size(fix_str{iC},1)
            if strcmpi(fix_str{iC}{iT,1},'empty')
                msc_cov_dta(cellfun(@isempty,msc_cov_dta(:,iC)),iC) = {fix_str{iC}{iT,2}};
            else
                msc_cov_dta(strcmpi(msc_cov_dta(:,iC),fix_str{iC}{iT,1}),iC) = {fix_str{iC}{iT,2}};
            end
        end
    end
end

%% Go through
out_dta_sbj = cell(size(cov_dta,1),1);
out_dta_col = col_hld(:,1)';
out_dta     = cell(size(cov_dta));
for iC = 1:size(col_hld,1)
   if strcmpi(col_hld{iC,end},'str')
       out_dta(:,iC) = repmat({''},size(out_dta,1),1);
   end
end

mss_dta     = ones(size(cov_dta));
msc_emy_sbj = cell(0);
for iS = 1:numel(cov_dta_sbj)
    switch cov_dta{iS,1}
        
        case {'UCSD' 'UCSF' 'Emory'}
            sbj_hld = strcmpi(thr_cov_dta_sbj,cov_dta_sbj{iS});
            use_dta = thr_cov_dta;
            use_col = thr_cov_dta_col;
            col_ind = 3;
            
            if sum(sbj_hld)==0 && strcmpi(cov_dta{iS,1},'Emory')
                sbj_hld = strcmpi(emy_cov_dta_sbj,cov_dta_sbj{iS});
                use_dta = emy_cov_dta;
                use_col = emy_cov_dta_col;
                col_ind = 4;
                
                if sum(sbj_hld)==0
                    sbj_hld = strcmpi(msc_cov_dta_sbj,cov_dta_sbj{iS});
                    use_dta = msc_cov_dta;
                    use_col = msc_cov_dta_col;
                    col_ind = 6;
                    
                    msc_emy_sbj = [ msc_emy_sbj cov_dta_sbj(iS) ];
                    
                    if sum(sbj_hld)==0
                        sbj_hld = zeros(1);
                        use_dta = [];
                        use_col = [];
                        col_ind = [];
                    end
                end
            elseif sum(sbj_hld)==0
                sbj_hld = zeros(1);
                use_dta = [];
                use_col = [];
                col_ind = [];
            end
            
        case { 'Rush' }
            sbj_hld = strcmpi(rsh_cov_dta_sbj,cov_dta_sbj{iS});
            use_dta = rsh_cov_dta;
            use_col = rsh_cov_dta_col;
            col_ind = 2;
            
        case { 'CCF' }
            sbj_hld = strcmpi(ccf_cov_dta_sbj,cov_dta_sbj{iS});
            use_dta = ccf_cov_dta;
            use_col = ccf_cov_dta_col;
            col_ind = 5;
            
        case { 'Bonn' 'MUSC' }
            sbj_hld = strcmpi(msc_cov_dta_sbj,cov_dta_sbj{iS});
            use_dta = msc_cov_dta;
            use_col = msc_cov_dta_col;
            col_ind = 6;
    end
    
    out_dta_sbj{iS,1} = cov_dta_sbj{iS};
    out_dta{iS,1} = cov_dta{iS,1};
    for iC = 2:size(col_hld,1)
        if ~(sum(sbj_hld)==0)
            if ~isempty(col_hld{iC,col_ind})
                out_dta{iS,iC} = use_dta{sbj_hld,strcmpi(use_col,col_hld{iC,col_ind})};
            else
                mss_dta(iS,iC) = -999;
            end
        else
            mss_dta(iS,iC) = NaN;
        end
    end
end

%% Get missing data report
mss_sbj.sbj     = [ out_dta_sbj(isnan(mss_dta(:,2))) out_dta(isnan(mss_dta(:,2)),1) ];
mss_sbj_dta.sbj = { '' '' };
for iC = 2:size(col_hld,1)
    mss_sbj.(col_hld{iC,1})     = [ out_dta_sbj(mss_dta(:,iC)==-999) out_dta(mss_dta(:,iC)==-999,1) ];
    mss_sbj_dta.(col_hld{iC,1}) = [ out_dta_sbj(cellfun(@isempty,out_dta(:,iC))) out_dta(cellfun(@isempty,out_dta(:,iC)),1) ];
end

% Make table
tbl_ste = unique(out_dta(:,1));
    tbl_ste_tbl = tabulate(out_dta(:,1));
tbl_col = fieldnames(mss_sbj);

mss_tbl     = cell(numel(tbl_ste),numel(tbl_col));
mss_tbl_dta = cell(numel(tbl_ste),numel(tbl_col));
for iC = 1:numel(tbl_col)
    tbl_hld     = tabulate(mss_sbj.(tbl_col{iC})(:,2));
    tbl_dta_hld = tabulate(mss_sbj_dta.(tbl_col{iC})(:,2));
    
    for iR = 1:numel(tbl_ste)
        %
        tbl_ind = strcmpi(tbl_hld(:,1),tbl_ste{iR});
        tot_ind = strcmpi(tbl_ste_tbl(:,1),tbl_ste{iR});
        if sum(tbl_ind)==0
            mss_tbl{iR,iC} = [ num2str(0) ' (' num2str(0) '%)' ];
        else
            mss_tbl{iR,iC} = [ num2str(tbl_hld{tbl_ind,2}) ' (' num2str(round(tbl_hld{tbl_ind,2}/tbl_ste_tbl{tot_ind,2}*100)) '%)' ];
        end
        
        %
        tbl_ind = strcmpi(tbl_dta_hld(:,1),tbl_ste{iR});
        tot_ind = strcmpi(tbl_ste_tbl(:,1),tbl_ste{iR});
        if sum(tbl_ind)==0
            mss_tbl_dta{iR,iC} = [ num2str(0) ' (' num2str(0) '%)' ];
        else
            mss_tbl_dta{iR,iC} = [ num2str(tbl_dta_hld{tbl_ind,2}) ' (' num2str(round(tbl_dta_hld{tbl_ind,2}/tbl_ste_tbl{tot_ind,2}*100)) '%)' ];
        end
        
    end    
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/missing_report.csv',[ {''} tbl_col'  ; tbl_ste mss_tbl])
cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/missing_report_data.csv',[ {''} tbl_col'  ; tbl_ste mss_tbl_dta])

%% Save out individual tables for sending to sites
ste_nme = unique(out_dta(:,strcmpi(out_dta_col,'Site')));
for iN = 1:numel(ste_nme)
    sve_dta = out_dta(strcmpi(out_dta(:,strcmpi(out_dta_col,'Site')),ste_nme{iN}),:);
    sve_sbj = out_dta_sbj(strcmpi(out_dta(:,strcmpi(out_dta_col,'Site')),ste_nme{iN}),:);
    cell2csv(['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/' ste_nme{iN} '_data_frame.csv'],[ [{''} out_dta_col] ; sve_sbj sve_dta ])
end

% Emory from MUSC
sve_dta = out_dta(strcmpi(out_dta(:,strcmpi(out_dta_col,'Site')),'Emory') & ismember(out_dta_sbj,msc_emy_sbj),:);
sve_sbj = msc_emy_sbj';
cell2csv(['/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/' 'Emory_from_MUSC' '_data_frame.csv'],[ [{''} out_dta_col] ; sve_sbj sve_dta ])

%% Add in stats
grp.L = find(strcmpi(out_dta(:,strcmpi(out_dta_col,'Side')),'L'));
grp.R = find(strcmpi(out_dta(:,strcmpi(out_dta_col,'Side')),'R'));

for iC = 1:size(out_dta,2)
    if isnumeric(out_dta{350,iC})
        out_dta(cellfun(@isempty,out_dta(:,iC)),iC) = {NaN};
    else
        out_dta(cellfun(@isempty,out_dta(:,iC)),iC) = {''};
    end
end

% ttest2 Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_dta_col = cellfun(@isnumeric, out_dta(350,:));

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = {{'L' 'R'}};
fcfg.grp_nme = {{'L' 'R'}};
fcfg.dta = out_dta(:,use_dta_col);
fcfg.sbj = out_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

fcfg = [];
fcfg.sbj_nme = grp_sbj{1};
fcfg.dta     = grp_dta{1};
fcfg.dta_nme = out_dta_col(:,use_dta_col);
fcfg.grp     = grp_typ{1};
fcfg.grp_nme = {'side_ttest2'};
fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/tables/demographics';
ejk_ttest2_independent( fcfg );

% fisher's Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_dta_col = ~cellfun(@isnumeric, out_dta(350,:));

fcfg = [];
fcfg.grp     = grp;
fcfg.grp_inc = {{'L' 'R'}};
fcfg.grp_nme = {{'L' 'R'}};
fcfg.dta = out_dta(:,use_dta_col);
fcfg.sbj = out_dta_sbj;
[ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
grp_dta{1}( cellfun(@isempty,grp_dta{1})) = {NaN};

fcfg = [];
fcfg.sbj = grp_sbj{1};
fcfg.dta_one = grp_dta{1};
fcfg.lbl_one = out_dta_col(:,use_dta_col);
fcfg.dta_two = repmat(grp_typ{1},1,sum(use_dta_col));
fcfg.lbl_two = strcat( 'group_', out_dta_col(use_dta_col));
fcfg.grp_nme = 'side_fisher';
fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/tables/demographics';
ejk_fisher_test( fcfg );

%% Make Table
grp_typ = {'L' 'R'};

tst_tbl = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/tables/demographics/side.ttest2/output_table.csv');
fsh_tbl = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/tables/demographics/side_fisher/output_table.csv');

% Setup table
clear tbl_dsg
for iR = 1:numel(col_hld(:,1))
    for iC = 1:numel(grp_typ)
        if strcmpi(col_hld{iR,end},'num')
            tbl_dsg{iR,iC} = [ 'mean/std' ','  '1' ',' grp_typ{iC} ',' col_hld{iR,1} ];
            row_lbl{iR} = col_hld{iR,1};
        elseif strcmpi(col_hld{iR,end},'str')
            out_dta(cellfun(@isnumeric,out_dta(:,strcmpi(out_dta_col,col_hld{iR,1}))),strcmpi(out_dta_col,col_hld{iR,1})) = {''};            
            hld_cat = unique(out_dta(:,strcmpi(out_dta_col,col_hld{iR,1})));
            hld_cat(cellfun(@isempty,hld_cat)) = [];
            hld_cat = strcat(hld_cat,'/');
            hld_cat = cat(2,hld_cat{:});
            hld_cat = hld_cat(1:end-1);
            tbl_dsg{iR,iC} = [ 'count'    ','  '1' ',' grp_typ{iC} ',' col_hld{iR,1} ','  hld_cat];
            row_lbl{iR} = [ col_hld{iR,1} ' (' hld_cat ')'];
        end
    end
    if strcmpi(col_hld{iR,end},'num')
        tbl_dsg{iR,iC+1} = ['copy' ',' '2' ',' col_hld{iR,1} ',' 'report'];
    elseif strcmpi(col_hld{iR,end},'str')
        tbl_dsg{iR,iC+1} = ['copy' ',' '3' ',' col_hld{iR,1} ',' 'report' ];
    end
end

fcfg = [];
fcfg.tbl = tbl_dsg;
fcfg.dta = {[ out_dta_col ; out_dta ] tst_tbl fsh_tbl };
fcfg.grp = grp;
tbl_out = ejk_create_table( fcfg );

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/current_table.csv',[ {'N'} {numel(grp.L)} {numel(grp.R)} {'Test'} ; row_lbl' tbl_out ]);

%% Save out

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/missing_data/demographic_table.csv',[ {'sbj_nme'} out_dta_col ; out_dta_sbj out_dta ]);















