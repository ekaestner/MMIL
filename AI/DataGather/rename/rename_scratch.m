rnm_fld = '/home/ekaestne/gitrep/MMIL/AI/DataGather/rename/';
rnm_fle = 'ENIGMAEpilepsy-RenameReport_DATA_2025-11-21_1548.csv';

onb_fld = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding/';

out_dir = '/space/mcdonald-syn01/1/projects/ekaestner/BIDS_dev/transfer_test/';

spc_ses_fle = 'specify_sessions.csv';

%%
cde_bok_site = {1  '001_CNH' ; ...
                2  '002_UniMoRe' ; ...
                3  'UDA_ENIGMA' ; ...
                4  '004_Rush' ; ...
                5  '005_UniCa' ; ...
                6  '006_USZ' ; ...
                7  '007_UNICAMP' ; ...
                8  '008_GWU' ; ...
                9  '009_UCT' ; ...
                10 '010_Nanjing' ; ...
                11 '011_UCSD' ; ...
                12 '012_ULB' ; ...
                13 '013_EPICZ' ; ...
                14 '014_FSM' ; ...
                15 '015_Gottingen' ; ...
                16 '016_UKBonn' ; ...
                17 '017_MNI' ; ...
                18 '018_NYU' ; ...
                19 '019_IGG' ; ...
                20 '020_UNAM' ; ...
                21 '021_UCLA' ; ...
                22 '022_UPenn' ; ...
                23 '023_MUSC' ; ...
                24 '024_LNF' ; ...
                25 '025_Emory' ; ...
                26 '026_MeyerClinic' ; ...
                27 'Tbingen_ENIGMA' ; ...
                28 '028_UCL' ; ...
                29 '029_MonashatAlfredInstitute' ; ...
                30 '030_UniMe' ; ...
                31 '031_UCSF' ; ...
                32 'Frankfurt_ENIGMA' ; ...
                33 'WULondonOntario_ENIGMA' ; ...
                34 'UCalgary_ENIGMA' ; ...
                35 '035_TLVMC' ; ...
                36 'BWH_ENIGMA' ; ...
                37 '037_HenryFord' ; ...
                99 'other_ENIGMA' };


%%
fcfg = [];
fcfg.dta_loc = [rnm_fld '/' spc_ses_fle];
[ spc_ses_dta, spc_ses_sbj, spc_ses_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [rnm_fld '/' rnm_fle];
[ key_dta, key_sbj, key_col ] = ejk_dta_frm(fcfg);
key_dta = [ key_sbj key_dta  ];
key_col = [ {'site'}  key_col];

new_sbj_col = strcmpi(key_col, 'study_id');
new_ses_col = strcmpi(key_col, 'redcap_repeat_instance');
old_sbj_col = strcmpi(key_col, 'subjid');
old_ses_col = strcmpi(key_col, 'sesid'); 
ste_col     = strcmpi(key_col, 'site');
img_col     = strcmpi(key_col, 'redcap_repeat_instrument');

for iS = 1:numel(key_sbj)
    if ~isempty(key_dta(iS,ste_col)) && ~((iS+1)>numel(key_sbj)) && strcmpi(key_dta(iS+1,img_col),'Imaging')
        key_dta{iS+1,ste_col} = key_dta{iS,ste_col};
    end
end

rmv_row = ~strcmpi(key_dta(:,img_col),'imaging');
key_dta(rmv_row,:) = [];

for iS = 1:size(key_dta,1)
    key_dta{iS,ste_col} = cde_bok_site{cell2mat(cde_bok_site(:,1))==key_dta{iS,ste_col},2};
end

for iS = 1:size(key_dta,1)
    ses_num_str = num2str(key_dta{iS,new_ses_col});
    if numel(ses_num_str)==1
        key_dta{iS,new_ses_col} = [ 'ses-0' ses_num_str];
    elseif numel(ses_num_str)==2
        key_dta{iS,new_ses_col} = [ 'ses-' ses_num_str];
    else
        error()
    end
end

[ numel(unique(key_dta(:,old_sbj_col))) numel(key_dta(:,old_sbj_col)) ]

key_dta(:,old_sbj_col) = strrep(key_dta(:,old_sbj_col),'SMS','sub-SMS');

% Save out subject list
% Get unique subjects and their corresponding site (first occurrence)
[unique_subs, ia] = unique(key_dta(:,2), 'first');
key_dta_site_only = key_dta(ia, [1 2 5]);
cell2csv( '/space/mcdonald-syn01/1/BIDS/ENIGMA_dataset/derivatives/rename_logs/20251124_182819/explore_missing_site_identifiers.csv', ...
    key_dta_site_only);

%%
pst_chk = repmat({'05_unknown_problem'},size(key_dta,1),1);
mss_ste = cell(0);
mss_nme = cell(0);
dbl_nme = cell(0);
ste_dir = dir(onb_fld); ste_dir = { ste_dir(:).name }; 
mask = ~cellfun(@isempty, regexp(ste_dir,'^\d{3}_','once')); % endsWith(ste_dir,'_bids') |
ste_dir = ste_dir(mask);

for iSD = 1:numel(ste_dir)
    sbj_dir = dir([ onb_fld '/' ste_dir{iSD} ]); sbj_dir = { sbj_dir(:).name }; sbj_dir = sbj_dir(startsWith(sbj_dir,'sub-'));

    for iF = 1:numel(sbj_dir)

        ses_fle_row = strcmpi( key_dta(:,old_sbj_col), sbj_dir{iF}) & ...
                      strcmpi(key_dta(:,ste_col), ste_dir{iSD});

        ses_dir = dir([ onb_fld '/' ste_dir{iSD} '/'  sbj_dir{iF}]); ses_dir = ses_dir([ses_dir(:).isdir]); ses_dir = { ses_dir(:).name }; ses_dir = ses_dir(startsWith(ses_dir,'ses-'));

        if ~isempty(ses_dir)
            if sum(ses_fle_row)==1
                key_dta{ses_fle_row,old_ses_col} = ses_dir{:};
                pst_chk(ses_fle_row) = {'01_A_single_assignment'};
            elseif sum(ses_fle_row)==0
                mss_nme = [ mss_nme ;  sbj_dir(iF)];
                mss_ste = [ mss_ste ;  ste_dir(iSD)];
            elseif sum(ses_fle_row) ~= numel(ses_dir)
                if numel(unique(key_dta(ses_fle_row,new_sbj_col))) > 1
                    dbl_nme = [ dbl_nme ;  key_dta(ses_fle_row,old_sbj_col) key_dta(ses_fle_row,new_sbj_col) key_dta(ses_fle_row,new_ses_col)];
                else
                    pst_chk(ses_fle_row) = {'04_Redcap_BIDS_session_mismatch'};
                end
                
            else
                if any(strcmpi(ses_dir,'ses-01')) || any(strcmpi(ses_dir,'ses-1'))
                    new_num = regexp(key_dta(ses_fle_row,new_ses_col), '\d+', 'match'); new_num = cellfun(@str2double,new_num); % Convert to number
                    ses_fle_ind = find(ses_fle_row);
                    for iN = 1:numel(ses_dir)                        
                        old_num = regexp(ses_dir{iN}, '\d+', 'match'); old_num = str2double(old_num{1}); % Convert to number
                        use_row = new_num==old_num;
                        key_dta{ses_fle_ind(use_row),old_ses_col} = ses_dir{use_row};
                        pst_chk(ses_fle_ind(use_row)) = {'01_B_numerical_rename'};
                    end         
                else
                    mul_ses_str = strjoin(ses_dir, ' | ');
                    ses_sbj_ind = find(strcmpi(spc_ses_sbj,unique(key_dta(ses_fle_row,new_sbj_col))));
                    if ~isempty(ses_sbj_ind)
                        ses_fle_ind = find(ses_fle_row);
                        for iSS = 1:numel(ses_fle_ind)
                            key_dta(ses_fle_ind(iSS),old_ses_col) = ...
                                spc_ses_dta( ses_sbj_ind(strcmpi( spc_ses_dta(ses_sbj_ind,strcmpi(spc_ses_col,'new_session')), ...
                                key_dta(ses_fle_ind(iSS),new_ses_col))), ...
                                strcmpi(spc_ses_col,'old_session'));                            
                        end
                        pst_chk(ses_fle_row) = {'03_A_successfully_user_defined_sessions'};
                    else
                        key_dta(ses_fle_row,old_ses_col) = {mul_ses_str} ;
                        pst_chk(ses_fle_row) = {'03_B_failed_user_defined_sessions'};
                    end
                     
                end
            end
        else
            if exist([ onb_fld '/' ste_dir{iSD} '/'  sbj_dir{iF} '/' 'anat' '/'], 'dir') || exist([ onb_fld '/' ste_dir{iSD} '/'  sbj_dir{iF} '/' 'dwi' '/'], 'dir') || exist([ onb_fld '/' ste_dir{iSD} '/'  sbj_dir{iF} '/' 'func' '/'], 'dir')
                if sum(ses_fle_row)==0
                    mss_nme = [ mss_nme ;  sbj_dir(iF)];
                    mss_ste = [ mss_ste ;  ste_dir(iSD)];
                else
                    key_dta{ses_fle_row,old_ses_col} = '999';
                    pst_chk(ses_fle_row) = {'02_no_session_directory'};
                end
            elseif sum(ses_fle_row)==0
                mss_nme = [ mss_nme ;  sbj_dir(iF)];
                mss_ste = [ mss_ste ;  ste_dir(iSD)];
            else
                pst_chk(ses_fle_row) = {'06_no_session_directories'};
                %error('ses_dir develop')
            end
        end
    end
end

dte_str = regexp(rnm_fle, '\d{4}-\d{2}-\d{2}', 'match', 'once');
col_nme = { 'site' 'new_sub' 'type' 'new_ses' 'old_sub' 'old_ses' 'assignment' };

cell2csv( [ out_dir '/' 'renaming_key' '_' dte_str '.csv'], [ col_nme ; key_dta pst_chk ])
cell2csv( [ out_dir '/' 'in_bids_not_in_redcap' '_' dte_str '.csv'], [mss_ste mss_nme] )
cell2csv( [ out_dir '/' 'bids_redcap_mismatch' '_' dte_str '.csv'], dbl_nme )








