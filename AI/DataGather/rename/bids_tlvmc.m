clear; clc

%%
onb_dir = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding/035_TLVMC_bids';
bid_dir = '/space/mcdonald-syn01/1/BIDS/ENIGMA_dataset/';

rnm_fld = '/home/ekaestne/gitrep/MMIL/AI/DataGather/rename/';
rnm_fle = 'TLVMC_rename.csv';

%%
fcfg = [];
fcfg.dta_loc = [rnm_fld '/' rnm_fle];
[ key_dta, key_sbj, key_col ] = ejk_dta_frm(fcfg);
key_dta = [ key_sbj key_dta  ];
key_col = [ {'Site'}  key_col];

bid_sbj_col     = strcmpi(key_col, 'BIDS_ID');
bid_ses_col     = strcmpi(key_col, 'BIDS_session');
ste_sbj_col     = strcmpi(key_col, 'Site_ID');
ste_ses_col     = strcmpi(key_col, 'Site_session');

%%
out_log = cell(size(key_dta,1),1);
for iS = 1:size(key_dta,1)
    
    ejk_chk_dir([onb_dir '/' key_dta{iS,ste_sbj_col}]);
    
    old_nme = [bid_dir '/' key_dta{iS,bid_sbj_col} '/' key_dta{iS,bid_ses_col} ];
    new_nme = [onb_dir '/' key_dta{iS,ste_sbj_col} '/' key_dta{iS,ste_ses_col} ];
    cmd = sprintf('cp -aL %s %s', old_nme, new_nme);

    if ~exist(old_nme,"dir")
        fprintf('No folder: %s \n',old_nme)
        out_log{iS} = sprintf('No folder: %s',old_nme);
    else
        out_log{iS} = sprintf('%s',cmd);
        status = system(cmd);
        fprintf('Completed: %s \n',new_nme)
    end
   
end

