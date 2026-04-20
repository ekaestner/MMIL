clear; clc

%%
org_dir = '/space/mcdonald-syn01/1/data/orig_enigma/uclh_bids/UCL_ENIGMA_BIDS';
bid_dir = '/space/mcdonald-syn01/1/onboarding/locipipe_onboarding/UCL_bids';

rnm_fld = '/home/ekaestne/gitrep/MMIL/AI/DataGather/rename/';
rnm_fle = 'ENIGMAEpilepsy-RenameReport_DATA_2025-04-30_1118.csv';

%%
fcfg = [];
fcfg.dta_loc = [rnm_fld '/' rnm_fle];
[ key_dta, key_sbj, key_col ] = ejk_dta_frm(fcfg);
key_dta = [ key_sbj key_dta  ];
key_col = [ {'site'}  key_col];

ste_col     = strcmpi(key_col, 'site');
img_col     = strcmpi(key_col, 'redcap_repeat_instrument');

for iS = 1:numel(key_sbj)
    if ~isempty(key_dta(iS,ste_col)) && ~((iS+1)>numel(key_sbj)) && strcmpi(key_dta(iS+1,img_col),'Imaging')
        key_dta{iS+1,ste_col} = key_dta{iS,ste_col};
    end
end

rmv_row = ~strcmpi(key_dta(:,img_col),'imaging');
key_dta(rmv_row,:) = [];

key_dta(cellfun(@isempty, key_dta(:,strcmpi(key_col,'site')))) = {NaN};
key_dta = key_dta(cell2mat(key_dta(:,strcmpi(key_col,'site')))==28,:);

%%
for iS = 2:size(key_dta,1)
    old_nme = [bid_dir '/' key_dta{iS,2}];
    new_nme = [bid_dir '/' 'sub-' key_dta{iS,5}];
    movefile(old_nme, new_nme);
end

