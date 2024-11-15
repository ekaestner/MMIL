%% Tidyverse setup 
reg_dta_row = [ 3:15];
img_dta_row = [ 1:2 16:19 ];

cmb_sbj = unique(red_cap_sbj);
cmb_dta = cell( numel(cmb_sbj), numel(red_cap_col) );
for iS = 1:numel(cmb_sbj)

    sbj_row = find(strcmpi(red_cap_sbj,cmb_sbj{iS}));
    if numel(sbj_row)>2
        sbj_row( find(cell2mat(cellfun(@(x) x>1,red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instance')),'uni',0 )))+1 ) = [];
    end
    reg_row = sbj_row( ~strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );
    img_row = sbj_row( strcmpi(red_cap_dta(sbj_row,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging') );

    cmb_dta(iS,reg_dta_row) = red_cap_dta(reg_row,reg_dta_row);
    if ~isempty(img_row); cmb_dta(iS,img_dta_row) = red_cap_dta(img_row,img_dta_row); end

    % fix ses BS
    ses_pos = string_find( cmb_dta(iS,strcmpi(red_cap_col,'subjid')), 'ses' );
    if ~isempty(ses_pos)
        ses_pos = strfind( cmb_dta{iS,strcmpi(red_cap_col,'subjid')}, 'ses' );
        cmb_dta{iS,strcmpi(red_cap_col,'subjid')}(ses_pos-1:end) = [];
    end

end

red_cap_dta = cmb_dta;
red_cap_sbj = cmb_sbj;







