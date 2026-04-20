%% Check
% [ {'sbj'} red_cap_col ; red_cap_sbj(1:6) red_cap_dta(1:6,:) ]
% 
% reg_dta_row = [ 3:20 25:41];
% img_dta_row = [ 1:2 21:24 ];
% 
% [ {'sbj'} red_cap_col(reg_dta_row) ; red_cap_sbj(2:2:6) red_cap_dta(2:2:6,reg_dta_row) ]
% 
% red_cap_col(img_dta_row)
% 
% ucl_ind = cellfun(@(x) x==28,red_cap_dta(:,strcmpi(red_cap_col,'site')),'uni',0);
%     ucl_ind(cellfun(@isempty,ucl_ind)) = {logical(0)};
%     ucl_ind = cell2mat(ucl_ind);
% img_ind = strcmpi(red_cap_dta(:,strcmpi(red_cap_col,'redcap_repeat_instrument')),'imaging');
% cov_ind = ~cellfun(@isempty,red_cap_dta(:,strcmpi(red_cap_col,'old_subjid'))) & ucl_ind;
% 
% [ num2cell(1:numel(red_cap_col))' ...
%   red_cap_col' ...
%   num2cell((sum(~cellfun(@isempty,red_cap_dta(img_ind,:))) / numel(red_cap_dta(img_ind,1)))') ...
%   num2cell((sum(~cellfun(@isempty,red_cap_dta(cov_ind,:))) / numel(red_cap_dta(cov_ind,1)))') ]


%% Tidyverse setup 
reg_dta_row = [ 3:20 25:41];
img_dta_row = [ 1:2 21:24 ];

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

%% Calculations
asm_col = strcmpi(red_cap_col,'current_asm');
sze_typ_col = string_find(red_cap_col,'severity_sz_type');
for iS = 1:numel(red_cap_sbj)

    % Current ASM
    if ~isempty(red_cap_dta{iS,asm_col})
        aed_hld = split(red_cap_dta{iS,asm_col},';');
        red_cap_dta{iS,asm_col} = numel(aed_hld);
    else
        red_cap_dta{iS,asm_col} = double([]);
    end

    % Seizure Severity Type
    if ~all(cellfun(@(x) x==0,red_cap_dta(iS,sze_typ_col)))
        gtc_num = find(cell2mat(red_cap_dta(iS,sze_typ_col)));
        if any(gtc_num==4)
            red_cap_dta(iS,sze_typ_col(1)) = {'GTC'};
        else
            red_cap_dta(iS,sze_typ_col(1)) = {'Focal'};
        end
    else
        red_cap_dta(iS,sze_typ_col(1)) = {cell(0)};
    end
end

% Rename GTC
red_cap_col{sze_typ_col(1)} = 'severity_sz_type';
red_cap_col(sze_typ_col(2:end)) = [];
red_cap_dta(:,sze_typ_col(2:end)) = [];



