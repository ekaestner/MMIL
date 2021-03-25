%% DTI
dti_dta_all = mmil_readtext('/home/ekaestne/PROJECTS/DATA/ROIHOLD_old/DTI_all_mmilmcdRSI.csv');

dti_col_nme = dti_dta_all(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), dti_col_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = dti_col_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = dti_col_nme{iC};
    end
end
col_typ = unique(col_typ);

%%
rsi_dta_all = mmil_readtext('/home/ekaestne/PROJECTS/DATA/ROIHOLD/RSI_res_ind_free_mmilmcdRSI.csv');

rsi_col_nme = rsi_dta_all(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), rsi_col_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = rsi_col_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = rsi_col_nme{iC};
    end
end
col_typ = unique(col_typ);

%%















