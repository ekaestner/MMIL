for iDT = 1:numel(dta_typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '04_labels' '/' typ_nme{iDT} '/'])

    %%
    fcfg = [];
    fcfg.dta_loc = [ new_dta_dir '/' '03_Datasets' '/' '00_' typ_nme{iDT} '_' 'dataset' '_' dte_str '.csv'];
    fcfg.dta_col = 2;
    [ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );


    %%  SETUP VARIABLES
    ses_var     = 'ses_nme';     var_col.ses_var_col = strcmpi(tot_col,ses_var);
    dta_set_var = 'dta_set_str'; var_col.dta_set_var_col = strcmpi(tot_col,dta_set_var);
    ste_var = 'Site';            var_col.ste_var_col = strcmpi(tot_col,ste_var);
    con_var = 'Diagnosis';       var_col.con_var_col = strcmpi(tot_col,con_var);
    loc_var = 'Localization';    var_col.loc_var_col = strcmpi(tot_col,loc_var);
    lat_var = 'Lateralization';  var_col.lat_var_col = strcmpi(tot_col,lat_var);
    dta_set_var = 'dataset';     var_col.dta_set_col = strcmpi(tot_col,dta_set_var);
    ste_var     = 'Site';        var_col.ste_var_col = strcmpi(tot_col,ste_var);
    mts_var = 'MTS';             var_col.mts_var_col = strcmpi(tot_col,mts_var);
    foc_var = 'Focal';           var_col.foc_var_col = strcmpi(tot_col,foc_var);

    %% SETUP CATEGORIES
    var_col_lop = fieldnames(var_col);
    for iV = 1:numel(var_col_lop)
        var_cat.(var_col_lop{iV}) = unique(tot_dta(:,var_col.(var_col_lop{iV})));
    end

    for iV = 1:numel(var_col_lop)
        for iC = 1:numel(var_cat.(var_col_lop{iV}))

            if strcmpi(var_cat.(var_col_lop{iV}){iC},'')
                grp.(var_col_lop{iV}).EMPTY = ...
                    find(strcmpi(tot_dta(:,var_col.(var_col_lop{iV})),var_cat.(var_col_lop{iV}){iC}));
            elseif ~isempty(strfind(var_cat.(var_col_lop{iV}){iC},'-'))
                grp.(var_col_lop{iV}).(strrep(var_cat.(var_col_lop{iV}){iC},'-','_')) = ...
                    find(strcmpi(tot_dta(:,var_col.(var_col_lop{iV})),var_cat.(var_col_lop{iV}){iC}));
            else
                grp.(var_col_lop{iV}).(var_cat.(var_col_lop{iV}){iC}) = ...
                    find(strcmpi(tot_dta(:,var_col.(var_col_lop{iV})),var_cat.(var_col_lop{iV}){iC}));
            end         

        end
    end

    %% Save
    save([ new_dta_dir '/' '04_labels' '/' typ_nme{iDT} '/' typ_nme{iDT} '_' 'group.mat'], 'grp')

end
