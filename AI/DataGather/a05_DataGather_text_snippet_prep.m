for iDT = 1:numel(dta_typ_nme)

    ejk_chk_dir([ new_dta_dir '/' '05_labels' '/' dta_typ_nme{iDT} '/'])

    %%
    fcfg = [];
    fcfg.dta_loc = [ new_dta_dir '/' '03_covariates' '/' dta_typ_nme{iDT} '/' '00_' dta_typ_nme{iDT} '_' 'covariates' '_' dte_str '.csv'];
    fcfg.dta_col = 2;
    [ cov_dta, cov_sbj, cov_col] = ejk_dta_frm( fcfg );

    %%
    ste_var = 'Site';           ste_var_col = strcmpi(cov_col,ste_var);
    age_var = 'Age';           age_var_col = strcmpi(cov_col,age_var);
    sex_var = 'Sex';           sex_var_col = strcmpi(cov_col,sex_var);
    hnd_var = 'Handedness';    hnd_var_col = strcmpi(cov_col,hnd_var);
    edu_var = 'Education';    edu_var_col = strcmpi(cov_col,edu_var);
    
    con_var = 'Diagnosis';  con_var_col = strcmpi(cov_col,con_var);    
    foc_var = 'Focal';      foc_var_col = strcmpi(cov_col,foc_var);
    loc_var = 'Localization';   loc_var_col = strcmpi(cov_col,loc_var);
    lat_var = 'Lateralization'; lat_var_col = strcmpi(cov_col,lat_var);
    
    mts_var     = 'MTS';        mts_var_col = strcmpi(cov_col,mts_var);
    
    age_ons_var = 'AgeOfOnset';  age_ons_var_col = strcmpi(cov_col,age_ons_var);   
    dur_var = 'Duration';  dur_var_col = strcmpi(cov_col,dur_var);  
    drg_res_var = 'DrugResistance';  drg_res_var_col = strcmpi(cov_col,drg_res_var);
    asm_num_var = 'ASM_number';   asm_num_var_col = strcmpi(cov_col,asm_num_var);
    
    ton_cln_var = 'TonicClonic'; ton_cln_var_col = strcmpi(cov_col,ton_cln_var);

    eng_out_var = 'EngelOutcome'; eng_out_var_col = strcmpi(cov_col,eng_out_var);


    %% 
    txt_snp = cell( size(cov_sbj,1), 2 );
    for iS = 1:numel(cov_sbj)
        if strcmpi(cov_dta{iS,sex_var_col},'Female')
            pro_non = 'She';
            pro_non_two = 'Her';
        elseif strcmpi(cov_dta{iS,sex_var_col},'Male')
            pro_non = 'He';
            pro_non_two = 'His';
        end
        
        txt_snp{iS,1} = cov_sbj{iS};

        txt_snp{iS,2} = [ 'Subject ' cov_sbj{iS} ' from site ' cov_dta{iS,ste_var_col} , ...
            ' is a ' num2str(round(cov_dta{iS,age_var_col})) ' year old, ' cov_dta{iS,hnd_var_col} '-handed ' ...
            cov_dta{iS,sex_var_col} ' with ' num2str(round(cov_dta{iS,edu_var_col})) ' years of education.'];


        txt_snp{iS,2} = [ txt_snp{iS,2} ' ' pro_non ' was diagnosed as ' cov_dta{iS,foc_var_col} ' ' ...
            cov_dta{iS,con_var_col} ', with a focus in the ' cov_dta{iS,loc_var_col} ' region on the brain ' ...
            'lateralized to the ' cov_dta{iS,lat_var_col} '.' ];

        txt_snp{iS,2} = [ txt_snp{iS,2} ' ' pro_non ' has ' ];;
        if strcmpi(cov_dta{iS,mts_var_col},'HS_mri') || strcmpi(cov_dta{iS,mts_var_col},'HS_pathology')
            txt_snp{iS,2} = [ txt_snp{iS,2} 'mesial temporal sclerosis (MTS).' ];
        elseif strcmpi(cov_dta{iS,mts_var_col},'noHS')
            txt_snp{iS,2} = [ txt_snp{iS,2} 'no mesial temporal sclerosis (no MTS).' ];
        else
            txt_snp{iS,2} = [ txt_snp{iS,2} 'an unknown MTS status.' ];
        end

        txt_snp{iS,2} = [ txt_snp{iS,2} ' ' pro_non_two ' epilepsy is ' cov_dta{iS,drg_res_var_col} '.' ...
            ' Seizures began at age ' num2str(round(cov_dta{iS,age_ons_var_col})) ' and has lasted for ' ...
            num2str(round(cov_dta{iS,dur_var_col})) ' years. ' pro_non ' is currently taking ' ... ...
            num2str(round(cov_dta{iS,asm_num_var_col})) ' anti-seizure medications. ' ];

        txt_snp{iS,2} = [ txt_snp{iS,2} ' ' pro_non ' has reported ' ];
        if strcmpi(cov_dta{iS,ton_cln_var_col},'Focal')
            txt_snp{iS,2} = [ txt_snp{iS,2} ' only localized seizures.' ];
        elseif strcmpi(cov_dta{iS,ton_cln_var_col},'GTC')
            txt_snp{iS,2} = [ txt_snp{iS,2} ' generalized tonic clonic seizures.' ];
        else
            txt_snp{iS,2} = [ txt_snp{iS,2} ' an seizure severity.' ];
        end

        txt_snp{iS,2} = [ txt_snp{iS,2} ' ' pro_non ' has undergone surgery with an outcome of Engel ' ...
            cov_dta{iS,eng_out_var_col} '.' ];

    end

end