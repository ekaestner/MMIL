

fcfg = [];
fcfg.sbj_nme = cln_dta(2:end,1); % N-x-1 column of subject names being put in
fcfg.dta     = cell2mat(cln_dta(2:end,2:end)); % N-x-M matrix of roi values ( rows = subjects, columns = ROIs)
fcfg.dta_lbl = cln_dta(1,2:end); % 1-x-M row of ROI names
fcfg.out_dir = [prj_dir '/' prj_nme '/' 'Data' '/' 'QC_cleaned' '/' 'Cognitive' '/']; % File location to save plots, etc.
fcfg.out_pre_fix = 'Cognitive'; % Prefix for the files to be saved
ejk_qc_roi(fcfg)