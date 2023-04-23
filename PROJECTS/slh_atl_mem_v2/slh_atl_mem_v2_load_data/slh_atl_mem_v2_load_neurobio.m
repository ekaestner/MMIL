
%%
% Demographics
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Data.csv' ];
fcfg.dta_col = 2;
[ dem_dta, dem_sbj, dem_col] = ejk_dta_frm( fcfg );

%%
inp_fle     = { mri_dev_fle dti_dev_fle };
inp_pre     = { mri_dev_pre dti_dev_pre };
inp_suf     = { mri_dev_suf dti_dev_suf };
inp_mse     = { mri_dev_mse dti_dev_mse };
inp_roi     = { mri_dev_roi dti_dev_roi };
inp_lat     = { mri_dev_lat dti_dev_lat };
inp_icv     = { mri_dev_icv dti_dev_icv };
inp_int     = { mri_var_int dti_var_int };

%%
neu_dta = cell(numel(dem_sbj),100);
neu_col = cell(1,100);
neu_cnt = 1;

%%
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
            
            if inp_roi{iF}{iM}(iR)==0
                roi_fle_nme = [ dta_dir '/' inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '.csv'];
            else
                roi_fle_nme = [ dta_dir '/' 'Data' '/' inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '.csv'];
                qal_dir     = [ dta_dir '/' 'Data' '/' 'QC' '/' inp_pre{iF} '/' inp_mse{iF}{iM}  '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '/' ];
            end
            
            if inp_roi{iF}{iM}(iR)==0
                roi_nme_hld = [];
                roi_hms_hld = [];
            else
                roi_nme_hld = { [ roi_loc '/' 'lh.' roi_nme{inp_roi{iF}{iM}(iR)}] [ roi_loc '/' 'rh.' roi_nme{inp_roi{iF}{iM}(iR)}] };
                roi_hms_hld = { 'lh' 'rh' };
            end
            
            % Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fcfg = [];
            fcfg.sbj_nme = dem_sbj;
            fcfg.fle_nme = inp_fle{iF};
            fcfg.mes_nme = inp_mse{iF}{iM};
            fcfg.rcn_nme = rcn_fle;
            fcfg.roi_nme = roi_nme_hld;
            fcfg.roi_hms = roi_hms_hld;
            
            [ neu_bio_dta, neu_bio_wrn ] = ejk_extract_mmps_roi( fcfg );
            neu_bio_dta(1,:) = cellfun(@(x) strrep(x,[inp_mse{iF}{iM} '_'],''),neu_bio_dta(1,:),'uni',0);
            
            % Save
%             cell2csv( roi_fle_nme, neu_bio_dta)
%             if ~isempty(neu_bio_wrn); cell2csv( [roi_fle_nme(1:end-4) '_warnings.txt' ], neu_bio_wrn); end
             for iRO = 1:numel(inp_int{iF}{iM})
                 roi_ind = string_find(neu_bio_dta(1,:),inp_int{iF}{iM}{iRO});
                 neu_dta(:,neu_cnt:neu_cnt+numel(roi_ind)-1) = neu_bio_dta(2:end,roi_ind);
                 neu_col(neu_cnt:neu_cnt+numel(roi_ind)-1) = strcat(neu_bio_dta(1,roi_ind),'_',inp_mse{iF}{iM});
                 neu_cnt = neu_cnt+numel(roi_ind);
             end
            

            % ICV normalize %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmpi(inp_icv{iF}{iM},'')
                
                fcfg = [];
                
                fcfg.dta     = neu_bio_dta;
                fcfg.cor_col = inp_icv{iF}{iM};
                
                neu_bio_dta_icv = ejk_cor_roi( fcfg );
                
                % Save
%                 cell2csv( [roi_fle_nme(1:end-4) '_' 'norm' '_' inp_icv{iF}{iM} '.csv'], neu_bio_dta_icv)
                 for iRO = 1:numel(inp_int{iF}{iM})
                     roi_ind = string_find(neu_bio_dta_icv(1,:),inp_int{iF}{iM}{iRO});
                     neu_dta(:,neu_cnt:neu_cnt+numel(roi_ind)-1) = neu_bio_dta_icv(2:end,roi_ind);
                     neu_col(neu_cnt:neu_cnt+numel(roi_ind)-1) = strcat(neu_bio_dta_icv(1,roi_ind),'_',inp_mse{iF}{iM},'_','ICVcor');
                     neu_cnt = neu_cnt+numel(roi_ind);
                 end


            end
            
            % Laterality Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if inp_lat{iF}(iM)==1
                
                
                [ dta_out , dta_lbl ] = ejk_create_laterality_index( neu_bio_dta(2:end,:) , neu_bio_dta(1,:));
                neu_bio_dta_lat = [ neu_bio_dta(:,1:4) [ dta_lbl ; num2cell(dta_out) ] ];
                
                % Save
%                 cell2csv( [roi_fle_nme(1:end-4) '_' 'LateralityIndex' '.csv'], neu_bio_dta_lat)
                 for iRO = 1:numel(inp_int{iF}{iM})
                     roi_ind = string_find(neu_bio_dta_lat(1,:),inp_int{iF}{iM}{iRO});
                     neu_dta(:,neu_cnt:neu_cnt+numel(roi_ind)-1) = neu_bio_dta_lat(2:end,roi_ind);
                     neu_col(neu_cnt:neu_cnt+numel(roi_ind)-1) = strcat(neu_bio_dta_lat(1,roi_ind),'_',inp_mse{iF}{iM},'_','Laterality');
                     neu_cnt = neu_cnt+numel(roi_ind);
                 end

            end
            
            % Clean up %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear neu_bio_dta neu_bio_dta_icv neu_bio_dta_lat neu_bio_wrn
            
        end
    end
end

neu_dta(:,cellfun(@isempty,neu_col)) = [];
neu_col(cellfun(@isempty,neu_col))   = [];

%% Combine
tot_dta = [ dem_dta neu_dta ];
tot_col = [ dem_col neu_col ];
tot_sbj = dem_sbj;

%% Save
cell2csv([ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data.csv'], [ {'sbj_nme'} tot_col ; tot_sbj tot_dta]);

