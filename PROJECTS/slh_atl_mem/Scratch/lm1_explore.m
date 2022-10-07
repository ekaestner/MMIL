
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

out_dir = [ prj_dir '/' prj_nme '/' 'FinalAnalysis' '/'];

%% Load & LM1
% Load Clinical data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
fcfg.dta_col = 2;
[ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );

% Load Redcap
fcfg = [];
fcfg.sbj_nme = cln_dta_sbj;
fcfg.red_fle = red_cap_fle;
fcfg.sep     = '|';
[sbj_dem , sbj_sze , sbj_scn , sbj_cog, sbj_emo, sbj_srg] = mmil_load_redcap(fcfg);

% Load Emory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Emory' '/' 'Emory_SLAH_Memory_11.22.21_cleanNONAMES.csv'];
fcfg.dta_col = 2;
[ emy_dta, emy_dta_sbj, emy_dta_col] = ejk_dta_frm( fcfg );

%% Get Data
nme_add = { 'sbj_cog' };
hld_add.sbj_cog = sbj_cog;

usd_emy_cph{1}  = { 'log_mem_raw_scr_one' 'LM I raw_pre' ; ...
                    'log_mem_nor_scr_one' 'LM I norm_ss_pre' ; ...
                    'log_mem_raw_scr_one_pst' 'LM I raw_6m.post' ; ...
                    'log_mem_nor_scr_one_pst' 'LM I norm_ss_6m.post' };

for iFN = 1:numel(nme_add)
    sbj_num = numel(hld_add.(nme_add{iFN}).sbj_nme);
    for iS = 1:size(emy_dta_sbj,1)
        hld_add.(nme_add{iFN}).sbj_nme{sbj_num+iS,1} = emy_dta_sbj{iS,1};
        for iAD = 1:size(usd_emy_cph{iFN},1)
           if isempty(usd_emy_cph{iFN}{iAD,2})
               hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
           elseif isempty(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})})
               hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
           elseif strcmpi(emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})},'n/a')
               hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = NaN;
           else
               hld_add.(nme_add{iFN}).( usd_emy_cph{iFN}{iAD,1} )(sbj_num+iS) = emy_dta{iS,strcmpi(emy_dta_col,usd_emy_cph{iFN}{iAD,2})};
           end
        end        
    end
end

kep_fld_nme = fieldnames([ hld_add.sbj_cog ]);
hld_add.sbj_cog = rmfield(hld_add.sbj_cog,kep_fld_nme(~ismember(kep_fld_nme,[ 'sbj_nme' ; usd_emy_cph{1}(:,1)])));

sbj_cog = hld_add.sbj_cog;

fcfg = [];
pst_cog_dta = ejk_post_cognitive(fcfg,sbj_cog);

cog_out_col = { 'lm1_pre' 'lm1_pst' 'lm1_chg' 'lm1_rci' 'lm1_pct' };
cog_out_dta = [ sbj_cog.log_mem_nor_scr_one sbj_cog.log_mem_nor_scr_one_pst pst_cog_dta.raw.log_mem_nor_scr_one_pst pst_cog_dta.rci.log_mem_nor_scr_one_pst pst_cog_dta.pct.log_mem_nor_scr_one_pst ];

% Save out
cell2csv( [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cog_LM1.csv'] , [ [{'sbj_nme'} cog_out_col] ; [pst_cog_dta.sbj_nme num2cell(cog_out_dta)] ])

%% Cognitive Analysis
% build table
tbl_typ = { 'mean/std' 'lm1_pre' '' 1; ...
            'mean/std' 'lm1_chg' '' 1 ; ...
            'mean/std' 'lm1_rci' '' 1};

% Load Final data
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv'];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] = ejk_dta_frm( fcfg );

fnl_dta     = [ fnl_dta     num2cell(cog_out_dta(:,ismember(cog_out_col,tbl_typ(:,2))))];
fnl_dta_col = [ fnl_dta_col cog_out_col(ismember(cog_out_col,tbl_typ(:,2)))];

dta_inp{1} = [ fnl_dta_col ; fnl_dta ];

% ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_nme     = { 'surgery' };
grp_nme_sub = { {'pst_cog' 'pst_cog_dti'} };

cog_var = tbl_typ(:,2);
col_ind = find(ismember(fnl_dta_col,cog_var));

for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        % get numeric group data
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));

        % ANOVA
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        fcfg.grp_inc = {fld_nme};
        fcfg.grp_nme = {fld_nme};
        fcfg.dta = fnl_dta(:,col_ind(anv_var));
        fcfg.sbj = fnl_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_LM1' '/'  ];
        ejk_1way_anova( fcfg )
                
    end
end

% 2-way ANOVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        % get numeric group data
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));

        % ANOVA
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        fcfg.grp_inc = {fld_nme};
        fcfg.grp_nme = {fld_nme};
        fcfg.dta = fnl_dta(:,col_ind(anv_var));
        fcfg.sbj = fnl_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
        
        grp_typ_hld = regexp(grp_typ{1},'_','split');
        grp_typ{1}  = cat(1,grp_typ_hld{:});
        
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = fnl_dta_col(col_ind(anv_var));
        fcfg.grp_one     = grp_typ{1}(:,1);
        fcfg.grp_nme_one = {'side'};
        fcfg.grp_two     = grp_typ{1}(:,2);
        fcfg.grp_nme_two = grp_nme(iG);
        fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_LM1' '/'  ];
        ejk_2way_anova( fcfg )
                
    end
end

% 1-sample ttest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        % get numeric group data
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));

        % ANOVA
        anv_var = cellfun(@isnumeric,fnl_dta(1,col_ind));
        
        for iGN = 1:numel(fld_nme)
            fcfg = [];
            fcfg.grp     = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
            fcfg.grp_inc = {fld_nme(iGN)};
            fcfg.grp_nme = {fld_nme(iGN)};
            fcfg.dta = fnl_dta(:,col_ind(anv_var));
            fcfg.sbj = fnl_dta_sbj;
            [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
            
            % Test
            fcfg = [];
            fcfg.sbj_nme = grp_sbj{1};
            fcfg.dta     = grp_dta{1};
            fcfg.dta_nme = fnl_dta_col(:,col_ind(anv_var));
            fcfg.men     = 0;
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'Cognitive_ttest1' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_LM1' '/' fld_nme{iGN} '/' ];
            ejk_ttest1( fcfg );
        end
        
    end
end

% Build Cognitive table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iG = 1:numel(grp_nme)
    for iGS = 1:numel(grp_nme_sub{iG})
        
        dta_inp{2} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_1wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_LM1' '/' grp_nme{iG} '/' 'output_table.csv']);
        dta_inp{3} = mmil_readtext([ out_dir '/' 'stats' '/' 'Cognitive_2wayANOVA' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_LM1' '/' 'side' '/' 'output_table.csv']);
                
        fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
        
        % Create Table
        fcfg = [];
        
        for iR = 1:size(tbl_typ)
            for iN = 1:numel(fld_nme)
                fcfg.tbl{iR,iN} = [ tbl_typ{iR,1} ',' num2str(tbl_typ{iR,4}) ',' fld_nme{iN} ',' tbl_typ{iR,2}];
            end
            fcfg.tbl{iR,iN+1} = ['copy,2,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report'];
            fcfg.tbl{iR,iN+2} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_one'];
            fcfg.tbl{iR,iN+3} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_two'];
            fcfg.tbl{iR,iN+4} = ['copy,3,' mmil_spec_char(tbl_typ{iR,2},{'_'},{'.'}) ',report_int'];
        end
        
        fcfg.dta = dta_inp;
        fcfg.grp = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS});
        tbl_out  = ejk_create_table( fcfg );
        
        % Add-ons
        num_sbj{1} = 'N';
        tbl_lbl{1} = '';
        for iN = 1:numel(fld_nme)
            num_sbj{iN+1} = numel(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iN}));
            tbl_lbl{iN+1} = fld_nme{iN};
        end
        num_sbj{iN+2} = '-';num_sbj{iN+3} = '-';num_sbj{iN+4} = '-';num_sbj{iN+5} = '-';
        tbl_lbl{iN+2} = 'Stats';
        tbl_lbl{iN+3} = 'Side';
        tbl_lbl{iN+4} = 'Surgery';
        tbl_lbl{iN+5} = 'Side*Surgery';
        
        % Out
        tbl_out = [ tbl_lbl; num_sbj; tbl_typ(:,2) tbl_out ];
        
        cell2csv([ out_dir '/' 'tables' '/' 'Table2'  '/' 'Cognitive' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '_LM1.csv'],tbl_out); clear fcfg tbl_out num_sbj tbl_lbl
        
    end
end

%% Correlation analysis
cor_typ = { 'spearman' }; % 'pearson' };

grp_nme     = { 'surgery' };
grp_nme_sub = { {'pst_cog' 'pst_cog_dti'} };

cog_var = tbl_typ(:,2);
cog_ind = find(ismember(fnl_dta_col,cog_var));

cor_ind = 1:size(fnl_dta_col,2); cor_ind = setxor(cor_ind,cog_ind); cor_ind = cor_ind(cellfun(@isnumeric,fnl_dta(1,cor_ind))); cor_ind = [ find(strcmpi(fnl_dta_col,'lm2_pre')) cor_ind ];
cor_var = fnl_dta_col(cor_ind);

for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            
            fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
            
            for iFN = 1:numel(fld_nme)
                grp_ind = grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}).(fld_nme{iFN});
                
                fcfg = [];
                
                fcfg.sbj_nme = fnl_dta_sbj( grp_ind, 1);
                
                fcfg.dta_two = cell2mat(fnl_dta( grp_ind, cog_ind));
                fcfg.lbl_two = strcat('cog_',fnl_dta_col(cog_ind));
                
                fcfg.cor_typ = cor_typ{iCT};
                
                fcfg.dta_one = cell2mat(fnl_dta( grp_ind, cor_ind));
                fcfg.lbl_one = fnl_dta_col(cor_ind);
                
                fcfg.pvl_cut = 0.05;
                fcfg.pvl_lib = 0.10;
                
                fcfg.force_plot = 0;
                
                fcfg.out_dir = [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' fld_nme{iFN}  '/' ];
                
                ejk_cross_cor( fcfg );
            end
            
        end
    end
end

% Fisher's Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grp_lod = { 'ltle_atl' 'ltle_slah' };
for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            
            for iGL = 1:numel(grp_lod)
                grp_rvl{iGL} = mmil_readtext( [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' grp_lod{iGL}  '/' 'cross_correlation_rvalues.csv'] );
                grp_num{iGL} = mmil_readtext( [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' grp_lod{iGL}  '/' 'cross_correlation_n.csv'] );
            end
            
            fcfg = [];
            
            fcfg.col_lbl = grp_rvl{1}(1,2:end);
            fcfg.row_lbl = grp_rvl{1}(2:end,1);
            
            fcfg.rvl_one = grp_rvl{1}(2:end,2:end);
            fcfg.rvl_two = grp_rvl{2}(2:end,2:end);
            
            fcfg.num_one = grp_num{1}(2:end,2:end);
            fcfg.num_two = grp_num{2}(2:end,2:end);
            
            fcfg.out_dir = [ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' 'FishersZ'  '/' ];
            fcfg.out_pre = 'ltle_surgery';
            
            ejk_fishersZ(fcfg);
            
        end
    end
end

% Make Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iCT = 1:numel(cor_typ)
    for iG = 1:numel(grp_nme)
        for iGS = 1:numel(grp_nme_sub{iG})
            for iCV = 1:numel(cog_var)
                
                fld_nme = fieldnames(grp.(grp_nme{iG}).(grp_nme_sub{iG}{iGS}));
                
                ref_tbl = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' fld_nme{1}  '/' 'cross_correlation_rvalues.csv' ]);
                
                bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+2 );
                bil_tbl_hld(:,1) = ref_tbl(:,1);
                
                for iFN = 1:numel(fld_nme)
                    
                    bil_tbl_hld{1,iFN+1} = fld_nme{iFN};
                    
                    cor_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' fld_nme{iFN}  '/' 'cross_correlation_rvalues.csv' ]);
                    cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                    pvl_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' fld_nme{iFN}  '/' 'cross_correlation_pvalues.csv' ]);
                    pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                    num_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' fld_nme{iFN}  '/' 'cross_correlation_n.csv' ]);
                    num_hld = num_hld(2:end,strcmpi(num_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                    num_hld(cellfun(@isstr,num_hld)) = {NaN};
                                                            
                    bil_tbl_hld(2:end,iFN+1) = cellfun(@(x,y,z) [ 'p=' num2str(roundsd(y,2)) ' ; r(' num2str(z) ')=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,num_hld,'uni',0);
                end
                
                fsh_hld = mmil_readtext([ out_dir '/' 'stats' '/' 'correlation_' '_' grp_nme{iG} '_'  grp_nme_sub{iG}{iGS} '_' cor_typ{iCT} '_LM1' '/' 'FishersZ'  '/' 'ltle_surgery_fishersZ.csv' ]);
                bil_tbl_hld(:,iFN+2) = fsh_hld(:,strcmpi(fsh_hld(1,:),strcat('cog.',mmil_spec_char(cog_var{iCV},{'_'},{'.'}))));
                
                cell2csv([ out_dir '/' 'tables' '/' 'Table3'  '/' 'Correlations' '_' grp_nme{iG} '_' grp_nme_sub{iG}{iGS} '_' cog_var{iCV} '_' cor_typ{iCT} '_LM1' '.csv'],bil_tbl_hld);
                
            end
        end
    end
end

%% Cognitive Examine 
% 2x2 ANOVA


% L-TLE ttest




