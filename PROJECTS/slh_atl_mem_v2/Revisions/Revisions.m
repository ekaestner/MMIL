out_dir = [ prj_dir '/' 'out' '/' 'Revisions' '/'];
stt_dir = [ out_dir '/' 'statistics' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(stt_dir);

load([ dta_dir '/' 'group_quality_check.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%% Site Group grp
tot_sbj_ind = [ grp.sde_ana.initial_QC.lft ; grp.sde_ana.initial_QC.rgh];
ste_nme = cell(numel(tot_sbj),1);

% Get Site
ste_nme(string_find(tot_sbj(:,1),{'epd\d'}),end)    = {'UCSD'};
ste_nme(string_find(tot_sbj(:,1),{'epd_ucsf\d'}),end) = {'UCSF'};
ste_nme(string_find(tot_sbj(:,1),{'fc\d'}),end)       = {'HC'};
ste_nme(string_find(tot_sbj(:,1),{'EPK\d'}),end)      = {'Emory'};

% Create Group
usd_ind = find(strcmpi(ste_nme,'UCSD'));
usf_ind = find(strcmpi(ste_nme,'UCSF'));
emy_ind = find(strcmpi(ste_nme,'Emory'));

grp.site.initial_QC.ucsd = intersect(usd_ind,tot_sbj_ind);
grp.site.initial_QC.ucsf = intersect(usf_ind,tot_sbj_ind);
grp.site.initial_QC.emory = intersect(emy_ind,tot_sbj_ind);

%% MTS stats grp
mts_col     = strcmpi(tot_col,'mts'); 
mts_sbj_ind = find(strcmpi(tot_dta(:,mts_col),'yes'));
non_sbj_ind = find(strcmpi(tot_dta(:,mts_col),'no'));

grp.srg_sde_ana.mts_QC.lft_atl = intersect(mts_sbj_ind,grp.srg_sde_ana.initial_QC.lft_atl);
grp.srg_sde_ana.mts_QC.lft_slh = intersect(mts_sbj_ind,grp.srg_sde_ana.initial_QC.lft_slh);
grp.srg_sde_ana.mts_QC.rgh_atl = intersect(mts_sbj_ind,grp.srg_sde_ana.initial_QC.rgh_atl);
grp.srg_sde_ana.mts_QC.rgh_slh = intersect(mts_sbj_ind,grp.srg_sde_ana.initial_QC.rgh_slh);

grp.srg_sde_ana.non_QC.lft_atl = intersect(non_sbj_ind,grp.srg_sde_ana.initial_QC.lft_atl);
grp.srg_sde_ana.non_QC.lft_slh = intersect(non_sbj_ind,grp.srg_sde_ana.initial_QC.lft_slh);
grp.srg_sde_ana.non_QC.rgh_atl = intersect(non_sbj_ind,grp.srg_sde_ana.initial_QC.rgh_atl);
grp.srg_sde_ana.non_QC.rgh_slh = intersect(non_sbj_ind,grp.srg_sde_ana.initial_QC.rgh_slh);

%% Fix Emory Date
dte_fld = { 'neu_psy_tst_dte_pst' 'srg_dte' };
tot_sbj_ind = [ grp.sde_ana.initial_QC.lft ; grp.sde_ana.initial_QC.rgh];

for iF = 1:numel(dte_fld)
    dte_col = strcmpi(tot_col,dte_fld{iF});
    for iS = 1:size(tot_sbj_ind,1)
        str_hld = regexpi(tot_dta{tot_sbj_ind(iS),dte_col},'/','split');
        if numel(str_hld)==3
            tot_dta{tot_sbj_ind(iS),dte_col} = [ str_hld{3} '-' str_hld{1} '-' str_hld{2} ];
        elseif numel(str_hld)==2
            tot_dta{tot_sbj_ind(iS),dte_col} = [ str_hld{2} '-' str_hld{1} '-' '1' ];
        end
    end
end

%% R1.5 Age
% age, lm2_chg
tot_sbj_ind = [ grp.sde_ana.initial_QC.lft ; grp.sde_ana.initial_QC.rgh];
age_col = find(strcmpi(tot_col,'age'));
chg_col = find(strcmpi(tot_col,'lm2_chg'));

[ age_cor_all, age_pvl_all ] = corr(cell2mat(tot_dta(tot_sbj_ind,age_col)),cell2mat(tot_dta(tot_sbj_ind,chg_col)),'Type','Pearson')
[ age_cor_lft, age_pvl_lft ] = corr(cell2mat(tot_dta(grp.sde_ana.initial_QC.lft,age_col)),cell2mat(tot_dta(grp.sde_ana.initial_QC.lft,chg_col)),'Type','Pearson')
[ age_cor_rgh, age_pvl_rgh ] = corr(cell2mat(tot_dta(grp.sde_ana.initial_QC.rgh,age_col)),cell2mat(tot_dta(grp.sde_ana.initial_QC.rgh,chg_col)),'Type','Pearson')
[ age_cor_lft_atl, age_pvl_lft_atl ] = corr(cell2mat(tot_dta(grp.srg_sde_ana.initial_QC.lft_atl,age_col)),cell2mat(tot_dta(grp.srg_sde_ana.initial_QC.lft_atl,chg_col)),'Type','Pearson')
[ age_cor_lft_slh, age_pvl_lft_slh ] = corr(cell2mat(tot_dta(grp.srg_sde_ana.initial_QC.lft_slh,age_col)),cell2mat(tot_dta(grp.srg_sde_ana.initial_QC.lft_slh,chg_col)),'Type','Pearson')

cell2csv()

%% R2.2 / R3.2 Site specific changes
% ILF_fiber_FA_Laterality, Unc_fiber_FA_Laterality Hippocampus_subcort_vol_Laterality
neu_var_nme = { 'ILF_fiber_FA_Laterality' 'Unc_fiber_FA_Laterality' 'Hippocampus_subcort_vol_Laterality' 'lm2_chg' };

neu_typ = [ repmat({'Emory'},1,numel(grp.site.initial_QC.emory)) repmat({'UCSD'},1,numel(grp.site.initial_QC.ucsd)) repmat({'UCSF'},1,numel(grp.site.initial_QC.ucsf)) ];

for iN = 1:numel(neu_var_nme)
    neu_col = strcmpi(tot_col,neu_var_nme{iN});
    neu_dta = cell2mat([  tot_dta(grp.site.initial_QC.emory,neu_col)' tot_dta(grp.site.initial_QC.ucsd,neu_col)' tot_dta(grp.site.initial_QC.ucsf,neu_col)' ]);
    fprintf('%s\n',neu_var_nme{iN})
    [p,tbl] = anova1(neu_dta,neu_typ,'off')
end

%% R3.1 Length of time
% srg_dte
% neu_psy_tst_dte_pst
tot_sbj_ind     = [ grp.sde_ana.initial_QC.lft ; grp.sde_ana.initial_QC.rgh];
srg_dte_col     = find(strcmpi(tot_col,'srg_dte'));
pst_cog_dte_col = find(strcmpi(tot_col,'neu_psy_tst_dte_pst'));

dur_hld = nan(numel(tot_sbj_ind),1);
for iS = 1:numel(tot_sbj_ind)
    % Surgery -to- PostNeuroPsych
    if ~isempty(tot_dta{tot_sbj_ind(iS),srg_dte_col}) && ~isempty(tot_dta{tot_sbj_ind(iS),pst_cog_dte_col})
        
        srg_dte = cellfun(@(x) str2num(x),regexp(tot_dta{tot_sbj_ind(iS),srg_dte_col},'-','split'));
        if numel(srg_dte)==1
            srg_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_srg_dte{iS,1},'/','split'));
            srg_dte = (srg_dte(3)*365) + (srg_dte(1)*30) + srg_dte(2);
        elseif numel(srg_dte)==3
            srg_dte = (srg_dte(1)*365) + (srg_dte(2)*30) + srg_dte(3);
        end
        
        cog_dte = cellfun(@(x) str2num(x),regexp(tot_dta{tot_sbj_ind(iS),pst_cog_dte_col},'-','split'));
        if numel(cog_dte)==1;
            cog_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'/','split'));
            cog_dte = (cog_dte(3)*365) + (cog_dte(1)*30) + cog_dte(2);
        elseif numel(cog_dte)==3
            cog_dte = (cog_dte(1)*365) + (cog_dte(2)*30) + cog_dte(3);
        end
        
        dur_hld(iS,1) = roundsd((cog_dte - srg_dte) / 30,3);
    end  
    
end

[ tot_sbj(tot_sbj_ind) num2cell(dur_hld) ]

dur_hld(dur_hld<0) = [];
dur_hld(dur_hld>100) = [];

min(dur_hld)
max(dur_hld)
median(dur_hld)
iqr(dur_hld)
mean(dur_hld)
std(dur_hld)

%% R4.1 MTS status
grp_use.nme     = { 'srg_sde_ana' };
grp_use.sub_nme = { { 'mts_QC' 'non_QC'} }; %

cog_var = { 'lm2_pre' 'lm2_chg' 'lm2_pst'};

cog_ind = find(ismember(tot_col,cog_var));

cor_ind = 1:size(tot_col,2); cor_ind = setxor(cor_ind,cog_ind); 
cor_ind = cor_ind(all(cellfun(@isnumeric,tot_dta(:,cor_ind)))); 
cor_ind = cor_ind(~all(cellfun(@isempty,tot_dta(:,cor_ind)))); 
cor_ind = [ find(strcmpi(tot_col,'lm2_pre')) cor_ind ];

cor_var = tot_col(cor_ind);

for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        
        for iFN = 1:numel(fld_nme)
            grp_ind = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iFN});
            
            fcfg = [];
            
            fcfg.sbj_nme = tot_sbj( grp_ind, 1);
            
            fcfg.dta_two = cell2mat(tot_dta( grp_ind, cog_ind));
            fcfg.lbl_two = strcat('cog_',tot_col(cog_ind));
            
            fcfg.cor_typ = 'spearman';
            
            fcfg.dta_one = cell2mat(tot_dta( grp_ind, cor_ind));
            fcfg.lbl_one = tot_col(cor_ind);
            
            fcfg.pvl_cut = 0.05;
            fcfg.pvl_lib = 0.10;
            
            fcfg.force_plot = 0;
            
            fcfg.out_dir = [ stt_dir '/' 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} '/' fld_nme{iFN}  '/' ];
            
            ejk_cross_cor( fcfg );
        end
        
    end
end
    
% Bivariate Table
for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        col_add = 0;
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        dta_tbl = cell(numel(cor_var), numel(cog_var)*numel(fld_nme)+numel(cog_var)-1 );
        col_tbl = cell(1,numel(cog_var)*numel(fld_nme)+numel(cog_var)-1);
        
        for iC = 1:numel(cog_var)
            for iFN = 1:numel(fld_nme)
                
                fdr_nme = [ 'correlation_' '_' grp_use.nme{iG} '_'  grp_use.sub_nme{iG}{iGS} ];
                              
                rvl_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_rvalues.csv' ] );
                rvl_hld(strcmpi(rvl_hld,'NA')) = {NaN};
                pvl_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_pvalues.csv' ] );
                pvl_hld(strcmpi(pvl_hld,'NA')) = {NaN};
                num_hld = mmil_readtext([ stt_dir '/' fdr_nme '/' fld_nme{iFN}  '/' 'cross_correlation_n.csv' ] );
                
                col_tbl{iFN+col_add} = [ cog_var{iC} '_' fld_nme{iFN}];
                for iV = 1:numel(cor_var)
                    
                    col_ind = strcmpi(rvl_hld(1,:),[ 'cog' '.' mmil_spec_char(cog_var{iC},{'_'},{'.'}) ]);
                    row_ind = strcmpi(rvl_hld(:,1),mmil_spec_char(cor_var{iV},{'_'},{'.'}) );
                    
                    dta_tbl{iV,iFN+col_add} = [ 'p = ' num2str(roundsd(pvl_hld{row_ind,col_ind},2)) '; r(' num2str(num_hld{row_ind,col_ind}-2) ') = ' num2str(roundsd(rvl_hld{row_ind,col_ind},2))];
                    
                end
            end
            col_add = col_add + numel(fld_nme) + 1;
            
        end
        
        cell2csv( [ stt_dir '/' fdr_nme '/' 'correlation_table.csv'], [ {''} col_tbl ; cor_var' dta_tbl ]);
        
    end
end

%% Surgical sites
slh_num = [ grp.srg_sde_ana.initial_QC.lft_slh ; grp.srg_sde_ana.initial_QC.rgh_slh ];
atl_num = [ grp.srg_sde_ana.initial_QC.lft_atl ; grp.srg_sde_ana.initial_QC.rgh_atl ];

fprintf('Emory: Slah - %i | ATL - %i \n', numel(intersect(slh_num,grp.site.initial_QC.emory)), numel(intersect(atl_num,grp.site.initial_QC.emory)) )
fprintf('UCSD: Slah - %i | ATL - %i \n', numel(intersect(slh_num,grp.site.initial_QC.ucsd)), numel(intersect(atl_num,grp.site.initial_QC.ucsd)) )
fprintf('UCSF: Slah - %i | ATL - %i \n', numel(intersect(slh_num,grp.site.initial_QC.ucsf)), numel(intersect(atl_num,grp.site.initial_QC.ucsf)) )

%% RCI output
tot_sbj_ind = [ grp.srg_sde_ana.initial_QC.lft_atl ; grp.srg_sde_ana.initial_QC.lft_slh ; grp.srg_sde_ana.initial_QC.rgh_atl ; grp.srg_sde_ana.initial_QC.rgh_slh ];
neu_typ = [ repmat({'lft_atl'},numel(grp.srg_sde_ana.initial_QC.lft_atl),1) ; ... 
            repmat({'lft_slh'},numel(grp.srg_sde_ana.initial_QC.lft_slh),1) ; ... 
            repmat({'rgh_atl'},numel(grp.srg_sde_ana.initial_QC.rgh_atl),1) ; ... 
            repmat({'rgh_slh'},numel(grp.srg_sde_ana.initial_QC.rgh_slh),1) ];
lm2_ver_col = find(strcmpi(tot_col,'lm2_ver'));
lm2_chg_col = find(strcmpi(tot_col,'lm2_chg'));
        
cell2csv([ '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/01_Lead/slh_atl_mem_v2/out/Revisions' '/' 'rci.csv'], [ tot_sbj(tot_sbj_ind) neu_typ tot_dta(tot_sbj_ind,[ lm2_ver_col lm2_chg_col ]) ])




