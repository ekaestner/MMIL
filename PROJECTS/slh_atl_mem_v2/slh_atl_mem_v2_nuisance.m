out_dir = [ prj_dir '/' 'out' '/' 'Nuisance' '/'];
cog_out_dir = [ out_dir '/' 'Cognitive' '/'];
cor_out_dir = [ out_dir '/' 'Correlation' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(cog_out_dir); ejk_chk_dir(cor_out_dir);

%%
cog_var = { 'lm2_pre' 'lm2_chg'};
cor_var = { 'fusiform_gwcsurf_FA_wm_ctx_Laterality' 'lateralorbitofrontal_gwcsurf_FA_wm_ctx_Laterality' 'ILF_fiber_FA_Laterality' 'Unc_fiber_FA_Laterality' ...
            'Hippocampus_subcort_vol_Laterality' 'Thalamus_subcort_vol_Laterality' };

grp_use.nme     = { 'srg_sde'              'sde' };
grp_use.sub_nme = { {'non_qc' 'non_qc_QC'} {'non_qc' 'non_qc_QC'} };

clear rcd_nme rcd_val
 
rcd_nme     = { 'eng_out' 'eng_out' 'lng_lat'};
rcd_val{1}  = { 'I' 'I' ; 'II' 'II+' ; 'III' 'II+' ; 'IV' 'II+'  };
rcd_val{2}  = { 'L' 'L' ; 'R' 'Atypical' ; 'B' 'Atypical' };

eng_val     = { 'I' 'II+' };
eng_val_col = { rgb('slate blue') rgb('brick red') };

lng_val     = { 'L' 'Atypical' };
lng_val_col = { rgb('forrest green') rgb('light brown') };

%% Load Data
load([ dta_dir '/' 'group_quality_check.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%% Check existence
[ ~, col_ind] = intersect(tot_col,{ 'eng_out' 'lng_lat' });
[ tot_sbj([ grp.sde.non_qc_QC.lft ; grp.sde.non_qc_QC.rgh ]) tot_dta( [ grp.sde.non_qc_QC.lft ; grp.sde.non_qc_QC.rgh ], col_ind) ]

%% Recode
% fcfg = [];
% 
% fcfg.dta     = tot_dta;
% fcfg.dta_col = tot_col;
% 
% fcfg.rcd_nme = rcd_nme;
% fcfg.rcd_val = rcd_val;
% 
% fcfg.swt_val = 1;
% 
% fcfg.frc_fll = { 'mts' };
% 
% tot_dta_rcd = ejk_recode(fcfg);

%% Cognitive Swarm Plots
% Engel
for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        
        for iC = 1:numel(cog_var)
            
            fcfg = [];
            
            cnt = 1;
            for iF = 1:numel(fld_nme)
                pos_num = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iF});
                for iE = 1:numel(eng_val)
                    
                    fcfg.ydt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'eng_out')),eng_val{iE}) ), strcmpi(tot_col,cog_var{iC})) );
                    fcfg.xdt{cnt}     = (iF-0.2) + (0.4*(iE-1));
                    fcfg.fce_col{cnt} = eng_val_col{iE};
                    fcfg.box_plt_col{cnt} = eng_val_col{iE};
                    fcfg.xlb{cnt} = [ fld_nme{iF} ' ' eng_val{iE}];
                    cnt = cnt + 1;
                end
            end
            
            fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
            
            fcfg.box_plt = ones(1,numel(fcfg.xdt));
            
            fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
            fcfg.ylb = cog_var(iC);
            
            fcfg.jtr = 1;
            fcfg.jtr_wdt = 0.15;
            fcfg.box_wdt = 0.15;
            
            fcfg.out_dir = [ cog_out_dir '/' grp_use.sub_nme{iG}{iGS} '/' ];
            fcfg.out_nme = [ 'Engel' '_' grp_use.nme{iG} '_'  cog_var{iC} ];
            
            ejk_scatter(fcfg)
            
        end
    end
end

% Lang Lat
for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        
        for iC = 1:numel(cog_var)
            
            fcfg = [];
            
            cnt = 1;
            for iF = 1:numel(fld_nme)
                pos_num = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iF});
                for iE = 1:numel(lng_val)
                    if ~isempty(pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'lng_lat')),lng_val{iE}) ))
                        fcfg.ydt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'lng_lat')),lng_val{iE}) ), strcmpi(tot_col,cog_var{iC})) );
                        fcfg.xdt{cnt}     = (iF-0.2) + (0.4*(iE-1));
                        fcfg.fce_col{cnt} = lng_val_col{iE};
                        fcfg.box_plt_col{cnt} = lng_val_col{iE};
                        fcfg.xlb{cnt} = [ fld_nme{iF} ' ' lng_val{iE}];
                        cnt = cnt + 1;
                    end
                end
            end
            
            fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
            
            fcfg.box_plt = ones(1,numel(fcfg.xdt));
            
            fcfg.xlm = [ 0.5 max([fcfg.xdt{:}])+0.5 ];
            fcfg.ylb = cog_var(iC);
            
            fcfg.jtr = 1;
            fcfg.jtr_wdt = 0.15;
            fcfg.box_wdt = 0.15;
            
            fcfg.out_dir = [ cog_out_dir '/' grp_use.sub_nme{iG}{iGS} '/' ];
            fcfg.out_nme = [ 'LanguageLaterality' '_' grp_use.nme{iG} '_'  cog_var{iC} ];
            
            ejk_scatter(fcfg)
            
        end
    end
end

%% Correlation Scatterplots
% Engel
for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        
        for iC = 1:numel(cog_var)
            for iCR = 1:numel(cor_var)               
                
                for iF = 1:numel(fld_nme)
                    
                    pos_num = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iF});
                    
                    fcfg = [];
                    
                    fcfg.xdt{1}     = cell2mat(tot_dta_rcd( pos_num, strcmpi(tot_col,cog_var{iC})) );
                    fcfg.ydt{1}     = cell2mat(tot_dta_rcd( pos_num, strcmpi(tot_col,cor_var{iCR})) );
                    fcfg.fce_col{1} = rgb('black');
                    
                    cnt = 2;                    
                    for iE = 1:numel(eng_val)
                        
                        fcfg.xdt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'eng_out')),eng_val{iE}) ), strcmpi(tot_col,cog_var{iC})) );
                        fcfg.ydt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'eng_out')),eng_val{iE}) ), strcmpi(tot_col,cor_var{iCR})) );
                        fcfg.fce_col{cnt} = eng_val_col{iE};
 
                        cnt = cnt + 1;
                    end
                    
                    fcfg.trd_lne = ones(1,numel(fcfg.xdt));
                    
                    fcfg.edg_col = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
                    
                    fcfg.xlb = cor_var(iCR);
                    fcfg.ylb = cog_var(iC);
                    
                    fcfg.ttl = fld_nme{iF};
                    
                    fcfg.out_dir = [ cor_out_dir '/' grp_use.sub_nme{iG}{iGS} '/' fld_nme{iF} '/'];
                    fcfg.out_nme = [ 'Engel' '_' grp_use.nme{iG} '_' fld_nme{iF} '_' cog_var{iC} '_' cor_var{iCR} ];
                    
                    ejk_scatter(fcfg)
                    
                end
            end
        end
    end
end

% Language Laterality
for iG = 1:numel(grp_use.nme)
    for iGS = 1:numel(grp_use.sub_nme{iG})
        
        fld_nme = fieldnames(grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}));
        
        for iC = 1:numel(cog_var)
            for iCR = 1:numel(cor_var)
                
                for iF = 1:numel(fld_nme)
                    
                    pos_num = grp.(grp_use.nme{iG}).(grp_use.sub_nme{iG}{iGS}).(fld_nme{iF});
                    
                    if ~isempty(pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'lng_lat')),lng_val{iE}) ))
                        fcfg = [];
                        
                        
                        fcfg.xdt{1}     = cell2mat(tot_dta_rcd( pos_num, strcmpi(tot_col,cog_var{iC})) );
                        fcfg.ydt{1}     = cell2mat(tot_dta_rcd( pos_num, strcmpi(tot_col,cor_var{iCR})) );
                        fcfg.fce_col{1} = rgb('black');
                        
                        cnt = 2;
                        for iE = 1:numel(lng_val)
                            
                            fcfg.xdt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'lng_lat')),lng_val{iE}) ), strcmpi(tot_col,cog_var{iC})) );
                            fcfg.ydt{cnt}     = cell2mat(tot_dta_rcd( pos_num( strcmpi(tot_dta_rcd(pos_num,strcmpi(tot_col,'lng_lat')),lng_val{iE}) ), strcmpi(tot_col,cor_var{iCR})) );
                            fcfg.fce_col{cnt} = lng_val_col{iE};
                            
                            cnt = cnt + 1;
                        end
                        
                        fcfg.trd_lne = ones(1,numel(fcfg.xdt));
                        
                        fcfg.edg_col = [ repmat({[0 0 0]},1,numel(fcfg.ydt)) ];
                        
                        fcfg.xlb = cor_var(iCR);
                        fcfg.ylb = cog_var(iC);
                        
                        fcfg.ttl = fld_nme{iF};
                        
                        fcfg.out_dir = [ cor_out_dir '/' grp_use.sub_nme{iG}{iGS} '/' fld_nme{iF} '/'];
                        fcfg.out_nme = [ 'LanguageLaterality' '_' grp_use.nme{iG} '_' fld_nme{iF} '_' cog_var{iC} '_' cor_var{iCR} ];
                        
                        ejk_scatter(fcfg)
                    end
                end
            end
        end
    end
end















