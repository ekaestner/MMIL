clear; clc;

slh_atl_mem_v2_constants

out_dir = [ prj_dir '/' 'out' '/' 'Correlations' '/'];
plt_dir = [ out_dir '/' 'plots' '/'];
ejk_chk_dir(out_dir); ejk_chk_dir(plt_dir);

%% Load Data
load([ dta_dir '/' 'group_quality_check.mat' ])

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' 'Total_Demographic_Clinical_Neurobiolgical_Data_recoded.csv' ];
fcfg.dta_col = 2;
[ tot_dta, tot_sbj, tot_col] = ejk_dta_frm( fcfg );

%% Setup figs
cog_nme = { 'lm2_chg'                         'lm2_pre'};
cog_col = [ find(strcmpi(tot_col,cog_nme{1}))  find(strcmpi(tot_col,cog_nme{2})) ];
    
mes_nme = { 'Hippocampus_subcort_vol_Laterality' 'Unc_fiber_FA_Laterality'         'ILF_fiber_FA_Laterality' };
mes_col = [ find(strcmpi(tot_col,mes_nme{1}))    find(strcmpi(tot_col,mes_nme{2})) find(strcmpi(tot_col,mes_nme{3})) ];

grp_typ     = { 'slah'                                      'atl'                                        'combined'};
grp_sub     = { 'srg_sde_ana'                               'srg_sde_ana'                                'sde_ana' };
grp_sub_sub = { 'initial_QC'                                'initial_QC'                                 'initial_QC' };
grp_nme     = { { 'lft_slh'        'rgh_slh' }              { 'lft_atl' 'rgh_atl' }                      { 'lft' 'rgh' } };
grp_col     = { { rgb('royal purple') rgb('faded purple') } { rgb('burnt orange' ) rgb('faded orange') } { rgb('light blue' ) rgb('light red')} };

%% Make scatters
for iC = 1:numel(cog_nme)
    for iN = 1:numel(mes_nme)
        for iT = 1:numel( grp_typ )
            
            % Data Gather
            for iG = 1:numel(grp_nme{iT})
                
                use_ind = grp.(grp_sub{iT}).(grp_sub_sub{iT}).(grp_nme{iT}{iG});
                
                ydt{iG} = cell2mat(tot_dta( use_ind, mes_col(iN) ));
                xdt{iG} = cell2mat(tot_dta( use_ind, cog_col(iC) ));
            end
            
            % Data Plot
            fcfg = [];
            
            fcfg.xdt     = xdt;
            fcfg.ydt     = ydt;
            
            fcfg.fce_col = grp_col{iT};
            fcfg.edg_col = { rgb('black') rgb('black') };
            
            fcfg.ylb = { mes_nme{iN}  };
            fcfg.xlb = { cog_nme{iC} };
            
            fcfg.trd_lne = ones(1,numel(grp_nme{iT}));
            
            fcfg.out_dir = plt_dir;
            fcfg.out_nme = [ grp_typ{iT} '_' mes_nme{iN} '_' cog_nme{iC} ];
            
            ejk_scatter(fcfg)
            
            % Data Plot
            fcfg = [];
            
            fcfg.xdt     = xdt(1);
            fcfg.ydt     = ydt(1);
            
            fcfg.fce_col = grp_col{iT}(1);
            fcfg.edg_col = { rgb('black') };
            
            fcfg.ylb = { mes_nme{iN}  };
            fcfg.xlb = { cog_nme{iC} };
            
            fcfg.trd_lne = [1];
            
            fcfg.out_dir = [ plt_dir '/'  cog_nme{iC} '/' ];
            fcfg.out_nme = [ grp_typ{iT} '_' mes_nme{iN} '_' cog_nme{iC} '_' 'left_only'];
            
            ejk_scatter(fcfg)
           
            clear ydt xdt
            
        end
    end
end