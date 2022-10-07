clear; clc;

slh_atl_mem_constants

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Manuscript/figures/Figure2';

%% Load Data
load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

% CHECK FOR NAN/NAN/0 RCI problem
fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'FinalSample.csv' ];
fcfg.dta_col = 2;
[ fnl_dta, fnl_dta_sbj, fnl_dta_col] =  ejk_dta_frm( fcfg );

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'additional_variables.csv'];
fcfg.dta_col = 2;
[ cln_nui_dta, cln_nui_dta_sbj, cln_nui_dta_col] = ejk_dta_frm( fcfg );

% Combine
fnl_dta     = [ fnl_dta     cln_nui_dta];
fnl_dta_col = [ fnl_dta_col cln_nui_dta_col ];

fnl_dta(cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_pre'))) & cellfun(@isnan,fnl_dta(:,strcmpi(fnl_dta_col,'lm2_chg'))),strcmpi(fnl_dta_col,'lm2_rci')) = {NaN};

%% Setup figs
cog_nme = { 'lm2_rci' };
cog_col = [ find(strcmpi(fnl_dta_col,cog_nme{1})) ];
    
mes_nme = { 'Unc'                                 'fusiform'                            'lateralorbitofrontal' };
mes_col = [ find(strcmpi(fnl_dta_col,mes_nme{1})) find(strcmpi(fnl_dta_col,mes_nme{2})) find(strcmpi(fnl_dta_col,mes_nme{3})) ];

grp_typ = { 'slah'                                 'atl' };
grp_nme = { { 'ltle_slah' 'rtle_slah' }            { 'ltle_atl' 'rtle_atl' } };
grp_col = { { rgb('royal purple') rgb('faded purple') } { rgb('burnt orange' ) rgb('faded orange') } };

%% Make scatters
for iC = 1:numel(cog_nme)
    for iN = 1:numel(mes_nme)
        for iT = 1:numel( grp_typ )
            
            % Data Gather
            for iG = 1:numel(grp_nme{iT})
                use_ind = grp.surgery.pst_cog_dti.(grp_nme{iT}{iG});
                ydt{iG} = cell2mat(fnl_dta( use_ind, mes_col(iN) ));
                xdt{iG} = cell2mat(fnl_dta( use_ind, cog_col(iC) ));
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
            
            fcfg.out_dir = out_dir;
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
            
            fcfg.out_dir = out_dir;
            fcfg.out_nme = [ grp_typ{iT} '_' mes_nme{iN} '_' cog_nme{iC} '_' 'left_only'];
            
            ejk_scatter(fcfg)
           
            clear ydt xdt
            
        end
    end
end
