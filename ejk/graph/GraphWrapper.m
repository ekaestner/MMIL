

function GraphWrapper(cfg)

if ~isfield( cfg , 'ovr_wrt' ); cfg.ovr_wrt = 0; end

ejk_chk_dir( cfg.out_dir )

%% Calculate Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cor_mtx con_smp dff_smp

if ~( exist( [ cfg.out_dir '/' 'GraphThickMatrix.mat' ], 'file') == 2 ) || cfg.ovr_wrt == 1
    
    fcfg = [];
    
    fcfg.sbj_grp = cfg.sbj_nme;
    
    fcfg.sbj_grp_col = cfg.sbj_grp_col;
    fcfg.sbj_grp_nme = cfg.sbj_grp_nme;
    
    fcfg.dta     = cfg.thk_dta;
    fcfg.roi_nme = cfg.thk_roi_nme;
    
    fcfg.clc_con_smp_mtx = 1;
    fcfg.con_smp_nme = cfg.pvl_con;
    fcfg.con_smp_rep = cfg.num_rep;
    fcfg.lve_num_out = cfg.lve_num_out;
    
    fcfg.clc_dff_mtx = 0;
    
    fcfg.out_dir = cfg.out_dir;
    fcfg.plt = 0;
    
    [ cor_mtx , con_smp , ~ ] = GraphThickMatrix(fcfg);
    
    save( [ cfg.out_dir '/' 'GraphThickMatrix.mat' ] , 'cor_mtx' , 'con_smp' )
    
else
    
    load( [ cfg.out_dir '/' 'GraphThickMatrix.mat' ] , 'cor_mtx' , 'con_smp' )
    
end

%% Binarize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.sbj_grp_col = cfg.sbj_grp_col;
fcfg.sbj_grp_nme = cfg.sbj_grp_nme;

fcfg.bin_dir     = cfg.bin_dir;
fcfg.bin_cut_off = cfg.bin_cut_off;

fcfg.cor_mtx         = cor_mtx;

fcfg.clc_con_smp_mtx = 1;
fcfg.con_smp         = con_smp;

fcfg.clc_dff_mtx = 0;

[ cor_mtx_bin , con_smp_bin , ~ ] = GraphBinzarize(fcfg);

%% Calculate and Plot Graph Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iGL = 1:numel(cfg.gbl_mes_nme)
    
    for iC = 1:numel(cfg.sbj_grp_col); out_dir{iC} = [ cfg.out_dir '/' cfg.sbj_grp_col{iC}]; ejk_chk_dir(out_dir{iC}); sve_nme{iC} = [ cfg.gbl_mes_plt{iGL} '_' cfg.gbl_mes_nme{iGL} '_' cfg.sbj_grp_col{iC} ]; end
        
    if ~( exist( [ cfg.out_dir '/' 'GraphData' '_' cfg.gbl_mes_nme{iGL} '.mat' ], 'file') == 2 ) || cfg.ovr_wrt == 1
        
        % Calculate
        fcfg = [];
        
        fcfg.mes = cfg.gbl_mes_nme{iGL};
        
        fcfg.dta = cor_mtx_bin;
        
        fcfg.con_spr_stt = 1;
        fcfg.con_spr_dta = con_smp_bin;
        
        fcfg.dff_cmp_stt = 0;
        
        plt_dta = GraphCalc(fcfg);
        
        save( [ cfg.out_dir '/' 'GraphData' '_' cfg.gbl_mes_nme{iGL} '.mat' ] , 'plt_dta' , '-v7.3');
        
    else
        
        load( [ cfg.out_dir '/' 'GraphData' '_' cfg.gbl_mes_nme{iGL} '.mat' ] );
        
    end
        
    % Plot    
    fcfg = [];
    
    fcfg.grp_ovr = cfg.sbj_grp_col;
    fcfg.grp_nme = cfg.sbj_grp_nme;
    fcfg.grp_col = cfg.sbj_grp_clr;
    
    fcfg.leg_loc = 'southwest';
    
    fcfg.bin_cut = cfg.bin_cut_off;
    
    fcfg.grp_dta = plt_dta;
    fcfg.thk_nme = cfg.thk_roi_nme;
    
    fcfg.con_typ = 'con_shf';
    fcfg.con_nme = cfg.pvl_con;
    fcfg.cmp_nme = cfg.pvl_cmp;
    
    fcfg.con_spr_plt = 1;
    
    fcfg.con_cmp_stt = 0;
    
    fcfg.pvl_lvl = cfg.pvl_lvl;
    
    fcfg.grp_plt     = cfg.grp_plt(iGL);
    fcfg.reg_lne_plt = cfg.reg_plt(iGL);
    
    fcfg.reg_srf_plt = cfg.reg_plt(iGL);
    fcfg.reg_col = cfg.reg_col;
    fcfg.fsr_dir = '/home/mmilmcdRSI/data/fsurf';
    fcfg.fsr_nme = 'fsaverage';
    fcfg.prj_dir = cfg.prj_dir;
    fcfg.sbj_nme = 'fsaverage';
    fcfg.prc_nme = cfg.prc_nme;
    
    fcfg.ttl     = cfg.gbl_mes_nme{iGL};
    fcfg.sve_loc = out_dir;
    fcfg.sve_nme = sve_nme;
    
    fcfg.hme_wrk = 1;
    
    GraphLinePlot(fcfg)
    
    %
    clear plt_dta
    
end

end