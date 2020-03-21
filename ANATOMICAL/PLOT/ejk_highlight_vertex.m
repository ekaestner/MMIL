function ejk_highlight_vertex(cfg)

top_pct = 1;

plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};
sph     = {'lh' 'rh'};
sph_vew = {'lat' 'med' 'ven'};

%% fsaverage load
lhs_surf_brain.surf_brain =  fs_read_surf(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'surf' '/' 'lh.inflated']);
lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
lhs_srf_lbl = fs_read_label(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'label' '/' 'lh.cortex.label']);

rhs_surf_brain.surf_brain =  fs_read_surf(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'surf' '/' 'rh.inflated']);
rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
rhs_srf_lbl = fs_read_label(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'label' '/' 'rh.cortex.label']);

brn_srf{1} = lhs_surf_brain;
brn_srf{2} = rhs_surf_brain;

brn_ctx{1} = lhs_srf_lbl;
brn_ctx{2} = rhs_srf_lbl;

%% Make Colormap
for iH = 1:numel(sph)
    cfg.vtx_col{iH} = [rgb('dark grey') cfg.vtx_col{iH}];
    fmr_col_map{iH} = [];
    for iC = 1:numel(cfg.vtx_col{iH})-1
        fmr_col_map{iH} = [fmr_col_map{iH} ; [linspace(cfg.vtx_col{iH}{iC}(1),cfg.vtx_col{iH}{iC+1}(1),ceil(1000*top_pct/(numel(cfg.vtx_col{iH})-1)))' linspace(cfg.vtx_col{iH}{iC}(2),cfg.vtx_col{iH}{iC+1}(2),ceil(1000*top_pct/(numel(cfg.vtx_col{iH})-1)))' linspace(cfg.vtx_col{iH}{iC}(3),cfg.vtx_col{iH}{iC+1}(3),ceil(1000*top_pct/(numel(cfg.vtx_col{iH})-1)))']; ];;
    end
    fmr_col_map{iH} = [ fmr_col_map{iH} ; cfg.vtx_col{iH}{end} ];
end

%% Make Vertices
for iH = 1:numel(sph)
    
    mem_grp{iH} = zeros( 1 , brn_srf{1}.surf_brain.nverts );
        
    for iV = 1:numel(cfg.vtx_cor{iH})
        
        pst_fun = @(x) pdist( [ x ; brn_srf{1}.surf_brain.coords(cfg.vtx_cor{iH}(iV),:) ] );
        app_to2_gvn_row = @(func,matrix) @(row) func(matrix(row,:));
        app_to2_row     = @(func,matrix) arrayfun(app_to2_gvn_row(func,matrix),1:size(matrix,1),'Uni',0)';
        tke_all = @(x) reshape([x{:}],size(x{1},2),size(x,1))';
        app_row = @(func,matrix) tke_all(app_to2_row(func,matrix));
        
        dst_hld = app_row(pst_fun,brn_srf{1}.surf_brain.coords);
        [ ~ , dst_hld_srt] = sort(dst_hld);
        
        mem_grp{iH}(dst_hld_srt(1:40)) = iV;
        
    end
    
    mem_grp{iH} = (mem_grp{iH} / iV) + .001;
    
end

%% Plot
fig_hld(1) = figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');

for iH = 1:numel(sph)
    for iV = 1:numel(sph_vew)
               
        pcfg = [];
        
        pcfg.surf_brain  = brn_srf{iH};
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iV};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iV},'visible','off','Parent',fig_hld(1));
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = fmr_col_map{iH};
        
        pcfg.tbl_pct = mem_grp{iH};
        pcfg.tbl_loc = 1:numel(mem_grp{iH});
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
    end
end

%% Save Figure
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '.png'],'-dpng','-r200')
close all


end