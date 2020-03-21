
function ejk_roi_plot(cfg)

ejk_chk_dir(cfg.out_dir);

%%
sph = { 'lh' 'rh' };
sph_vew = { 'lat' 'ven' 'med' };

%%
[~,albl,~]=fs_read_annotation([ cfg.prj_dir '/' 'DATA' '/'  cfg.fsr_nme '/' 'ROIs' '/' sph{1} '.aparc' cfg.prc_nme '.annot']);
pct_hld = [albl num2cell(zeros(size(albl,1),1))];

ind = 1:size(albl,1);
pct_hld = pct_hld(ind,:);
for iR = ind
    pct_hld(ind(iR),2)     = { iR };
end

if ~strcmpi( cfg.fsr_nme , 'fsaverage' )
    sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.fsr_dir));
    sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.fsr_nme '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
else
    sbj_dir_lst = dir(sprintf('%s/*',cfg.fsr_dir));
    sbj_dir_lst = regexp({sbj_dir_lst.name},['' cfg.fsr_nme],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
end

srf_brn{1} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

num_col = 200;
dst_col = distinguishable_colors(num_col);

col{1} = rgb('grey');
for iR = 1:size(pct_hld,1)
    col{iR + 1} = dst_col(iR,:);
    col{iR + 1} = col{iR + 1} + [ 0.00 0.05 0.125 ];
    if col{iR + 1}(1) > 1; col{iR + 1}(1) = 1; end
    if col{iR + 1}(2) > 1; col{iR + 1}(2) = 1; end
    if col{iR + 1}(3) > 1; col{iR + 1}(3) = 1; end
end

col_map = [];
top_pct = 1;
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
end
top_pct = numel(col)-1;

for iH = 1:2
    for iSP = 1:numel(sph_vew)
        
        fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
        
        pcfg = [];
        
        pcfg.surf_brain  = srf_brn{iH};
        pcfg.aparc       = [ cfg.prj_dir '/' 'DATA' '/'  cfg.fsr_nme '/' 'ROIs' '/' sph{iH} '.aparc' cfg.prc_nme '.annot'];
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iSP};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        if ~strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 1],'visible','off','Parent',fig_hld(1));
        elseif strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 0.6],'visible','off','Parent',fig_hld(1));
        end
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = col_map;
        
        pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
        pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
        
        pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
        print([ cfg.out_dir '/' 'ROI_guide' '_' 'aparc_' mmil_spec_char(cfg.out_nme,{'.'}) '_' sph{iH} '_' sph_vew{iSP} '.png'],'-dpng','-r200')
        close all
        
    end
end

end







