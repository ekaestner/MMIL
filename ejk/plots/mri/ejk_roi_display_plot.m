function ejk_roi_display_plot(cfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cfg = [];
% 
% cfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
% cfg.fsr_nme = 'fsaverage';
% 
% cfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
% cfg.prc_nme = '.aparc.annot';
% 
% cfg.col_map = { rgb('purplish blue') rgb('blue') rgb('grayish blue') rgb('light grey') rgb('reddish grey') rgb('red') rgb('orangish red') };
% cfg.low_rng_num = [ -0.20 0.20 ];
% cfg.hgh_rng_num = [ -0.50 0.50 ];
% 
% cfg.roi_tbl_loc = [ {plt_tbl{1}(:,1)}           {plt_tbl{2}(:,1)}    ];
% cfg.roi_tbl = [ {cell2mat(plt_tbl{1}(:,2))} {cell2mat(plt_tbl{2}(:,2))}    ];
% 
% cfg.sph = { 'lh' 'rh' };
% cfg.sph_vew = { 'lat' 'ven' 'med' };
% 
% cfg.out_dir = '/home/ekaestner/Downloads/';
% cfg.out_nme  = 'test_roi';
% 
% ejk_roi_display_plot(fcfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 2/3/2021 by ejk 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
if ~isfield(cfg,'sph');     cfg.sph = { 'lh' 'rh' }; end
if ~isfield(cfg,'sph_vew'); cfg.sph_vew = { 'lat' 'ven' 'med' }; end
if numel(cfg.low_rng_num)==1; if cfg.low_rng_num>=0; cfg.low_rng_num = [ 0 cfg.low_rng_num]; elseif cfg.low_rng_num<0; cfg.low_rng_num = [ cfg.low_rng_num 0]; end; end
if numel(cfg.hgh_rng_num)==1; if cfg.hgh_rng_num>=0; cfg.hgh_rng_num = [ 0 cfg.hgh_rng_num]; elseif cfg.hgh_rng_num<0; cfg.hgh_rng_num = [ cfg.hgh_rng_num 0]; end; end
    
%% Load Pieces
% Get fsurf direction
if ~strcmpi( cfg.fsr_nme , 'fsaverage' )
    sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.fsr_dir));
    sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.fsr_nme '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
else
    sbj_dir_lst = dir(sprintf('%s/*',cfg.fsr_dir));
    sbj_dir_lst = regexp({sbj_dir_lst.name},['' cfg.fsr_nme],'match'); sbj_dir_lst = [sbj_dir_lst{:}];
end

% .pial files
srf_brn{1} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.coords = srf_brn{1}.vertices;
srf_brn{1}.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.coords = srf_brn{2}.vertices;
srf_brn{2}.faces = srf_brn{2}.faces;

%% Put together vertex color
col_map = [];
for iC = 1:numel(cfg.col_map)-1
    col_map = [col_map ; [linspace(cfg.col_map{iC}(1),cfg.col_map{iC+1}(1),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(2),cfg.col_map{iC+1}(2),ceil(1000/(numel(cfg.col_map)-1)))' linspace(cfg.col_map{iC}(3),cfg.col_map{iC+1}(3),ceil(1000/(numel(cfg.col_map)-1)))']; ];
end

for iH = 1:numel(srf_brn)
    
    [col_loc,albl,actbl] = fs_read_annotation([ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ]);
    
    % color vertices
    if isfield(cfg,'inc_reg'); [ ~ , lbl_ind ] = intersect(ejk_switch_names(albl,cfg.clr_fld,'split'),cfg.inc_reg); else lbl_ind = 1:numel(albl); end
    
    use_plt = cfg.roi_tbl{iH};
    cst_adj = ( 1000 - 1 ) / ( cfg.hgh_rng_num(2) - cfg.hgh_rng_num(1) );
    for iC = 1:size(use_plt,1)
        if use_plt(iC) > cfg.low_rng_num(1) && use_plt(iC) < cfg.low_rng_num(2); use_plt(iC) = mean(cfg.hgh_rng_num); end
        if use_plt(iC) < cfg.hgh_rng_num(1); use_plt(iC) = cfg.hgh_rng_num(1); end
        if use_plt(iC) > cfg.hgh_rng_num(2); use_plt(iC) = cfg.hgh_rng_num(2); end
        use_plt(iC) = round(cst_adj * ( use_plt(iC) - cfg.hgh_rng_num(1)) + 1);    
    end
    
    for iLC = 1:numel(albl)        
        tbl_ind = find(strcmpi(cellfun(@(x) x(4:end),cfg.roi_tbl_loc{iH},'uni',0),albl{iLC}));
        if isfield(cfg,'inc_reg'); if ~any(tbl_ind==lbl_ind); tbl_ind = []; end; end
        if ~isempty(tbl_ind) && sum(col_loc==iLC)>0
            col_hld(iLC,:) = col_map(use_plt(tbl_ind),:);
        else
            col_hld(iLC,:) = rgb('light grey');
        end
    end
    
    col{iH} = repmat(rgb('white'), [size(srf_brn{iH}.coords, 1) 1]);
    
    if size(albl,1) == size(col_loc,1)
        error('Fix Me')
    else
        for iCL = 1:size(col_hld,1)
            if any(iCL==col_loc)
                col{iH}(iCL==col_loc,:) = repmat(col_hld(iCL,:),sum(iCL==col_loc),1);
            end
        end
    end
    
end

%% Plot
fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);

plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.0 0.4 0.2] [0.0 0.2 0.4 0.4]} {[0.4 0.6 0.4 0.4] [0.4 0.0 0.4 0.2] [0.4 0.2 0.4 0.4]}};

for iH = 1:numel(cfg.sph)
    for iSP = 1:numel(cfg.sph_vew)       
        
        pcfg = [];
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.surf_brain  = srf_brn{iH};
        pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
        
        pcfg.sph         = cfg.sph{iH};
        pcfg.sph_vew     = cfg.sph_vew{iSP};
        
        pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
                
        pcfg.col_map = col{iH};
                        
        pcfg.top_pct = 1;
        
        if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
        
        pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
        
        pcfg.sve_img = 1;
        
        ejk_surf_plot(pcfg);
                
    end
end

%% Colorbar
ax1 = axes('OuterPosition',[0.85 0.20 0.04 0.60],'visible','off','Parent',fig_hld(1));

colormap(ax1,col_map)
clb = colorbar('west','Position',[0.92 0.20 0.02 0.60]);
clb.TickLength = 0;
col_num = [ num2cell(roundsd(linspace(cfg.hgh_rng_num(1),cfg.low_rng_num(1),5),2)) ...
            {0} ...
            num2cell(roundsd(linspace(cfg.low_rng_num(2),cfg.hgh_rng_num(2),5),2)) ];
clb.TickLabels = cellfun(@num2str,col_num,'uni',0);


%% Save
ejk_chk_dir(cfg.out_dir)
print([ cfg.out_dir '/' cfg.out_nme  '.png'],'-dpng','-r200')
close all

end







