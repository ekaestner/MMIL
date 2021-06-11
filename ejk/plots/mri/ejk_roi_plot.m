
function ejk_roi_plot(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fcfg = [];
% 
% fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
% fcfg.fsr_nme = 'fsaverage';
% 
% fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD';
% fcfg.prc_nme = '.aparc.annot';
% 
% fcfg.inc_reg = { 'Entorhinal' 'Parahippocampal' };
% 
% fcfg.sph = { 'lh' 'rh' };
% fcfg.sph_vew = { 'lat' 'ven' 'med' };
% 
% fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena';
% fcfg.out_nme  = 'roi_for_alena';
% 
% ejk_roi_plot(fcfg);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 2/3/2021 by ejk 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
if ~isfield(cfg,'sph');     cfg.sph = { 'lh' 'rh' }; end
if ~isfield(cfg,'sph_vew'); cfg.sph_vew = { 'lat' 'ven' 'med' }; end
if ~isfield(cfg,'sph_vew'); cfg.sph_vew = { 'lat' 'ven' 'med' }; end

ejk_chk_dir(cfg.out_dir)

%% Load Pieces
% load annotations
[~,albl,~] = fs_read_annotation([ cfg.roi_loc '/' cfg.sph{1} cfg.prc_nme ]);
pct_hld    = [albl num2cell(zeros(size(albl,1),1))];

ind = 1:size(albl,1);
pct_hld = pct_hld(ind,:);
for iR = ind; pct_hld(ind(iR),2)     = { iR }; end

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
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf([ cfg.fsr_dir '/' sbj_dir_lst{:} '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

% Set up colormap
num_col = numel(albl);
dst_col = distinguishable_colors(num_col);

if isfield(cfg,'roi_col')
    for iRO = 1:numel(cfg.inc_reg)
        roi_ind = find(strcmpi(albl,cfg.inc_reg{iRO}));
        dst_col(roi_ind,:) = cfg.roi_col{iRO};
    end
end

col{1} = rgb('light grey');
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

%% Put together vertex color
tbl_pct                      = cell2mat(pct_hld(:,2)) ./ top_pct;
tbl_pct(isnan(tbl_pct)) = 0;

num_col = numel(albl);
dst_col = distinguishable_colors(num_col);

if isfield(cfg,'roi_col')
    for iRO = 1:numel(cfg.inc_reg)
        roi_ind = find(strcmpi(albl,cfg.inc_reg{iRO}));
        dst_col(roi_ind,:) = cfg.roi_col{iRO};
    end
end

col{1} = rgb('light grey');
for iR = 1:size(pct_hld,1)
    col{iR + 1} = dst_col(iR,:);
    col{iR + 1} = col{iR + 1} + [ 0.00 0.05 0.125 ];
    if col{iR + 1}(1) > 1; col{iR + 1}(1) = 1; end
    if col{iR + 1}(2) > 1; col{iR + 1}(2) = 1; end
    if col{iR + 1}(3) > 1; col{iR + 1}(3) = 1; end
end

col_map = [];
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000/(numel(col)-1)))']; ];
end

for iH = 1:numel(srf_brn)
    
    tbl_loc = strcat(cfg.sph{iH},'_',pct_hld(:,1));
    
    [col_loc,albl,actbl] = fs_read_annotation([ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ]);
    if isfield(cfg,'inc_reg'); [ ~ , lbl_ind ] = intersect(albl,cfg.inc_reg); else lbl_ind = 1:numel(albl); end
     
    % color vertices     
    for iLC = 1:numel(albl)        
        tbl_ind = find(strcmpi(cellfun(@(x) x(4:end),tbl_loc,'uni',0),albl{iLC}));
        if isfield(cfg,'inc_reg'); if ~any(tbl_ind==lbl_ind); tbl_ind = []; end; end
        if ~isempty(tbl_ind) && sum(col_loc==iLC)>0
            col_hld(iLC,:) = col_map(round(tbl_pct(tbl_ind)*1000),:);
        else
            col_hld(iLC,:) = rgb('light grey');
        end
    end
    
    col{iH} = repmat(rgb('white'), [size(srf_brn{iH}.nverts, 1) 1]);
    
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
for iH = 1:numel(cfg.sph)
    for iSP = 1:numel(cfg.sph_vew)
        
        fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
        
        pcfg = [];
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.surf_brain  = srf_brn{iH};
        pcfg.aparc       = [ cfg.roi_loc '/' cfg.sph{iH} cfg.prc_nme ];
        
        pcfg.sph         = cfg.sph{iH};
        pcfg.sph_vew     = cfg.sph_vew{iSP};
        
        if ~strcmpi(cfg.sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 1],'visible','off','Parent',fig_hld(1));
        elseif strcmpi(cfg.sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 0.6],'visible','off','Parent',fig_hld(1));
        end
                
        pcfg.col_map = col{iH};
               
        pcfg.tbl_loc = strcat(cfg.sph{iH},'_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        if isfield(cfg,'inc_reg'); pcfg.inc_reg = cfg.inc_reg; end
        
        pcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';
        
        pcfg.sve_img = 1;
        
        ejk_surf_plot(pcfg);
        
        print([ cfg.out_dir '/' 'ROI_guide' '_' 'aparc_' mmil_spec_char(cfg.out_nme,{'.'}) '_' cfg.sph{iH} '_' cfg.sph_vew{iSP} '.png'],'-dpng','-r200')
        close all
        
    end
end

end







