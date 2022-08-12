% cfg = [];
% 
% cfg.prj_dir = '/home/ekaestne/PROJECTS/';
% cfg.prj_nme = 'CognitivePhenotype';
% cfg.out_dir = '';
% 
% cfg.grp_dir = '/space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/';
% cfg.grp_nme = {'ALL_42controls_BETAS_Aug31+tlrc' ...
%                 'ALL_15NO_IMP_BETAS_Aug31+tlrc' ...
%                 'ALL_16LANG_BETAS_Aug31+tlrc' ...
%                 'ALL_10MEMORY_BETAS_Aug31+tlrc' ...
%                 'ALL_10LMs_BETAS_Aug31+tlrc'};
% cfg.plt_nme = {'Control' ...
%                 'NoImpairment' ...
%                 'LanguageImpairment' ...
%                 'MemoryImpairment' ...
%                 'LanguageAndMemoryImpairment'};
% 
% cfg.fmr_col_map = {'bright blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright red'};
% cfg.fmr_rng_num = [-10 10];
% 
% cfg.sph     = {'lh' 'rh'};
% cfg.sph_vew = {'lat' 'med' 'ven'};
%
%  
% 3dVol2Surf  -out_1D /home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf/AllControl_lhs.1D
% 
% 3dVol2Surf -spec /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/suma_TT_N27_ROI/TT_N27_both.spec -surf_A rh.pial.gii -sv /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/suma_TT_N27_ROI/TT_N27+tlrc -grid_parent /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/'ALL_42controls_BETAS_Aug31+tlrc[7]' -map_func mask -out_1D /home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf/AllControl_rhs.1D

function mmil_anat_surf_plot(cfg)

if ~isfield(cfg,'hme_wrk'); cfg.hme_wrk=0; end
if ~isfield(cfg,'lhs_nme'); cfg.lhs_nme = 'lh'; end
if ~isfield(cfg,'rhs_nme'); cfg.rhs_nme = 'rh'; end

%% Load Data
hms = {'lhs' 'rhs'};

top_pct = 1;

plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};
sph     = {'lh' 'rh'};
sph_vew = {'lat' 'med' 'ven'};

% plt_dta{1} = fs_load_mgh([cfg.dta_dir '/' cfg.srf_typ_nme{1}]);
% plt_dta{2} = fs_load_mgh([cfg.dta_dir '/' cfg.srf_typ_nme{2}]);

cfg.plt_dta{1}(isnan(cfg.plt_dta{1})) = 0;
cfg.plt_dta{2}(isnan(cfg.plt_dta{2})) = 0;

hld_val = [cfg.plt_dta{1} cfg.plt_dta{2}]; if size(hld_val,2)==2; hld_val = [cfg.plt_dta{1} ; cfg.plt_dta{2}]; end
std_val = std(hld_val);
men_val = mean(hld_val);
max_val = max(abs(men_val-std_val*3),abs(men_val+std_val*3));

%% fsaverage load
if isfield(cfg,'hme_wrk') && cfg.hme_wrk == 1
    lhs_surf_brain.surf_brain =  fs_read_surf(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' cfg.lhs_nme '.inflated']);
    lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
    lhs_srf_lbl = fs_read_label(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/label/' '/' cfg.lhs_nme  '.cortex.label']);
    
    rhs_surf_brain.surf_brain =  fs_read_surf(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' cfg.rhs_nme '.pial']);
    rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
    rhs_srf_lbl = fs_read_label(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/label/' '/' cfg.rhs_nme '.cortex.label']);
elseif isfield(cfg,'hme_wrk') && cfg.hme_wrk == 2
    lhs_surf_brain.surf_brain =  fs_read_surf(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' cfg.lhs_nme '.inflated']);
    lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
    lhs_srf_lbl = fs_read_label(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/label/' '/' cfg.lhs_nme  '.cortex.label']);
    
    rhs_surf_brain.surf_brain =  fs_read_surf(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' cfg.rhs_nme '.inflated']);
    rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
    rhs_srf_lbl = fs_read_label(['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/label/' '/' cfg.rhs_nme '.cortex.label']);
elseif isfield(cfg,'fsr_dir') && ~isempty(cfg.fsr_dir)
    lhs_surf_brain.surf_brain =  fs_read_surf([cfg.fsr_dir '/' 'surf' '/' cfg.lhs_nme  '.pial']);
    lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
    lhs_srf_lbl = fs_read_label([cfg.fsr_dir '/' 'label' '/' cfg.lhs_nme  '.cortex.label']);
    
    rhs_surf_brain.surf_brain =  fs_read_surf([cfg.fsr_dir '/' 'surf' '/' cfg.rhs_nme '.pial']);
    rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
    rhs_srf_lbl = fs_read_label([ cfg.fsr_dir '/' 'label' '/' cfg.rhs_nme '.cortex.label']);
else
    lhs_surf_brain.surf_brain =  fs_read_surf(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'surf' '/' cfg.lhs_nme '.pial']);
    lhs_surf_brain.surf_brain.coords = lhs_surf_brain.surf_brain.vertices;
    lhs_srf_lbl = fs_read_label(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'label' '/' cfg.lhs_nme '.cortex.label']);
    
    rhs_surf_brain.surf_brain =  fs_read_surf(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'surf' '/' cfg.rhs_nme '.pial']);
    rhs_surf_brain.surf_brain.coords = rhs_surf_brain.surf_brain.vertices;
    rhs_srf_lbl = fs_read_label(['/home/mmilmcd/data/FSRECONS/fsaverage/' '/' 'label' '/' cfg.rhs_nme '.cortex.label']);
end

brn_srf{1} = lhs_surf_brain;
brn_srf{2} = rhs_surf_brain;

brn_ctx{1} = lhs_srf_lbl;
brn_ctx{2} = rhs_srf_lbl;

%% Make Colormap
cfg.fmr_col_map = cellfun(@rgb ,cfg.fmr_col_map,'uni',0);

fmr_col_map = [];
for iC = 1:numel(cfg.fmr_col_map)-1
    fmr_col_map = [fmr_col_map ; [linspace(cfg.fmr_col_map{iC}(1),cfg.fmr_col_map{iC+1}(1),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))' linspace(cfg.fmr_col_map{iC}(2),cfg.fmr_col_map{iC+1}(2),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))' linspace(cfg.fmr_col_map{iC}(3),cfg.fmr_col_map{iC+1}(3),ceil(1000*top_pct/(numel(cfg.fmr_col_map)-1)))']; ];
end

%% Plot
% Plot
fig_hld(1) = figure('units','pixels','outerposition',[0 0 1800 1200],'Visible','off');

for iH = 1:numel(sph)   
    
    if isfield(cfg,'deg_fre')
        pvl_grp_low = -4.0; % tinv(10^-fdr_hld(iH),cfg.deg_fre);
        pvl_grp_hgh = 4.0; %tinv(1-10^-fdr_hld(iH),cfg.deg_fre);
        
        mem_grp = cfg.plt_dta{iH};
        mem_grp( mem_grp >  pvl_grp_low & mem_grp <  pvl_grp_hgh ) = 0;
        mem_grp = mem_grp ./ max_val;
        mem_grp( mem_grp <= -1) = -0.99;
        mem_grp( mem_grp >= 1)  = 0.99;
        mem_grp = (mem_grp + 1) / 2;
    elseif isfield(cfg,'low_rng_num') && ~isempty(cfg.low_rng_num)
        
        pvl_grp_low = cfg.low_rng_num(1);
        pvl_grp_hgh = cfg.low_rng_num(2);
        
        mem_grp = cfg.plt_dta{iH};
        
        if isnan(pvl_grp_low)
            error('')
        elseif isnan(pvl_grp_hgh)
            mem_grp( mem_grp>0 ) = 0;
            mem_grp( mem_grp >  pvl_grp_low ) = 0;
            mem_grp( mem_grp <= pvl_grp_low ) = mem_grp( mem_grp <= pvl_grp_low ) ./ abs(cfg.hgh_rng_num(1)) ;
            mem_grp( mem_grp <= -1) = -0.99;
            mem_grp = (mem_grp + 1) / 2;
        else
            mem_grp( mem_grp >  pvl_grp_low & mem_grp <  pvl_grp_hgh ) = 0;
            mem_grp( mem_grp <= pvl_grp_low | mem_grp >= pvl_grp_hgh ) = mem_grp(mem_grp <= pvl_grp_low | mem_grp >= pvl_grp_hgh) ./ cfg.hgh_rng_num(2);
            mem_grp( mem_grp <= -1) = -0.99;
            mem_grp( mem_grp >= 1)  = 0.99;
            mem_grp = (mem_grp + 1) / 2;
        end
        
    elseif isfield(cfg,'pvl_rng_num') && ~isempty(cfg.pvl_rng_num)
        mem_grp = cfg.plt_dta{iH};
        mem_grp = (mem_grp-1)*-1;
        mem_grp(mem_grp < 1-cfg.pvl_rng_num(2)) = 0.5;
    else
        mem_grp = repmat(0.5,1,numel(grp{iG}{iH}));
    end
       
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
        
        pcfg.col_map = fmr_col_map;
        
        pcfg.tbl_pct = mem_grp;
        pcfg.tbl_loc = 1:numel(mem_grp);
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
    end
    
end

% Colorbar
if ~isfield(cfg,'fmr_rng_num') || isempty(cfg.fmr_rng_num)
    cfg.fmr_rng_num = [-roundsd(max_val,3) roundsd(max_val,3)];
end

if isfield(cfg,'hme_wrk') && cfg.hme_wrk == 1

end
    
ax1 = axes('OuterPosition',[0.85 0.20 0.04 0.60],'visible','off','Parent',fig_hld(1));

colormap(ax1,fmr_col_map)
clb = colorbar('west','Position',[0.92 0.20 0.02 0.60]);
clb.TickLength = 0;
col_num = [ num2cell(roundsd(linspace(cfg.hgh_rng_num(1),cfg.low_rng_num(1),5),2)) ...
            {0} ...
            num2cell(roundsd(linspace(cfg.low_rng_num(2),cfg.hgh_rng_num(2),5),2)) ];
clb.TickLabels = cellfun(@num2str,col_num,'uni',0);

% Save Figure
print(gcf,[cfg.out_dir '/' cfg.out_pre_fix '.png'],'-dpng','-r200')
close all
    
end
