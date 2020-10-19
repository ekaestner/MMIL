clear; clc;

%% Colors
mci_clr = [ 163 86 148 ] / 255;
tle_clr = [ 86 152 163 ] / 255;
con_clr = [ 0 33 87 ] / 255;

mci_hgh_clr = [ 255 0   205 ] / 255;
tle_hgh_clr = [ 0   218 255 ] / 255;
con_hgh_clr = [ 229 239 255 ] / 255;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision';

%% Figure 1C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.sbj_grp_col     = { 'Diagnosis'                                    };
cfg.sbj_grp_nme     = { { 'HC'        'MCI'       'EPD_Old' }         };
cfg.sbj_grp_clr     = { { con_clr     mci_clr     tle_clr  }     };
cfg.grp_col_hgh_lgh = { { con_hgh_clr mci_hgh_clr tle_hgh_clr }     };

%% Top - Path Efficiency
load( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'TestResiduals' '/' 'Thickness_Desikan_Residuals' '/' 'GraphData' '_' 'PathEfficiency' '.mat' ] );

plt_dta_trm = plt_dta;
plt_dta_trm.dta     = plt_dta_trm.dta(2);
plt_dta_trm.con_spr = plt_dta_trm.con_spr(2);

cfg.bin_cut_off = 10:5:50;

%
fcfg = [];

fcfg.grp_ovr         = cfg.sbj_grp_col;
fcfg.grp_nme         = cfg.sbj_grp_nme;
fcfg.grp_col         = cfg.sbj_grp_clr;
fcfg.grp_col_hgh_lgh = cfg.grp_col_hgh_lgh;

fcfg.bin_cut = cfg.bin_cut_off;

fcfg.grp_dta = plt_dta_trm;

fcfg.con_typ = 'con_shf';
fcfg.con_nme = { 'HC' };
fcfg.cmp_nme = { { 'MCI' 'EPD_Old' } };

fcfg.con_spr_plt = 1;

fcfg.pvl_lvl = .001;

fcfg.grp_plt     = 1;
fcfg.reg_lne_plt = 0;
cfg.reg_srf_plt  = 0;

fcfg.ttl     = 'PathEfficiency';
fcfg.sve_loc = {'/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision'};
fcfg.sve_nme = {'PathEfficiency'};

fcfg.ylm     = [.10 .80];

fcfg.sig_loc = 'bottom';

fcfg.hme_wrk = 1;

GraphLinePlot(fcfg)

%% Middle - Transitivity
load( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'TestResiduals' '/' 'Thickness_Desikan_Residuals' '/' 'GraphData' '_' 'Transitivity' '.mat' ] );

plt_dta_trm = plt_dta;
plt_dta_trm.dta     = plt_dta_trm.dta(2);
plt_dta_trm.con_spr = plt_dta_trm.con_spr(2);

cfg.bin_cut_off = 10:5:50;

%
fcfg = [];

fcfg.grp_ovr         = cfg.sbj_grp_col;
fcfg.grp_nme         = cfg.sbj_grp_nme;
fcfg.grp_col         = cfg.sbj_grp_clr;
fcfg.grp_col_hgh_lgh = cfg.grp_col_hgh_lgh;

fcfg.bin_cut = cfg.bin_cut_off;

fcfg.grp_dta = plt_dta_trm;

fcfg.con_typ = 'con_shf';
fcfg.con_nme = { 'HC' };
fcfg.cmp_nme = { { 'MCI' 'EPD_Old' } };

fcfg.con_spr_plt = 1;

fcfg.pvl_lvl = .001;

fcfg.grp_plt     = 1;
fcfg.reg_lne_plt = 0;
cfg.reg_srf_plt  = 0;

fcfg.ttl     = 'Transitivity';
fcfg.sve_loc = {'/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision'};
fcfg.sve_nme = {'Transitivity'};

fcfg.ylm     = [.35 .85];

fcfg.sig_loc = 'bottom';

fcfg.hme_wrk = 1;

GraphLinePlot(fcfg)

%% Bottom - Modularity
load( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'TestResiduals' '/' 'Thickness_Desikan_Residuals' '/' 'GraphData' '_' 'Modularity' '.mat' ] );

plt_dta_trm = plt_dta;
plt_dta_trm.dta     = plt_dta_trm.dta(2);
plt_dta_trm.con_spr = plt_dta_trm.con_spr(2);

cfg.bin_cut_off = 10:5:50;

%
fcfg = [];

fcfg.grp_ovr         = cfg.sbj_grp_col;
fcfg.grp_nme         = cfg.sbj_grp_nme;
fcfg.grp_col         = cfg.sbj_grp_clr;
fcfg.grp_col_hgh_lgh = cfg.grp_col_hgh_lgh;

fcfg.bin_cut = cfg.bin_cut_off;

fcfg.grp_dta = plt_dta_trm;

fcfg.con_typ = 'con_shf';
fcfg.con_nme = { 'HC' };
fcfg.cmp_nme = { { 'MCI' 'EPD_Old' } };

fcfg.con_spr_plt = 1;

fcfg.pvl_lvl = .001;

fcfg.grp_plt     = 1;
fcfg.reg_lne_plt = 0;
cfg.reg_srf_plt  = 0;

fcfg.ttl     = 'Modularity';
fcfg.sve_loc = {'/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision'};
fcfg.sve_nme = {'Modularity'};

fcfg.ylm     = [.00 .35];

fcfg.sig_loc = 'bottom';

fcfg.hme_wrk = 1;

GraphLinePlot(fcfg)

%% Figure 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TLE
thk_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
cfg.thk_nme = thk_dta(1,2:end);

load( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'TestResiduals' '/' 'Thickness_Desikan_Residuals' '/' 'GraphData' '_' 'LocalPathEfficiency' '.mat' ] );
cfg.grp_dta = plt_dta;

iGO = 2;
cmp_loc = 3;
iBC = 4;

cfg.pvl_lvl = .01;

% %%%%%%%%
plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};

% Regional Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Views
sph     = {'lh'  'rh'};
sph_vew = {'lat' 'med' 'ven'};

% Load Surface
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
srf_brn{1} = fs_read_surf([ '/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage' '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf([ '/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage' '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

% Colorbar
top_pct = 1;
col{1} = hsv2rgb([ .27  .02  .65 ]); % rgb('orange');
col{2} = hsv2rgb([ .27  .02  .65 ]); % rgb('neon red');
col{3} = hsv2rgb([ .27  .02  .65 ]); %rgb('red');
col{4} = hsv2rgb([ .27  .02  .65 ]); % rgb('dark red');
col{5} = hsv2rgb([ .27  .02  .65 ]); % rgb('medium grey')-0.15;
col{6} = hsv2rgb([ .00  .45  .80 ]); % rgb('dark blue');
col{7} = hsv2rgb([ .98  .60  .70 ]); % rgb('blue');
col{8} = hsv2rgb([ .96  .75  .60 ]); % rgb('bright blue');
col{9} = hsv2rgb([ .92  1.0  .35 ]); % rgb('cyan');
col_map = [];
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
end
col_map = flipud(col_map);

% Plot
fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);

for iH = 1:numel(sph)
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,albl,~]=fs_read_annotation([ '/home/ekaestne/PROJECTS/DATA/fsaverage/ROIs/' sph{iH} '.aparc'  '.split' '.annot']);
    albl_hld = albl;
    albl = cellfun( @(x) mmil_spec_char( x , {'-'}) , albl , 'uni' , 0 );
    pct_hld = [albl num2cell(nan(size(albl,1),1))];
    
    for iF = 1:size(pct_hld,1)
        
        roi_loc = find(strcmpi( cfg.thk_nme , [ sph{iH} 's' '_'  pct_hld{iF,1} ] ));
        
        con_hld = nanmean( squeeze(cfg.grp_dta.con_spr_reg{iGO}( : , roi_loc , iBC )),2); % cfg.reg_col
        dta_hld = nanmean( cfg.grp_dta.dta{iGO}{cmp_loc}( roi_loc , iBC )); % cfg.reg_col
        
        if ~isempty(roi_loc)
            
            pct_plt(1) = prctile( con_hld , 100-((cfg.pvl_lvl/2)*100));
            pct_plt(2) = prctile( con_hld , (cfg.pvl_lvl/2)*100);
            
            if pct_plt(1) < dta_hld || ...
                    pct_plt(2) > dta_hld
                
                men_roi = nanmean( con_hld );
                std_roi = nanstd( con_hld );
                pct_hld{iF,2} = ( dta_hld - men_roi ) / std_roi;
                
            else
                
                pct_hld{iF,2} = 0;
                
            end
            
        else
            pct_hld{iF,2} = nan;
        end
    end
    
    pct_hld( cell2mat(pct_hld(:,2))>0 , 2 )  = {0};
    pct_hld( cell2mat(pct_hld(:,2))<-5 , 2 ) = {-10};
    
    pct_hld(:,2) = num2cell( (cell2mat(pct_hld(:,2)) / 20) + .5);
    
    pct_hld(:,1) = albl_hld;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(cfg,'reg_inc_plt') && ~isempty(cfg.reg_inc_plt)
        iIN = ismember(pct_hld(:,1),unique(cellfun(@(x) x(4:end),cfg.reg_inc_plt,'uni',0)) );
        pct_hld = pct_hld(iIN,:);
    end
    
    for iSP = 1:numel(sph_vew)
               
        pcfg = [];
        
        pcfg.surf_brain  = srf_brn{iH};
        
        pcfg.aparc       = [ '/home/ekaestne/PROJECTS' '/' 'DATA' '/' 'fsaverage' '/' 'ROIs' '/' sph{iH} '.aparc'  '.split' '.annot'];
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iSP};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = col_map;
        
        pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
        pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
        
        pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
    end
end

% Add Colorbar
axes('OuterPosition',[.9 .05 .05 .9],'visible','off')
colormap(col_map)
clb = colorbar('v','Position',[.9 .05 .05 .9]);
clb.TickLength = 0;
clb.TickLabels = cellfun(@num2str,num2cell(roundsd(linspace(-5,5,11),2)),'uni',0);
ylabel(clb,'LocalPathEfficiency','FontSize',20)

% Save
print(gcf,[ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision' '/' 'tle' '_' 'region_surface_plot' '_' 'density' '25' '_colorplay' '.png'],'-dpng','-r200')
close all

%% MCI
thk_dta = mmil_readtext( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'ROI' '/' 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv' ] );
cfg.thk_nme = thk_dta(1,2:end);

load( [ '/home/ekaestne/PROJECTS' '/' 'OUTPUT' '/' 'Epilepsy_and_Aging' '/' 'TestResiduals' '/' 'Thickness_Desikan_Residuals' '/' 'GraphData' '_' 'LocalPathEfficiency' '.mat' ] );
cfg.grp_dta = plt_dta;

iGO = 2;
cmp_loc = 2;
iBC = 4;

cfg.pvl_lvl = .01;

% %%%%%%%%
plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};

% Regional Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Views
sph     = {'lh'  'rh'};
sph_vew = {'lat' 'med' 'ven' };

% Load Surface
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
srf_brn{1} = fs_read_surf([ '/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage' '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf([ '/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage' '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

% Colorbar
top_pct = 1;
col{1} = hsv2rgb([ .27  .02  .65 ]); % rgb('orange');
col{2} = hsv2rgb([ .27  .02  .65 ]); % rgb('neon red');
col{3} = hsv2rgb([ .27  .02  .65 ]); %rgb('red');
col{4} = hsv2rgb([ .27  .02  .65 ]); % rgb('dark red');
col{5} = hsv2rgb([ .27  .02  .65 ]); % rgb('medium grey')-0.15;
col{6} = hsv2rgb([ .00  .45  .80 ]); % rgb('dark blue');
col{7} = hsv2rgb([ .98  .60  .70 ]); % rgb('blue');
col{8} = hsv2rgb([ .96  .75  .60 ]); % rgb('bright blue');
col{9} = hsv2rgb([ .92  1.0  .35 ]); % rgb('cyan');
col_map = [];
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
end
col_map = flipud(col_map);

% Plot
fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);

for iH = 1:numel(sph)
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,albl,~]=fs_read_annotation([ '/home/ekaestne/PROJECTS/DATA/fsaverage/ROIs/' sph{iH} '.aparc'  '.split' '.annot']);
    albl_hld = albl;
    albl = cellfun( @(x) mmil_spec_char( x , {'-'}) , albl , 'uni' , 0 );
    pct_hld = [albl num2cell(nan(size(albl,1),1))];
    
    for iF = 1:size(pct_hld,1)
        
        roi_loc = find(strcmpi( cfg.thk_nme , [ sph{iH} 's' '_'  pct_hld{iF,1} ] ));
        
        con_hld = nanmean( squeeze(cfg.grp_dta.con_spr_reg{iGO}( : , roi_loc , iBC )),2); % cfg.reg_col
        dta_hld = nanmean( cfg.grp_dta.dta{iGO}{cmp_loc}( roi_loc , iBC )); % cfg.reg_col
        
        if ~isempty(roi_loc)
            
            pct_plt(1) = prctile( con_hld , 100-((cfg.pvl_lvl/2)*100));
            pct_plt(2) = prctile( con_hld , (cfg.pvl_lvl/2)*100);
            
            if pct_plt(1) < dta_hld || ...
                    pct_plt(2) > dta_hld
                
                men_roi = nanmean( con_hld );
                std_roi = nanstd( con_hld );
                pct_hld{iF,2} = ( dta_hld - men_roi ) / std_roi;
                
            else
                
                pct_hld{iF,2} = 0;
                
            end
            
        else
            pct_hld{iF,2} = nan;
        end
    end
    
    pct_hld( cell2mat(pct_hld(:,2))>0 , 2 )  = {0};
    pct_hld( cell2mat(pct_hld(:,2))<-5 , 2 ) = {-10};
    
    pct_hld(:,2) = num2cell( (cell2mat(pct_hld(:,2)) / 20) + .5);
    
    pct_hld(:,1) = albl_hld;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(cfg,'reg_inc_plt') && ~isempty(cfg.reg_inc_plt)
        iIN = ismember(pct_hld(:,1),unique(cellfun(@(x) x(4:end),cfg.reg_inc_plt,'uni',0)) );
        pct_hld = pct_hld(iIN,:);
    end
    
    for iSP = 1:numel(sph_vew)
               
        pcfg = [];
        
        pcfg.surf_brain  = srf_brn{iH};
        
        pcfg.aparc       = [ '/home/ekaestne/PROJECTS' '/' 'DATA' '/' 'fsaverage' '/' 'ROIs' '/' sph{iH} '.aparc'  '.split' '.annot'];
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iSP};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iSP},'visible','off','Parent',fig_hld(1));
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = col_map;
        
        pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
        pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
        
        pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
    end
end

% Add Colorbar
axes('OuterPosition',[.9 .05 .05 .9],'visible','off')
colormap(col_map)
clb = colorbar('v','Position',[.9 .05 .05 .9]);
clb.TickLength = 0;
clb.TickLabels = cellfun(@num2str,num2cell(roundsd(linspace(-5,5,11),2)),'uni',0);
ylabel(clb,'LocalPathEfficiency','FontSize',20)

% Save
print(gcf,[ '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/epd_age_cnn/figures_tables/revision' '/' 'mci' '_' 'region_surface_plot' '_' 'density' '25' '_colorplay' '.png'],'-dpng','-r200')
close all
