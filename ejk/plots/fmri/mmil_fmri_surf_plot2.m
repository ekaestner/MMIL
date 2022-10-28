% fcfg = [];
% 
% fcfg.prj_dir = '/home/ekaestne/PROJECTS/';
% fcfg.prj_nme = 'CognitivePhenotype';
% 
% fcfg.grp_dir = '/space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/';
% fcfg.grp_nme = {'ALL_42controls_BETAS_Aug31+tlrc' ...
%                 'ALL_15NO_IMP_BETAS_Aug31+tlrc' ...
%                 'ALL_16LANG_BETAS_Aug31+tlrc' ...
%                 'ALL_10MEMORY_BETAS_Aug31+tlrc' ...
%                 'ALL_10LMs_BETAS_Aug31+tlrc'};
% fcfg.plt_nme = {'Control' ...
%                 'NoImpairment' ...
%                 'LanguageImpairment' ...
%                 'MemoryImpairment' ...
%                 'LanguageAndMemoryImpairment'};
% 
% fcfg.fmr_col_map = {'bright blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright red'};
% fcfg.fmr_rng_num = [-10 10];
% 
% fcfg.sph     = {'lh' 'rh'};
% fcfg.sph_vew = {'lat' 'med' 'ven'};
%
%  
% 3dVol2Surf  -out_1D /home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf/AllControl_lhs.1D
% 
% 3dVol2Surf -spec /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/suma_TT_N27_ROI/TT_N27_both.spec -surf_A rh.pial.gii -sv /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/suma_TT_N27_ROI/TT_N27+tlrc -grid_parent /space/syn09/1/data/MMILDB/MCD_BOLD/GroupAnalysis/08312018/'ALL_42controls_BETAS_Aug31+tlrc[7]' -map_func mask -out_1D /home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf/AllControl_rhs.1D
% %% Create Files
% for iH = 1:numel(hms)
%    
%     for iG = 1:numel(fcfg.grp_nme)
%         
%         cmd = '3dVol2Surf';
%         cmd = [cmd ' ' '-spec' ' ' fcfg.grp_dir '/' 'suma_TT_N27_ROI' '/' 'TT_N27_both.spec'];
%         cmd = [cmd ' ' '-surf_A ' hms{iH} '.pial.gii'];
%         cmd = [cmd ' ' '-sv' ' ' fcfg.grp_dir '/' 'suma_TT_N27_ROI' '/' 'TT_N27+tlrc'];
%         cmd = [cmd ' ' '-grid_parent' ' ' fcfg.grp_dir '/' '''' fcfg.grp_nme{iG} '[7]' '''' ];
%         cmd = [cmd ' ' '-map_func mask'];
%         cmd = [cmd ' ' '-out_1D' ' '  fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf' '/' 'data_hold' '/' fcfg.plt_nme{iG} '_' hms{iH} 's.1D'];
%         
%         out = unix(cmd);
%         
%     end
%     
%     cmd = ['mris_convert' ' ' fcfg.grp_dir '/' 'suma_TT_N27_ROI' '/' hms{iH} '.pial.gii'  ' ' fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf' '/' 'data_hold' '/' hms{iH} '.pial'];
%     out = unix(cmd);
%     
% end

function mmil_fmri_surf_plot2(fcfg)

%%
hms = {'lh' 'rh'};

if exist([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' ])~=7;                                 mkdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' ]); end
if exist([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf'])~=7;                 mkdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf']); end
if exist([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf' '/' 'data_hold'])~=7; mkdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf' '/' 'data_hold']); end

%% Plot
plt_loc = {{[0.0 0.6 0.4 0.4] [0.0 0.2 0.4 0.4] [0.0 0.0 0.4 0.2]} {[0.4 0.6 0.4 0.4] [0.4 0.2 0.4 0.4] [0.4 0.0 0.4 0.2]}};

% Make Colormap
fcfg.fmr_col_map = cellfun(@rgb ,fcfg.fmr_col_map,'uni',0);

fmr_col_map = [];
for iC = 1:numel(fcfg.fmr_col_map)-1
    fmr_col_map = [fmr_col_map ; [linspace(fcfg.fmr_col_map{iC}(1),fcfg.fmr_col_map{iC+1}(1),ceil(1000*top_pct/(numel(fcfg.fmr_col_map)-1)))' linspace(fcfg.fmr_col_map{iC}(2),fcfg.fmr_col_map{iC+1}(2),ceil(1000*top_pct/(numel(fcfg.fmr_col_map)-1)))' linspace(fcfg.fmr_col_map{iC}(3),fcfg.fmr_col_map{iC+1}(3),ceil(1000*top_pct/(numel(fcfg.fmr_col_map)-1)))']; ];
end

for iG = 1:numel(fcfg.grp_nme)
    
    fig_hld(1) = figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
    
    for iH = 1:numel(hms)
        
        % Get Data
        fmr_dat.(hms{iH}) = mmil_readtext(['/home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf' '/' 'data_hold' '/' fcfg.plt_nme{iG} '_' hms{iH} 's.1D'],' ','#');
        
        fmr_dat_hld.(hms{iH}) = fmr_dat.(hms{iH});
        fmr_dat_cln_sze.(hms{iH}) = size(fmr_dat_hld.(hms{iH}));
        
        fmr_dat_cln.(hms{iH}) = cell(fmr_dat_cln_sze.(hms{iH})(1)-4,7);
        for iR = 1:fmr_dat_cln_sze.(hms{iH})(1)-4
            fmr_dat_cln.(hms{iH})(iR,:) = fmr_dat_hld.(hms{iH})(iR+4,~cellfun(@isempty,fmr_dat_hld.(hms{iH})(iR+4,:)));
        end
        fmr_dat_cln.(hms{iH}) = cell2mat(fmr_dat_cln.(hms{iH}));
        
        min_loc = find(fmr_dat_cln.(hms{iH})(:,end)<fcfg.fmr_rng_num(1)); fmr_dat_cln.(hms{iH})(min_loc,end) = fcfg.fmr_rng_num(1)+abs(fcfg.fmr_rng_num(1)/100);
        max_loc = find(fmr_dat_cln.(hms{iH})(:,end)>fcfg.fmr_rng_num(2)); fmr_dat_cln.(hms{iH})(max_loc,end) = fcfg.fmr_rng_num(2);
        
        fmr_dat_cln.(hms{iH})(fmr_dat_cln.(hms{iH})(:,end)<2.8 & fmr_dat_cln.(hms{iH})(:,end)>-2.8,end) = 0.1;
        
        fmr_dat_cln.(hms{iH})(:,1) = fmr_dat_cln.(hms{iH})(:,1) + 1;
        
        fmr_val.(hms{iH}) = ( fmr_dat_cln.(hms{iH})(:,end) + abs(fcfg.fmr_rng_num(1)) ) / sum(abs(fcfg.fmr_rng_num));
        fmr_loc.(hms{iH}) = fmr_dat_cln.(hms{iH})(:,1);
        
        % Load Brain
        srf_brn{iH}.surf_brain =  fs_read_surf(['/home/ekaestne/PROJECTS/OUTPUT/CognitivePhenotype/FunctionalSurf' '/' 'data_hold' '/' hms{iH} '.pial']);
        srf_brn{iH}.surf_brain.coords = srf_brn{iH}.surf_brain.vertices;
        
        %
        for iP = 1:numel(plt_loc{iH})
            
            pcfg = [];
            
            pcfg.surf_brain  = srf_brn{iH};
            
            pcfg.sph         = fcfg.sph{iH};
            pcfg.sph_vew     = fcfg.sph_vew{iP};
            
            pcfg.label       = 0;
            pcfg.radius      = [];
            pcfg.alpha       = 1;
            
            pcfg.non_ele     = [];
            pcfg.sve_img     = 0; % ###
            
            pcfg.axe_hnd = axes('OuterPosition',plt_loc{iH}{iP},'visible','off','Parent',fig_hld(1));
            pcfg.fig_hdl = fig_hld(1);
            
            pcfg.col_map = fmr_col_map;
            
            pcfg.tbl_pct = fmr_val.(hms{iH});
            pcfg.tbl_loc = fmr_loc.(hms{iH});
            
            pcfg.top_pct = 1;
            
            nyu_plot2(pcfg);
            
        end
        
    end
    
    % Colorbar
    ax1 = axes('OuterPosition',[0.85 0.20 0.04 0.60],'visible','off','Parent',fig_hld(1));
    colormap(ax1,fmr_col_map)
    clb = colorbar('west','Position',[0.92 0.20 0.02 0.60]);
    clb.TickLength = 0;
    clb.TickLabels = cellfun(@num2str,num2cell(linspace(fcfg.fmr_rng_num(1),fcfg.fmr_rng_num(2),9)),'uni',0);
    
    % Save Figure
    print([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'FunctionalSurf' '/' fcfg.plt_nme{iG} '.png'],'-dpng','-r300')
    close all
    
end

end





