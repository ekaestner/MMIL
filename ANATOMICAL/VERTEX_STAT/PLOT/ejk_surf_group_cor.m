% % Setup
% clear; clc;
%
% % Folders
% prj_dir = '/home/ekaestne/PROJECTS/';
% prj_nme = 'Tests';
% fsr_nme = 'sbj001_Freesurfer_Recons.csv';
%
% % Freesurfer Locations
% fsr_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' fsr_nme]);
%     fsr_fle = fsr_fle(:,1:2);
%     fsr_fle(cellfun(@isempty,fsr_fle)) = {''};
%
% % Participant Group
% grp_fle = 'allsubjects.csv';
% grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
%     grp_fle(cellfun(@isempty,grp_fle)) = {''};
%
% % Get group
% cfg = [];
% cfg.grp_fle = grp_fle;
% grp_str = ejk_collate_groups(cfg);
%
% % REDCAP
% red_fle = 'sbj000_total_2019_03_27.csv';
%
% cfg = [];
% cfg.prj_dir = prj_dir;
% cfg.red_fle = red_fle;
% cfg.sbj_nme = grp_fle;
% [sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(cfg);
%
% % Call function
% cfg = [];
%
% cfg.prj_dir = prj_dir;
% cfg.prj_nme = prj_nme;
%
% cfg.sbj_grp     = grp_str;
% cfg.sbj_grp_col = { 'phenotype' };
% cfg.sbj_grp_ind = { [1] };
%
% cfg.mes_typ     = 'wmparc_md'; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'
%
% cfg.cov         = { sbj_dem           sbj_dem }; % { sbj_cog      sbj_cog           sbj_cog        sbj_cog               sbj_cog };
% cfg.cov_nme     = { 'sbj_age' 'sbj_edu' }; % {'vp2_nor_scr' 'cat_flu_nor_scr' 'bnt_nor_scr'  'log_mem_nor_scr_one' 'ltr_tot_nor_scr' };
%
% cfg.hms     = {'lhs' 'rhs'};

function ejk_surf_group_cor(cfg) % ProjID,analysis_outdir,cfg.mes_typ,analysis_dtidir)

%% FUNCTION BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg,'hms'); cfg.hms = {'lhs' 'rhs'}; end

% Get group
fcfg = [];
fcfg.grp_fle = cfg.sbj_grp;
cfg.grp_str  = ejk_collate_groups(fcfg);

if ~isdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'correlation' '/' cfg.mes_typ]); mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'correlation' '/' cfg.mes_typ]); end
for iGR = 1:numel(cfg.sbj_grp_col); if ~isdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'correlation' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})]); mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'correlation' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})]); end; end

%% Get and Collate Data
%
for iHM = 1:numel(cfg.hms)
    load([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' cfg.mes_typ '_' cfg.hms{iHM}])
    ind = nan(size(cfg.grp_str.sbj_nme,1),1);
    for iSB = 1:size(cfg.grp_str.sbj_nme,1)
        ind(iSB) = find(strcmpi(srf_dta_sbj(:,1),cfg.grp_str.sbj_nme{iSB,1}));
    end
    dta_hld.(cfg.hms{iHM}).dta     = srf_dta(ind,:);
    dta_hld.(cfg.hms{iHM}).sbj_nme = srf_dta_sbj(ind,:);
end

% Gather Data
for iCV = 1:numel(cfg.cov)
    cov_fld_nme = fieldnames(cfg.cov{iCV});
    for iCN = 1:numel(cov_fld_nme)
        if isnumeric(cfg.cov{iCV}.(cov_fld_nme{iCN}))
            for iSB = 1:size(dta_hld.(cfg.hms{iHM}).sbj_nme,1)
                cov_dta.(cov_fld_nme{iCN})(iSB,1) = cfg.cov{iCV}.(cov_fld_nme{iCN})(strcmpi(cfg.cov{iCV}.sbj_nme,dta_hld.(cfg.hms{1}).sbj_nme{iSB,1}),1);
            end
        end
    end
end

%
plt_dta = cell(0);
plt_nme = cell(0);
cov_fld_nme = fieldnames(cov_dta);

for iHM = 1:numel(cfg.hms)
    hms_nme = cfg.hms{iHM};
    
    for iGR = 1:numel(cfg.sbj_grp_col)
        sbj_grp_col = cfg.sbj_grp_col{iGR};
        
        for iCV = 1:numel(cov_fld_nme)
            cov_nme = cov_fld_nme{iCV};
            
            for iCD = 1:numel(cfg.sbj_grp_ind{iGR})
                
                grp_cde_nme = mmil_spec_char(cfg.grp_str.code.(cfg.sbj_grp_col{iGR}){cfg.sbj_grp_ind{iGR}(iCD)},{' ','&'});
                
                sbj_ind = find(cfg.grp_str.(sbj_grp_col)==cfg.sbj_grp_ind{iGR}(iCD));
                sbj_ind( isnan(dta_hld.(hms_nme).dta(sbj_ind,1)) | isnan(cov_dta.(cov_fld_nme{iCV})(sbj_ind,1)) ) = [];
                
                grp_dta.(hms_nme).(sbj_grp_col).(grp_cde_nme).cor = zeros(1,size(dta_hld.(hms_nme).dta,2));
                grp_dta.(hms_nme).(sbj_grp_col).(grp_cde_nme).pvl = zeros(1,size(dta_hld.(hms_nme).dta,2));
                
                for iC = 1:size(dta_hld.(hms_nme).dta,2)
                    
                    [cor_hld,pvl_hld] = corrcoef( [ dta_hld.(hms_nme).dta(sbj_ind,iC)  ...
                        dta_hld.cov(sbj_ind,iCV) ]);
                    
                    grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).cor(1,iC) = cor_hld(1,2);
                    grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).pvl(1,iC) = pvl_hld(1,2);
                    
                end
                
                
                if iHM == 1
                    plt_nme{end+1} = [ grp_cde_nme '_FOR_' cfg.mes_typ '_BY_' cfg.cov_nme{iCV}];
                    plt_dta{end+1}{iHM} = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).pvl;
                    plt_dta{end}{iHM}(isnan(plt_dta{end}{iHM})) = 1;
                else
                    plt_dta{end}{iHM} = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).pvl;
                    plt_dta{end}{iHM}(isnan(plt_dta{end}{iHM})) = 1;
                end
                
            end
        end
        
    end
end

%% Plots
for iP = 1:numel(plt_nme)
    
    pcfg = [];
    
    pcfg.out_dir     = [cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'correlation' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})];
    pcfg.out_pre_fix = plt_nme{iP};
    
    pcfg.plt_dta = { plt_dta{iP}{1} plt_dta{iP}{2} };
    
    pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    pcfg.pvl_rng_num = [0 0.05];
    
    mmil_anat_surf_plot(pcfg)
    
end

end



