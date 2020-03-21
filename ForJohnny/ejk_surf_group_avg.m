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
% cfg.prj_dir = '/home/ekaestne/PROJECTS/';
% cfg.prj_nme = 'Tests';
% 
% cfg.mes_typ     = 'wmparc_md'; 
% 
% cfg.sbj_grp     = grp_fle;
% cfg.sbj_grp_col = { 'phenotype' };
% 
% cfg.sbj_grp_ind = { [1 2 3 4 5] };
% cfg.sbj_grp_cmp = { {[1 2] [1 3] [1 4] [1 5]} };
% 
% cfg.hms     = {'lhs' 'rhs'};

%% FUNCTION BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ejk_surf_group_avg(cfg)

if ~isfield(cfg,'hms'); cfg.hms = {'lhs' 'rhs'}; end

% Get group
fcfg = [];
fcfg.grp_fle = cfg.sbj_grp;
cfg.grp_str  = ejk_collate_groups(fcfg);

if ~isdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'group' '/' cfg.mes_typ]); mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'group' '/' cfg.mes_typ]); end
for iGR = 1:numel(cfg.sbj_grp_col); if ~isdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'group' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})]); mkdir([cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'group' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})]); end; end

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

%
for iHM = 1:numel(cfg.hms)
    
    hms_nme = cfg.hms{iHM};
    
    for iGR = 1:numel(cfg.sbj_grp_col)
        
        grp_ovr_nme = mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'});
        
        for iCD = 1:numel(cfg.sbj_grp_ind{iGR})
            
            grp_cde_nme = mmil_spec_char(cfg.grp_str.code.(cfg.sbj_grp_col{iGR}){cfg.sbj_grp_ind{iGR}(iCD)},{' ','&'});
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum    = sum(cfg.grp_str.(grp_ovr_nme)==cfg.sbj_grp_ind{iGR}(iCD));
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum  = nansum(dta_hld.(hms_nme).dta(cfg.grp_str.(grp_ovr_nme)==cfg.sbj_grp_ind{iGR}(iCD),:),1);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum2 = nansum(dta_hld.(hms_nme).dta(cfg.grp_str.(grp_ovr_nme)==cfg.sbj_grp_ind{iGR}(iCD),:).^2,1);
            
            fprintf('N for diagnosis %s, hemi %s: %d\n',grp_cde_nme,hms_nme,grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecmean   = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum / (eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd    = sqrt((grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum*grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum2 - grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum.^2)./(eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum*(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum-1)));
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstderr = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd/sqrt(eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
        end
    end
end

%% Differences
plt_nme_ind = 1;
for iHM = 1:numel(cfg.hms)
    
    plt_dat_ind = 1;
    
    hms_nme = cfg.hms{iHM};
    
    for iGR = 1:numel(cfg.sbj_grp_col)
        
        out_dir = [cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' 'group' '/' cfg.mes_typ '/' mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'})];
        
        grp_ovr_nme = mmil_spec_char(cfg.sbj_grp_col{iGR},{' ','&'});
        
        for iC = 1:numel(cfg.sbj_grp_cmp{iGR})
            
            i1 = cfg.sbj_grp_cmp{iGR}{iC}(1); 
                grp_cde_one_nme = mmil_spec_char(cfg.grp_str.code.(cfg.sbj_grp_col{iGR}){i1},{' ','&'});
            i2 = cfg.sbj_grp_cmp{iGR}{iC}(2);
                grp_cde_two_nme = mmil_spec_char(cfg.grp_str.code.(cfg.sbj_grp_col{iGR}){i2},{' ','&'});
            
            % absolute difference
            vecdiff   = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
            vecstderr = sqrt(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
            vectval   = vecdiff./vecstderr;
            
            plt_dta{plt_dat_ind}{iHM} = vectval;
            plt_dta{plt_dat_ind}{iHM}(isnan(plt_dta{plt_dat_ind}{iHM})) = 0;
            plt_dat_ind = plt_dat_ind + 1;
            
            if iHM==1
                srf_typ_nme{plt_nme_ind}     = sprintf('%s/Cort_AbsDff_%s_diff_tval_%s_%s-%s.mgh',out_dir,cfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
                srf_typ_num(plt_nme_ind)     = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).nsum + grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).nsum - 2;
                srf_typ_out_pre{plt_nme_ind} = sprintf('Cort_AbsDff_%s_diff_tval_%s_%s.mgh',cfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme);
                srf_typ_out_dir{plt_nme_ind} = out_dir;
                plt_nme_ind = plt_nme_ind + 1;
            end
            
        end
        
    end
    
end

%% FDR
% for iP = 1:numel(plt_dta)
%     pvl_hld = 1-tcdf(abs(plt_dta{iP}{1}),srf_typ_num(iP));    
% 	[pID,pN] = FDR(pvl_hld(pvl_hld~=.50),.05);
%     tinv(pID,srf_typ_num(iP))
% end

%% Plots
for iP = 1:numel(plt_dta)

    pcfg = [];
    
    pcfg.out_dir = srf_typ_out_dir{iP};
    pcfg.dta_dir = [];
    pcfg.out_pre_fix = srf_typ_out_pre{iP};
    
    pcfg.plt_dta = { plt_dta{iP}{1} plt_dta{iP}{2} };
    
    pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    pcfg.fmr_rng_num = [];
        
    pcfg.deg_fre = srf_typ_num(iP);
       
    mmil_anat_surf_plot(pcfg)
    
end

%%
% fname = sprintf('%s/Cort_AbsDff_%s_diff_mean_%s_%s-%s.mgh',out_dir,cfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
% fs_save_mgh(vecdiff,fname,eye(4));
% 
% fname = sprintf('%s/Cort_AbsDff_%s_diff_tval_%s_%s-%s.mgh',out_dir,cfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
% fs_save_mgh(vectval,fname,eye(4));
% 
% fname = sprintf('%s/Cort_AbsDff_%s_diff_stderr_%s_%s-%s.mgh',out_dir,cfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
% fs_save_mgh(vecstderr,fname,eye(4));

















