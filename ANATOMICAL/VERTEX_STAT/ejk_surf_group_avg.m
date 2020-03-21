% % Setup
% clear; clc;
% 
% % Folders
% prj_dir = '/home/ekaestne/PROJECTS/';
% prj_nme = 'AnatomicalPhenotype';
% fsr_nme = 'sbj001_Freesurfer_Recons_uptodate.csv';
% 
% % Freesurfer Locations
% fsr_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' fsr_nme]);
%     fsr_fle = fsr_fle(:,1:2);
%     fsr_fle(cellfun(@isempty,fsr_fle)) = {''};
% 
% % Participant Group
% grp_fle = 'subjects_phenotype_update2_update.csv';
% grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
%     grp_fle(cellfun(@isempty,grp_fle)) = {''};
% 
% % Get group
% cfg = [];
% cfg.grp_fle = grp_fle;
% grp_str = ejk_collate_groups(cfg);
% 
% % Call function
% cfg = [];
% 
% cfg.prj_dir = '/home/ekaestne/PROJECTS/';
% cfg.prj_nme = 'AnatomicalPhenotype';
% 
% cfg.sbj_grp = grp_str;
% cfg.sbj_grp_col = { 'Phenotype' };
% cfg.sbj_grp_ind = { [2 3 4 5 6] };
% cfg.sbj_grp_cmp = { {[2 3] [2 4] [2 5] [2 6]} };
% 
% cfg.mes_typ     = 'wmparc_md'; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'
% 
% cfg.hms     = {'lhs' 'rhs'};

function ejk_surf_group_avg(cfg) % ProjID,analysis_outdir,cfg.mes_typ,analysis_dtidir)

%% FUNCTION BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if ~isfield(cfg,'hms'); cfg.hms = {'lhs' 'rhs'}; end

%%
fcfg = [];
fcfg.grp_fle = cfg.sbj_grp;
grp_fle      = ejk_collate_groups(fcfg);

%% Get and Collate Data
for iHM = 1:numel(cfg.hms)
    
    hms_nme = cfg.hms{iHM};
    
    for iGR = 1:numel(cfg.sbj_grp_col)
        
        grp_ovr_nme = cfg.sbj_grp_col{iGR};
        
        for iCD = 1:numel(cfg.sbj_grp_nme{iGR})
            
            grp_cde_nme = cfg.sbj_grp_nme{iGR}{iCD};
            grp_cde_ind = find( strcmpi( grp_fle.table.(grp_ovr_nme)(:,1) , cfg.sbj_grp_nme{iGR}{iCD} ) );
            grp_cde_sbj = grp_fle.sbj_nme( grp_fle.(grp_ovr_nme) == grp_cde_ind );            
            [ ~ , grp_dta_ind ] = intersect( cfg.dta{iHM}.srf_dta_sbj , grp_cde_sbj );
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).sbj_nme = cfg.dta{iHM}.srf_dta_sbj(grp_dta_ind);
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta  = cfg.dta{iHM}.srf_dta(grp_dta_ind,:);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum    = numel(grp_dta_ind) - sum(isnan(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta(:,1))); 
                         
            fprintf('N for diagnosis %s, hemi %s: %d\n',grp_cde_nme,hms_nme,grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
            vec_sum = nansum( grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta , 1 );
            vec_sm2 = nansum( grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).dta.^2 , 1 );
            vec_num = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum ;
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecmean   = vec_sum / vec_num;
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd    = sqrt( ( vec_num * vec_sm2 - vec_sum.^2 ) ./ ( vec_num * (vec_num - 1) ) );
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstderr = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd / sqrt(vec_num);
            
        end
    end
end

%% Differences
plt_nme_ind = 1;

for iGR = 1:numel(cfg.sbj_grp_col)
    
    grp_ovr_nme = cfg.sbj_grp_col{iGR};
    
    out_dir = [ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'surface' '/' cfg.dta_lbl{iGR} '/' grp_ovr_nme ];
    ejk_chk_dir( out_dir );
    
    for iC = 1:numel(cfg.sbj_grp_cmp{iGR})
        
        grp_cde_one_nme = cfg.sbj_grp_nme{iGR}{ cfg.sbj_grp_cmp{iGR}{iC}(1) };
        grp_cde_two_nme = cfg.sbj_grp_nme{iGR}{ cfg.sbj_grp_cmp{iGR}{iC}(2) };
        
        % absolute difference
        vecdiff   = grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
        vecstderr = sqrt(grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
        plt_men_dat{1}{plt_nme_ind}   = vecdiff;
        plt_tvl_dat{1}{plt_nme_ind}   = vecdiff ./ vecstderr;
        
        vecdiff   = grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
        vecstderr = sqrt(grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(cfg.hms{2}).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
        plt_men_dat{2}{plt_nme_ind}   = vecdiff;
        plt_tvl_dat{2}{plt_nme_ind}   = vecdiff ./ vecstderr;
        
        srf_typ_num(plt_nme_ind)     = grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_one_nme).nsum + grp_dta.(cfg.hms{1}).(grp_ovr_nme).(grp_cde_two_nme).nsum - 2;
        srf_typ_out_men_pre{plt_nme_ind} = sprintf('Cort_%s_diff_mean_%s_%s.mgh',cfg.dta_lbl{iGR},grp_cde_one_nme,grp_cde_two_nme);
        srf_typ_out_tvl_pre{plt_nme_ind} = sprintf('Cort_%s_diff_tval_%s_%s.mgh',cfg.dta_lbl{iGR},grp_cde_one_nme,grp_cde_two_nme);
        srf_typ_out_dir{plt_nme_ind} = out_dir;
        plt_nme_ind = plt_nme_ind + 1;
        
    end
    
end

%% Mean-Diff Plots
for iP = 1:numel(srf_typ_out_men_pre)
    
    pcfg = [];
    
    pcfg.out_dir     = srf_typ_out_dir{iP};
    pcfg.out_pre_fix = srf_typ_out_men_pre{iP};
    
    pcfg.plt_dta = { plt_men_dat{1}{iP} plt_men_dat{2}{iP} };
    
    pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
    pcfg.low_rng_num = [ -0.05 0.05 ];
    pcfg.hgh_rng_num = [ -0.25 0.25 ];
       
    mmil_anat_surf_plot(pcfg)
    
end

%% T-val Plots
for iP = 1:numel(srf_typ_out_tvl_pre)
    
    pcfg = [];
    
    pcfg.out_dir     = srf_typ_out_dir{iP};
    pcfg.out_pre_fix = srf_typ_out_tvl_pre{iP};
    
    pcfg.plt_dta = { plt_tvl_dat{1}{iP} plt_tvl_dat{2}{iP} };
    
    pcfg.fmr_col_map = {'cyan' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
        
    pcfg.low_rng_num = [ -2.365 2.365 ];
    pcfg.hgh_rng_num = [ -5.5 5.5 ];
       
    mmil_anat_surf_plot(pcfg)
    
end

%%
for iP = 1:numel(srf_typ_out_tvl_pre)

ejk_chk_dir( [ srf_typ_out_dir{iP} '/' 'mgz'] )


end

end

















