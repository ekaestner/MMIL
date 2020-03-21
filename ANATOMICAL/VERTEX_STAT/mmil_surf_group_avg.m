%

function mmil_surf_group_avg(fcfg) % ProjID,analysis_outdir,fcfg.mes_typ,analysis_dtidir)

% Setup
clear; clc;

% Folders
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'AnatomicalPhenotype';
fsr_nme = 'sbj001_Freesurfer_Recons_uptodate.csv';

% Freesurfer Locations
fsr_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' fsr_nme]);
    fsr_fle = fsr_fle(:,1:2);
    fsr_fle(cellfun(@isempty,fsr_fle)) = {''};

% Participant Group
grp_fle = 'subjects_phenotype_update2_update.csv';
grp_fle = mmil_readtext([prj_dir '/' 'SUBJECTS' '/' 'tmp' '/' prj_nme '/' grp_fle]);
    grp_fle(cellfun(@isempty,grp_fle)) = {''};

% Get group
fcfg = [];
fcfg.grp_fle = grp_fle;
grp_str = ejk_collate_groups(fcfg);

% Call function
fcfg = [];

fcfg.prj_dir = '/home/ekaestne/PROJECTS/';
fcfg.prj_nme = 'AnatomicalPhenotype';

fcfg.sbj_grp = grp_str;
fcfg.sbj_grp_col = { 'Phenotype' };
fcfg.sbj_grp_ind = { [2 3 4 5 6] };
fcfg.sbj_grp_cmp = { {[2 3] [2 4] [2 5] [2 6]} };

fcfg.mes_typ     = 'wmparc_md'; % 'aMRI_thickness' 'wmparc_md' 'rsfMRI_variance'

fcfg.hms     = {'lhs' 'rhs'};

%% FUNCTION BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if ~isfield(fcfg,'hms'); fcfg.hms = {'lhs' 'rhs'}; end

if ~isdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'surface' '/' fcfg.mes_typ]); mkdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'surface' '/' fcfg.mes_typ]); end
for iGR = 1:numel(fcfg.sbj_grp_col); if ~isdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'surface' '/' fcfg.mes_typ '/' mmil_spec_char(fcfg.sbj_grp_col{iGR},{' ','&'})]); mkdir([fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'surface' '/' fcfg.mes_typ '/' mmil_spec_char(fcfg.sbj_grp_col{iGR},{' ','&'})]); end; end

%% Get and Collate Data
%
for iHM = 1:numel(fcfg.hms)
   load([fcfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms{iHM}])
   ind = nan(size(fcfg.sbj_grp.sbj_nme,1),1);
   for iSB = 1:size(fcfg.sbj_grp.sbj_nme,1)
       ind(iSB) = find(strcmpi(srf_dta_sbj(:,1),fcfg.sbj_grp.sbj_nme{iSB,1}));
   end
   dta_hld.(fcfg.hms{iHM}).dta     = srf_dta(ind,:);
   dta_hld.(fcfg.hms{iHM}).sbj_nme = srf_dta_sbj(ind,:);
end

%
for iHM = 1:numel(fcfg.hms)
    
    hms_nme = fcfg.hms{iHM};
    
    for iGR = 1:numel(fcfg.sbj_grp_col)
        
        grp_ovr_nme = mmil_spec_char(fcfg.sbj_grp_col{iGR},{' ','&'});
        
        for iCD = 1:numel(fcfg.sbj_grp_ind{iGR})
            
            grp_cde_nme = mmil_spec_char(fcfg.sbj_grp.code.(fcfg.sbj_grp_col{iGR}){fcfg.sbj_grp_ind{iGR}(iCD)},{' ','&'});
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum    = sum(fcfg.sbj_grp.(grp_ovr_nme)==fcfg.sbj_grp_ind{iGR}(iCD));
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum  = nansum(dta_hld.(hms_nme).dta(fcfg.sbj_grp.(grp_ovr_nme)==fcfg.sbj_grp_ind{iGR}(iCD),:),1);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum2 = nansum(dta_hld.(hms_nme).dta(fcfg.sbj_grp.(grp_ovr_nme)==fcfg.sbj_grp_ind{iGR}(iCD),:).^2,1);
            
            fprintf('N for diagnosis %s, hemi %s: %d\n',grp_cde_nme,hms_nme,grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecmean   = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum / (eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd    = sqrt((grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum*grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum2 - grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecsum.^2)./(eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum*(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum-1)));
            grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstderr = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).vecstd/sqrt(eps+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_nme).nsum);
            
        end
    end
end

%% Differences
plt_nme_ind = 1;
for iHM = 1:numel(fcfg.hms)
    
    hms_nme = fcfg.hms{iHM};
    
    for iGR = 1:numel(fcfg.sbj_grp_col)
        
        out_dir = [fcfg.prj_dir '/' 'OUTPUT' '/' fcfg.prj_nme '/' 'surface' '/' fcfg.mes_typ '/' mmil_spec_char(fcfg.sbj_grp_col{iGR},{' ','&'})];
        
        grp_ovr_nme = mmil_spec_char(fcfg.sbj_grp_col{iGR},{' ','&'});
        
        for iC = 1:numel(fcfg.sbj_grp_cmp{iGR})
            
            i1 = fcfg.sbj_grp_cmp{iGR}{iC}(1); 
                grp_cde_one_nme = mmil_spec_char(fcfg.sbj_grp.code.(fcfg.sbj_grp_col{iGR}){i1},{' ','&'});
            i2 = fcfg.sbj_grp_cmp{iGR}{iC}(2);
                grp_cde_two_nme = mmil_spec_char(fcfg.sbj_grp.code.(fcfg.sbj_grp_col{iGR}){i2},{' ','&'});
            
            % absolute difference
            vecdiff   = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).vecmean - grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).vecmean;
            vecstderr = sqrt(grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).vecstderr.^2+grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).vecstderr.^2);
            vectval   = vecdiff./vecstderr;
            
            fname = sprintf('%s/Cort_AbsDff_%s_diff_mean_%s_%s-%s.mgh',out_dir,fcfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
            fs_save_mgh(vecdiff,fname,eye(4));
            
            fname = sprintf('%s/Cort_AbsDff_%s_diff_tval_%s_%s-%s.mgh',out_dir,fcfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
            fs_save_mgh(vectval,fname,eye(4));
            
            fname = sprintf('%s/Cort_AbsDff_%s_diff_stderr_%s_%s-%s.mgh',out_dir,fcfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
            fs_save_mgh(vecstderr,fname,eye(4));
            
            if iHM==1
                srf_typ_nme{plt_nme_ind}     = sprintf('%s/Cort_AbsDff_%s_diff_tval_%s_%s-%s.mgh',out_dir,fcfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme,hms_nme);
                srf_typ_num(plt_nme_ind)     = grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_one_nme).nsum + grp_dta.(hms_nme).(grp_ovr_nme).(grp_cde_two_nme).nsum - 2;
                srf_typ_out_pre{plt_nme_ind} = sprintf('Cort_AbsDff_%s_diff_tval_%s_%s.mgh',fcfg.mes_typ,grp_cde_one_nme,grp_cde_two_nme);
                srf_typ_out_dir{plt_nme_ind} = out_dir;
                plt_nme_ind = plt_nme_ind + 1;
            end
            
        end
        
    end
    
end

%% Plots
for iP = 1:numel(srf_typ_nme)

    pcfg = [];
    
    pcfg.out_dir = srf_typ_out_dir{iP};
    pcfg.dta_dir = [];
    pcfg.out_pre_fix = srf_typ_out_pre{iP};
    
    pcfg.srf_typ_nme = [ srf_typ_nme(iP) regexprep(srf_typ_nme(iP),'lhs','rhs') ];
    
    pcfg.fmr_col_map = {'neon blue' 'blue' 'dark blue' 'grey' 'dark red' 'red' 'bright yellow'};
    pcfg.fmr_rng_num = [];
        
    pcfg.deg_fre = srf_typ_num(iP);
       
    mmil_anat_surf_plot(pcfg)
    
end

end

















