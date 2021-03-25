%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define
cfg = [];

cfg.prj_dir = '/home/ekaestne/PROJECTS/';

cfg.sbj_nme_fle = 'sbj000_total_2019_03_27.csv';
cfg.sbj_fsr_fle = 'sbj001_Freesurfer_Recons.csv';
cfg.bld_sbj_dir = 'sbj006_Bold_Data_ACM.csv';

cfg.mri_fsr_dir = '/home/mmilmcd/data/FSRECONS/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
cfg.dti_fsr_dir = '/home/mmilmcdRSI/data/proc_dti/';
cfg.fmr_fsr_dir = '/home/mmilmcdRSI/data/fsurf/';

cfg.bld_dir     = '/home/mmilmcd/data/MCD_BOLD/subjects/';

cfg.hms     = {'lhs' 'rhs'};

cfg.lbl_nme = {'.a2009s' }; % ''

%% List of Subjects Put Together
sbj_nme_hld = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_nme_fle]);
sbj_fsr_hld = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_fsr_fle]);
sbj_fsr_hld(1,cellfun(@isempty,sbj_fsr_hld(1,:))) = {'PotentialRecon'};

dta_typ     = { 'T1' 'DTI' 'fMRI' };
dta_typ_col = find(ismember(sbj_fsr_hld(1,:),dta_typ));

fmr_sbj_dir = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.bld_sbj_dir]);

sbj_nme = sbj_nme_hld(2:end,1);
sbj_fsr = cell(numel(sbj_nme),numel(dta_typ)+1);
for iSB = 1:numel(sbj_nme)
    
    sbj_fsr{iSB,1} = sbj_nme{iSB};
    
    row_ind = find(strcmpi(sbj_fsr_hld(:,1),sbj_nme{iSB}));
    if ~isempty(row_ind)
        sbj_fsr(iSB,2:end) = sbj_fsr_hld(row_ind,dta_typ_col);
    end
end

cfg.sbj_fsr = sbj_fsr(:,1:2);
cfg.sbj_fsr(cellfun(@isempty,sbj_fsr(:,2)),2) = {''};
cfg.sbj_fsr(:,2) = regexprep(cfg.sbj_fsr(:,2),'FSURF_','');
cfg.sbj_fsr(:,2) = regexprep(cfg.sbj_fsr(:,2),'\..*','');

%% Cortical Thickness
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

fcfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 2819; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

fcfg.anl_dir = '';

for iH = 1:numel(cfg.hms)
    srf_dta = nan(size(cfg.sbj_fsr,1),163842);
    for iS = 1:size(cfg.sbj_fsr,1)
        fcfg.hms         = cfg.hms{iH}(1:2);
        fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
        srf_dta_hld = ejk_extract_vertices(fcfg);
        if ~isempty(srf_dta_hld); srf_dta(iS,:) = srf_dta_hld; end
    end
    srf_dta_sbj = cfg.sbj_fsr;
%     save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end

