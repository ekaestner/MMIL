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

cfg.fmr_stt = { 'N-FF_GLT#0_Tstat'  'N#0_Tstat'         'FF#0_Tstat' };
cfg.fmr_nme = { 'N_FF'              'NW'                'FF' };
cfg.fmr_typ = {'-c'};

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

cell2csv([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_fsr_fle(1:end-4) '_' 'uptodate' '.csv'],cfg.sbj_fsr)

%% Cortical Thickness
prc_nme = { '' '.a2009s' }; % '' '.a2009s'

tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.ovr_dir = cfg.mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iPR = 1:numel(cfg.lbl_nme)
    fcfg.prc_nme = prc_nme{iPR};
    
    for iS = 1:size(cfg.sbj_fsr,1)
        fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
        fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
        [ gry_thk_dta_hld , tot_lbl ]= ejk_extract_grey_thickness(fcfg);
        if ~isempty(gry_thk_dta_hld); gry_thk_dta(iS,:) = gry_thk_dta_hld; else gry_thk_dta(iS,:) = nan(1,size(gry_thk_dta,2)); end
    end
    
    sve_sbj_nme = cfg.sbj_fsr(:,1);
    cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'MRI' '_' 'thickness' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(gry_thk_dta,1))] [tot_lbl ; num2cell(gry_thk_dta)] ]);
    clear gry_thk_dta tot_lbl
    
end
toc