function ejk_extract_data

%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define
cfg = [];

cfg.prj_dir = '/home/ekaestne/PROJECTS/';

cfg.sbj_nme_fle = 'sbj000_total_2019_03_27.csv';
cfg.sbj_fsr_fle = 'sbj001_Freesurfer_Recons.csv';
cfg.bld_sbj_dir = 'sbj002_BOLD.csv';

cfg.mri_fsr_dir = '/home/mmilmcd/data/FSRECONS/'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
cfg.dti_fsr_dir = '/home/mmilmcdRSI/data/proc_dti/';
cfg.fmr_fsr_dir = '/home/mmilmcdRSI/data/fsurf/';
cfg.cnn_fsr_dir = '/space/syn09/1/data/MMILDB/ABALA/Connectome/connectomes';

cfg.bld_dir     = '/home/mmilmcd/data/MCD_BOLD/subjects/';

cfg.hms     = {'lhs' 'rhs'};

cfg.lbl_nme = {'.a2009s' }; % ''

cfg.fmr_stt = { 'N-FF_GLT#0_Tstat'  'N#0_Tstat'         'FF#0_Tstat' };
cfg.fmr_nme = { 'N_FF'              'NW'                'FF' };
cfg.fmr_typ = {'-c'};

%% List of Subjects Put Together
sbj_nme_hld = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_nme_fle]);
sbj_fsr = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_fsr_fle]);
sbj_fsr = sbj_fsr(:,[1 2]); sbj_fsr(cellfun(@isempty,sbj_fsr(:,2)),2) = {''};

dta_typ     = { 'T1' 'DTI' 'fMRI' };
dta_typ_col = find(ismember(sbj_fsr(1,:),dta_typ));

fmr_sbj_dir = mmil_readtext([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.bld_sbj_dir]);

sbj_nme = sbj_nme_hld(2:end,1);
% sbj_fsr = cell(numel(sbj_nme),numel(dta_typ)+1);
% for iSB = 1:numel(sbj_nme)
%     
%     sbj_fsr{iSB,1} = sbj_nme{iSB};
%     
%     row_ind = find(strcmpi(sbj_fsr_hld(:,1),sbj_nme{iSB}));
%     if ~isempty(row_ind)
%         sbj_fsr(iSB,2:end) = sbj_fsr_hld(row_ind,dta_typ_col);
%     end
% end

cfg.sbj_fsr = sbj_fsr(:,1:2);
cfg.sbj_fsr(cellfun(@isempty,sbj_fsr(:,2)),2) = {''};
cfg.sbj_fsr(:,2) = regexprep(cfg.sbj_fsr(:,2),'FSURF_','');
cfg.sbj_fsr(:,2) = regexprep(cfg.sbj_fsr(:,2),'\..*','');
% 
% cell2csv([cfg.prj_dir '/' 'SUBJECTS' '/' cfg.sbj_fsr_fle(1:end-4) '_' 'uptodate' '.csv'],cfg.sbj_fsr)

%% ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROI Volume Size
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iP = 1:numel(cfg.lbl_nme)
    fcfg.prc_nme = cfg.lbl_nme{iP};
    for iS = 1:size(cfg.sbj_fsr,1)
        fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
        vox_dta_hld = ejk_extract_roi_voxels(fcfg);
        if ~isempty(vox_dta_hld); vox_dta(iS,:) = vox_dta_hld; else vox_dta(iS,:) = nan(1,size(vox_dta,2)); end
    end
    prc_nme = mmil_readtext(['/home/ekaestne/PROJECTS/SCRIPTS/streamint' '/' 'parc.aparc' fcfg.prc_nme '.annot']);
    prc_nme = prc_nme(:,2)';
    sve_sbj_nme = cfg.sbj_fsr(:,1);
    save([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'vox' '_' 'aparc' mmil_spec_char(fcfg.prc_nme,{'.'}) '.mat'],'sve_sbj_nme','prc_nme','vox_dta');
    clear vox_dta
end
toc

%% fMRI
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
fcfg.fmr_dir = cfg.bld_dir;
fcfg.fmr_fsr_dir = cfg.fmr_fsr_dir;

for iP = 1:numel(cfg.lbl_nme)
    fcfg.prc_nme = cfg.lbl_nme{iP};
    for iFS = 1:numel( cfg.fmr_stt)
        fcfg.fmr_stt = cfg.fmr_stt{iFS};
        for iFN = 1:numel(cfg.fmr_nme)
            fcfg.fmr_nme = cfg.fmr_nme{iFN};
            for iFT = 1:numel(cfg.fmr_typ)
                fcfg.fmr_typ = cfg.fmr_typ{iFT};
                
                for iS = 1:size(cfg.sbj_fsr,1)
                    if any(strcmpi(fmr_sbj_dir(:,1),cfg.sbj_fsr{iS,1}))
                        
                    fcfg.sbj_fmr_dir = fmr_sbj_dir{strcmpi(fmr_sbj_dir(:,1),cfg.sbj_fsr{iS,1}),2};
                    
                    fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
                    fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
                    [fmr_dta_hld , tot_lbl]= ejk_extract_fmri_roi_voxels(fcfg);
                    if ~isempty(fmr_dta_hld)
                        if size(prc_nme,2) == size(fmr_dta_hld,1)
                            if ~isempty(fmr_dta_hld)
                                fmr_dta(iS,:) = cell2mat(fmr_dta_hld(:,2))'; prc_nme = fmr_dta_hld(:,1)';
                            else
                                fmr_dta(iS,:) = nan(1,size(tot_lbl,1)-2); fmr_dta(iS,:) = nan(1,size(fmr_dta,2));
                            end
                        else
                            for iRW = 1:size(prc_nme,2)
                                if any(strcmpi(fmr_dta_hld(:,1),prc_nme{iRW}))
                                    fmr_dta(iS,iRW) = fmr_dta_hld{strcmpi(fmr_dta_hld(:,1),prc_nme{iRW}),2};
                                else
                                    fmr_dta(iS,iRW) = 0;
                                end
                            end
                        end
                    end
                    else
                        fmr_dta(iS,:) = nan(1,size(fmr_dta,2));
                    end
                end
                sve_sbj_nme = cfg.sbj_fsr(:,1);
                save([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'fMRI' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '_' mmil_spec_char(fcfg.fmr_nme,{'.' '-' ' '}) '_' mmil_spec_char(fcfg.fmr_typ,{'.' '-' ' '}) '.mat'],'sve_sbj_nme','prc_nme','fmr_dta');
                cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'fMRI' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '_' mmil_spec_char(fcfg.fmr_nme,{'.' '-' ' '}) '_' mmil_spec_char(fcfg.fmr_typ,{'.' '-' ' '}) '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(fmr_dta,1))] [prc_nme ; num2cell(fmr_dta)] ])                
                clear fmr_dta
            end
        end
    end
end
toc

%% Volumes
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.ovr_dir = cfg.mri_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iP = 1:numel(cfg.lbl_nme)
    fcfg.prc_nme = cfg.lbl_nme{iP};
    for iS = 1:size(cfg.sbj_fsr,1)
        fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
        [ vol_dta_hld , tot_lbl] = ejk_extract_volumes(fcfg);
        if ~isempty(vol_dta_hld); vol_dta(iS,:) = vol_dta_hld; else vol_dta(iS,:) = nan(1,size(vol_dta,2)); end
    end
    
    sve_sbj_nme = cfg.sbj_fsr(:,1);
    cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Volumes' '_' 'aparc' '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(vol_dta,1))] [tot_lbl ; num2cell(vol_dta)] ]);
    clear vol_dta
    
end
toc

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

%% Fibers
fib_typ = { 'FA' 'MD' };
scl_fct = [ 1    1000 ];

tic;
for iFT = 1:numel(fib_typ)
    for iS = 1:size(cfg.sbj_fsr,1)
        
        fcfg = [];
        
        fcfg.ovr_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/proc_dti';
        
        fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
        fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
        
        fcfg.atl_dir = 'AtlasTrack';
        fcfg.map_dir = 'fiber_maps';
        fcfg.fib_typ = fib_typ{iFT};
        fcfg.scl_fct = scl_fct(iFT);
        
        fcfg.thr_prb = 0.08;
        
        fcfg.fib_cde = [101 102 103 104 105 106 107 108 109 110 115 116 117 118 119 120 121 122 123 133 134 135 136 137 138 141 142 143 144 145 146 147 148 149 150 1014 1024 2000 2001 2002 2003 2004];
        
        fcfg.min_val = 1e-06;
        
        [ fib_dta_hld , tot_lbl] = ejk_extract_fibers(fcfg);
        if ~isempty(fib_dta_hld); fib_dta(iS,:) = fib_dta_hld; else fib_dta(iS,:) = nan(1,size(fib_dta,2)); end
        
    end
    
    sve_sbj_nme = cfg.sbj_fsr(:,1);
    cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'Fibers' '_' 'aparc' '_' fib_typ{iFT} '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(fib_dta,1))] [tot_lbl ; num2cell(fib_dta)] ]);
    
end
toc;

%% Superficial White Matter
fib_typ = { 'FA' 'MD' }; % 'FA' 'MD'
scl_fct = [ 1    1000 ];
prc_nme = { '' }; % '' '.a2009s'

tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.dti_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iPR = 1:numel(cfg.lbl_nme)
    fcfg.prc_nme = prc_nme{iPR};
    
    for iFT = 1:numel(fib_typ)
        fcfg.fib_typ = fib_typ{iFT};
        fcfg.scl_fct = scl_fct(iFT);
        
        for iS = 1:size(cfg.sbj_fsr,1)
            fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
            fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
            [ wmp_dta_hld , tot_lbl ]= ejk_extract_wmparc(fcfg);
            if ~isempty(wmp_dta_hld); wmp_dta(iS,:) = wmp_dta_hld; else wmp_dta(iS,:) = nan(1,size(tot_lbl,2)); end
        end
        
        sve_sbj_nme = cfg.sbj_fsr(:,1);
        cell2csv([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'WMParc' '_' 'aparc' '_' mmil_spec_char(fcfg.prc_nme,{'.'}) '_' fib_typ{iFT} '.csv'],[ ['SubjId' ; sve_sbj_nme(1:size(wmp_dta,1))] [tot_lbl ; num2cell(wmp_dta)] ]);
        clear wmp_dta
        
    end
end
toc

%% Calculate Alicia
fcfg = [];
fcfg.prj_dir = cfg.prj_dir;
ejk_extract_alicia(fcfg)

%% Calculate LI

%% NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONNECTOME - NORM
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.ovr_dir = cfg.cnn_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iS = 1:size(cfg.sbj_fsr,1)
    fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
    fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
    [ cnn_dta_hld , row_lbl_hld , col_lbl_hld ]= ejk_extract_connectome(fcfg);
    if ~isempty(cnn_dta_hld); cnn_dta(iS,:,:) = cnn_dta_hld; row_lbl = row_lbl_hld; col_lbl = col_lbl_hld; else cnn_dta(iS,:,:) = nan(1,size(cnn_dta,2),size(cnn_dta,3)); end
end

sve_sbj_nme = cfg.sbj_fsr(:,1);

save([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'DTI' '_' 'Connectome' '_' 'norm' '.mat'],'sve_sbj_nme','cnn_dta','row_lbl','col_lbl');
clear cnn_dta row_lbl col_lbl

toc

%% GRAPH THEORY
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.ovr_dir = cfg.cnn_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

for iS = 1:size(cfg.sbj_fsr,1)
    fcfg.sbj_fsr_dir = cfg.sbj_fsr{iS,2};
    fcfg.sbj_nme     = cfg.sbj_fsr{iS,1};
    [ cnn_dta_hld , row_lbl_hld , col_lbl_hld ]= ejk_extract_connectome(fcfg);
    if ~isempty(cnn_dta_hld); cnn_dta(iS,:,:) = cnn_dta_hld; row_lbl = row_lbl_hld; col_lbl = col_lbl_hld; else cnn_dta(iS,:,:) = nan(1,size(cnn_dta,2),size(cnn_dta,3)); end
end

sve_sbj_nme = cfg.sbj_fsr(:,1);

save([cfg.prj_dir '/' 'DATA' '/' 'ROIHOLD' '/' 'DTI' '_' 'Connectome' '_' 'norm' '.mat'],'sve_sbj_nme','cnn_dta','row_lbl','col_lbl');
clear cnn_dta row_lbl col_lbl

toc

%% fMRI CONNECTIVITY

%% SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MRI Thickness Surf Data
tic;
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
    save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% DTI FA Surf Data
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.dti_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

fcfg.mes_typ = 'wmparc_fa'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 313; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

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
    save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% DTI MD Surf Data
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.dti_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

fcfg.mes_typ = 'wmparc_md'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 313; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

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
    save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% rsfMRI Variance Surf Data
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.fmr_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

fcfg.mes_typ = 'rsfMRI_variance'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 256; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

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
    save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

%% rsfMRI Surf Data
tic;
fcfg = [];

fcfg.prj_dir = cfg.prj_dir;
fcfg.prc_dir = cfg.fmr_fsr_dir; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf

fcfg.mes_typ = 'rsfMRI'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
fcfg.smt_stp  = 256; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;

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
    save([cfg.prj_dir '/' 'DATA' '/' 'SRFHOLD' '/' 'surf' '_' fcfg.mes_typ '_' fcfg.hms 's.mat'],'srf_dta_sbj','srf_dta');
    clear srf_dta
end
toc

end






