%% Testing
%
% cfg = [];
%
% cfg.prj_dir = '/home/ekaestne/PROJECTS/';
% cfg.prc_dir = '/space/syn09/1/data/MMILDB/MCD_RSI/fsurf'; % /space/syn09/1/data/MMILDB/MCD_RSI/proc_dti ; /space/syn09/1/data/MMILDB/MCD_RSI/proc_bold ; /space/syn09/1/data/MMILDB/MCD_RSI/fsurf
%
% cfg.sbj_fsr_dir = 'epd090_fmri_170515';
%
% cfg.hms     = 'lh'; %{'lh' 'rh'};
%
% cfg.mes_typ = 'aMRI_thickness'; % 'rsfMRI_variance' 'aMRI_thickness'; wmparc_fa; wmparc_md
% cfg.smt_stp  = 313; % cfg.sph_smt = 313; % rsfMRI = 256 ; aMRI = 313;
%
% cfg.anl_dir = '';

%
%
%

function srf_dta = ejk_extract_vertices(cfg)

%% Setup
switch cfg.mes_typ
    case 'rsfMRI'
        sbj_dir_lst = dir(sprintf('%s/BOLDPROC*',cfg.prc_dir));
        n = regexp({sbj_dir_lst.name},['BOLDPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}]; % '_\d{8}.+_1$'
        if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
        if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end
    case 'rsfMRI_variance'
        sbj_dir_lst = dir(sprintf('%s/BOLDPROC*',cfg.prc_dir));
        n = regexp({sbj_dir_lst.name},['BOLDPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];
        if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
        if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end
    case 'aMRI_thickness'
        sbj_dir_lst = dir(sprintf('%s/FSURF*',cfg.prc_dir));
        n = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];
        if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
        if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end
    case 'wmparc_fa'
        sbj_dir_lst = dir(sprintf('%s/DTIPROC*',cfg.prc_dir));
        n = regexp({sbj_dir_lst.name},['DTIPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];
        if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
        if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end
    case 'wmparc_md'
        sbj_dir_lst = dir(sprintf('%s/DTIPROC*',cfg.prc_dir));
        n = regexp({sbj_dir_lst.name},['DTIPROC_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];
        if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
        if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end
end

if ~isfield(cfg,'anl_dir') || strcmpi(cfg.anl_dir,'')
    switch cfg.mes_typ
        case 'rsfMRI';          anl_dir = 'rsBOLD_analysis';
        case 'rsfMRI_variance'; anl_dir = 'rsBOLD_analysis';
        case 'aMRI_thickness';  anl_dir = 'analysis';
        case 'wmparc_fa';       anl_dir = 'DTanalysis';
        case 'wmparc_md';       anl_dir = 'DTanalysis';
    end
else
    anl_dir = cfg.anl_dir;
end

if isvar('con_pth')
    switch cfg.mes_typ
        case 'rsfMRI',          fle_nme = sprintf('%s/%s/rsBOLD-sm%d-ico7-%s.mgz',        con_pth,anl_dir,cfg.smt_stp,cfg.hms);
        case 'rsfMRI_variance', fle_nme = sprintf('%s/%s/rsBOLD-sm%d-ico7-var-%s.mgz',    con_pth,anl_dir,cfg.smt_stp,cfg.hms);
        case 'aMRI_thickness',  fle_nme = sprintf('%s/%s/thickness-sphere-sm%d-%s.mgz',   con_pth,anl_dir,cfg.smt_stp,cfg.hms);
        case 'wmparc_fa',       fle_nme = sprintf('%s/%s/FA_pdist-1.0-sphere-sm%d-%s.mgz',con_pth,anl_dir,cfg.smt_stp,cfg.hms);
        case 'wmparc_md',       fle_nme = sprintf('%s/%s/MD_pdist-1.0-sphere-sm%d-%s.mgz',con_pth,anl_dir,cfg.smt_stp,cfg.hms);
    end
end

if isvar('fle_nme')
    
    if ~exist(fle_nme,'file')
        fprintf('%s: WARNING: %s not found\n',mfilename,fle_nme);
        not_her = 1;
    else
        not_her = 0;
    end
    
    if not_her
        srf_dta = [];
    else
        
        srf_dta = mmil_rowvec(fs_load_mgh(fle_nme));
        if isempty(srf_dta)
            print('%s: WARNING: unable to correctly read %s\n',mfilename,fle_nme);
            not_her = 1;
        end
        
        fprintf('%s: %s: loading data...\n',mfilename,con_pth);
        
        if not_her
            srf_dta = [];
        else
            srf_dta = mmil_rowvec(fs_load_mgh(fle_nme));
        end
    end
    
else
    srf_dta = [];
end

end
