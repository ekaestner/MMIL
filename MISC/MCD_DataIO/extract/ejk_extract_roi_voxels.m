function vox_dta  = ejk_extract_roi_voxels(cfg)

%% Parcellation Table
prc_nme = mmil_readtext(['/home/ekaestne/PROJECTS/SCRIPTS/streamint' '/' 'parc.aparc' cfg.prc_nme '.annot']);
vox_dta = nan(1,size(prc_nme,1));

try
%% Load
sbj_dir_lst = dir(sprintf('%s/FSURF*',cfg.prc_dir));
n = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); n = [n{:}];
if isempty(n); fprintf('%s: WARNING: no Recon container found for subject %s\n',mfilename,cfg.sbj_fsr_dir); end
if ~isempty(n); con_dir = n{1}; con_pth = sprintf('%s/%s',cfg.prc_dir,con_dir); end

anl_dir = 'mri';

fle_nme = sprintf('%s/%s/aparc%s+aseg.mgz',con_pth,anl_dir,cfg.prc_nme);

vol_dta = fs_load_mgh(fle_nme);

%% Get data
for iRW = 1:size(prc_nme,1)
    vox_dta(1,iRW) = sum(vol_dta(:)==prc_nme{iRW,1});
end
catch; end

end