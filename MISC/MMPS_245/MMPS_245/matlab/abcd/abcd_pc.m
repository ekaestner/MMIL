function abcd_pc(indir,outdir)
%function abcd_pc(indir,outdir)
%
% Created:  08/18/16 by Jose Antolin
% Last Mod: 09/29/16 by Don Hagler
%

% NOTE: based on QC_PC_local, last mod 08/18/16, created by Jose Antolin

if ~mmil_check_nargs(nargin,2), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create output folder
mmil_mkdir(outdir);

% get all folders
se_list = dir(indir);
se_list = se_list([se_list.isdir] & ~strncmpi('.', {se_list.name}, 1)); 

% select one series
for i=1:length(se_list)
  se_indir = fullfile(indir,se_list(i).name);
  se_outdir = fullfile(outdir,se_list(i).name);
  fprintf('%s: checking series %d in %s...',mfilename,i,se_indir);
  tic;
  proc_series(se_indir,se_outdir);
  toc;
end;

fname_check = sprintf('%s/.finished',outdir);
fid = fopen(fname_check,'wt');
if fid<0
  error('failed to open file %s for writing',fname_check);
end;
fprintf(fid,indir);
fclose(fid);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function proc_series(se_indir,se_outdir)
  mmil_mkdir(se_outdir);
  fid = fopen(fullfile(se_outdir,'PC.log'),'w');
  if fid==-1
    error('failed to write in output directory %s',se_outdir);
  end;
  if ~exist(se_indir)
    fprintf(fid,'ERROR: input directory %s does not exist\n',se_indir);
    fprintf('%s: WARNING: input directory %s does not exist\n',mfilename,se_indir);
    return;
  end;
  fnames = dir(se_indir);
  fnames = setdiff({fnames.name},{'.','..'});
  if isempty(fnames)
    fprintf(fid,'ERROR: input directory %s is empty\n',se_indir);
    fprintf('%s: WARNING: input directory %s is empty\n',mfilename,se_indir);
    return
  end
  
  n_files = length(fnames);
  for i=1:n_files
    fname = fullfile(se_indir,fnames{i});
    if ~file_is_dicom(fname)
      fprintf(fid,'ERROR: input directory must contain only dicom files\n');
      fprintf('%s: WARNING: input directory must contain only dicom files\n');
      return;
    end
    metadata = dicominfo(fname);
    if ~isempty(mmil_getfield(metadata,'AcquisitionMatrix'));
      break;
    end  
  end

  %%% Add the json file for phase encoding direcion in SIEMENS FM %%
  %json_file = dir(fullfile(se_indir, '/../*.json'));
  %se_json = fullfile(se_indir, '/..', json_file(1).name);
  %json_info = loadjson(se_json);
  %%% ===================================================%
  
  if strfind(metadata.Manufacturer,'GE')
    abcd_pc_GE(metadata,n_files,se_outdir);
  elseif strfind(metadata.Manufacturer,'SIEMENS')
    csadata = spm_dicom_headers(fname);
    abcd_pc_Siemens(metadata,n_files,se_outdir,csadata);
  elseif strfind(metadata.Manufacturer, 'Philips')
    abcd_pc_Philips(metadata,n_files,se_outdir);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isdicom=file_is_dicom(filename)
  isdicom=false;
  try
    isdicom = mmil_isdicomfile(filename);
  catch ME
  end
return;

