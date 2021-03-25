function mb_recon(fname_pfile,outdir,outstem,forceflag)
%function mb_recon(fname_pfile,outdir,outstem,forceflag)
%
% Purpose: multi-band recon
%
% Required Input:
%   fname_pfile: P-file containing raw k-space multi-band data
%
% Optional Input:
%   outdir: output directory
%     {default = [pwd '/mb_recon']}
%   outstem: output file stem
%     {default = 'mb_recon'}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  02/11/17 by Don Hagler
% Last Mod: 02/11/17 by Don Hagler
%
% Copied from call_mux_epi_main
% Prev Mod:  02/03/16 by ?
% Last Mod:  02/11/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist(fname_pfile,'file')
  error('file %s not found',fname_pfile);
end;
if ~exist('outdir','var') || isempty(outdir)
  outdir = [pwd '/mb_recon'];
end;
if ~exist('outstem','var') || isempty(outstem)
  outstem = 'mb_recon';
end;
if ~exist('forceflag','var') || isempty(forceflag)
  forceflag = 0;
end;
mmil_mkdir(outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mat = sprintf('%s/%s_header.mat',outdir,outstem);
if ~exist(fname_mat,'file') || forceflag
  dataFileHeader = read_MR_headers(fname_pfile,'all','raw');
  nslices = dataFileHeader.rdb_hdr.nslices/dataFileHeader.rdb_hdr.reps;
  ncoil = (dataFileHeader.rdb_hdr.dab(2)-dataFileHeader.rdb_hdr.dab(1))+1;
  % coil compression: number of virtual coils
  n_vcoils = ncoil/2; 
  save(fname_mat,'dataFileHeader','nslices','n_vcoils');
else
  load(fname_mat);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mat = sprintf('%s/%s_maxim.mat',outdir,outstem);
if ~exist(fname_mat,'file') || forceflag
  % check for pepolar
  fprintf('%s: getting mux epi parameters...\n',mfilename);
  p = mux_epi_params(fname_pfile,[],[],n_vcoils,1,'1Dgrappa',1,1,1,[]);
  if bitand(p.dacq_ctrl, 4)
    % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
    nt_to_recon = []; % safer to be = [1,2]?
  else
    nt_to_recon = 1;
  end
  % initial reconstruction to calculate maximum image intensity
  fprintf('%s: performing initial reconstruction...\n',mfilename);
  d = mb_mux_epi(fname_pfile,[],[],[],nt_to_recon,n_vcoils,1,'1Dgrappa',1,1,1);
  T2 = d(:,:,:,2);
  maxim = real(max(T2(:)));
  clear T2 d
  % save maxim for later use for scaling dicoms
  save(fname_mat,'maxim');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% repeat reconstruction, one slice at a time
%warning off;
for sl = 1:nslices
  outfile = sprintf('%s/%s_s%05d.mat',outdir,outstem,sl);
  if ~exist(outfile,'file') || forceflag
    fprintf('%s: processing slice #%d/%d...\n',mfilename,sl,nslices);
    mb_mux_epi(fname_pfile,outfile,[],sl,[],n_vcoils,1,'1Dgrappa',1,1,1);
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

