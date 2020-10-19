function mb_write_dicoms(indir,fname_dicom,outdir,instem)
%function mb_write_dicoms(indir,fname_dicom,[outdir],[instem])
%
% Purpose: write dicoms after multi-band recon
%
% Required Input:
%   indir: input directory containing mb_recon mat files
%     outdir produced by mb_recon
%   fname_dicom: single dicom file containing dicom metadata
%
% Optional Input:
%   outdir: output directory
%     {default = [pwd '/mb_dicoms'];
%   instem: input file stem for mb_recon mat files
%     {default = 'mb_recon'}
%
% Created:  02/11/17 by Don Hagler
% Last Mod: 02/14/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist(indir,'dir')
  error('directory %s not found',indir);
end;
if ~exist(fname_dicom,'file')
  error('file %s not found',fname_dicom);
end;
if ~exist('outdir','var') || isempty(outdir)
  outdir = [pwd '/mb_dicoms'];
end;
if ~exist('instem','var') || isempty(instem)
  instem = 'mb_recon';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find file containing header
fname_header = sprintf('%s/%s_header.mat',indir,instem);
if ~exist(fname_header,'file')
  error('file %s not found',fname_header);
end;

% find file containing scaling factor
fname_maxim = sprintf('%s/%s_maxim.mat',indir,instem);
if ~exist(fname_maxim,'file')
  error('file %s not found',fname_maxim);
end;

% find mat files in indir
flist = dir(sprintf('%s/%s_s*.mat',indir,instem));
if isempty(flist)
  error('indir %s contains no mat files with instem %s',indir,instem);
end;
fnames = {flist.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(outdir);

% load header
tmp = load(fname_header);
dataFileHeader = tmp.dataFileHeader;

% load maxim
tmp = load(fname_maxim);
maxim = tmp.maxim;

if maxim<0
  error('maximum image intensity is %0.1f (should be > 0)',maxim);
end;

% get new series UID for dicom write
%% todo: use original SeUID?
new_se_uid = get_new_se_uid(fname_dicom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over mat files, creating dicoms for each slice
warning off;
for i=1:length(fnames)
  fname_mat = sprintf('%s/%s',indir,fnames{i});
  data = [];
  data = load(fname_mat);
  slice_demux_coord = calc_demux_coord(dataFileHeader,data.p);
  for t = 1:size(data.d,4)
    for sl = 1:size(data.d,3)
      d_(:,:,sl,t) = imrotate(data.d(:,:,sl,t),...
                              dataFileHeader.rdb_hdr.rotation*90);
    end
  end
  t_array = 1:size(d_(:,:,:,2:end),4);
  num_t2 = 1;
  % write dicoms
  for sl_ind = 1:data.p.mux
    write_dicoms(real(squeeze(d_(:,:,sl_ind,2:end))),...
                 data.sl_loc(sl_ind),...
                 slice_demux_coord(data.sl_loc(sl_ind)),...
                 dataFileHeader,...
                 t_array,...
                 num_t2,...
                 fname_dicom,...
                 maxim,...
                 new_se_uid,...
                 outdir);
  end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slice_demux_coord = calc_demux_coord(dataFileHeader,p)
  % MJM: Information needed for DICOM writing
  slice_thk   = dataFileHeader.image.slthick;
  if mod(p.mux,2) == 1 || mod(p.num_slices,2) == 0, % Odd MB factor or even number of slices
    slice_add = (p.num_slices*p.mux-p.num_slices)/2;
    if p.descend_acq,
      slice_demux_coord = horzcat(linspace(p.start_loc+(slice_thk*slice_add),p.start_loc+slice_thk,slice_add),...
                          linspace(p.start_loc,p.end_loc,p.num_slices),...
                          linspace(p.end_loc-slice_thk,p.end_loc-(slice_add*slice_thk),slice_add));
    else
      slice_demux_coord = horzcat(linspace(p.start_loc-(slice_thk*slice_add),p.start_loc-slice_thk,slice_add),...
                          linspace(p.start_loc,p.end_loc,p.num_slices),...
                          linspace(p.end_loc+slice_thk,p.end_loc+(slice_add*slice_thk),slice_add));
    end
  else  % Even MB facor with odd number of slices
    slice_add = (p.mux*p.num_slices - 2*p.num_slices)/2;
    if p.descend_acq,
      p.start_loc = p.start_loc + p.num_slices/2*slice_thk;
      p.end_loc = p.end_loc - p.num_slices/2*slice_thk;
      slice_demux_coord = horzcat(linspace(p.start_loc+(slice_thk*slice_add),p.start_loc+slice_thk,slice_add),...
                          linspace(p.start_loc,p.end_loc,2*p.num_slices),...
                          linspace(p.end_loc-slice_thk,p.end_loc-(slice_add*slice_thk),slice_add));
    else
      p.start_loc = p.start_loc - p.num_slices/2*slice_thk;
      p.end_loc = p.end_loc + p.num_slices/2*slice_thk;
      slice_demux_coord = horzcat(linspace(p.start_loc-(slice_thk*slice_add),p.start_loc-slice_thk,slice_add),...
                          linspace(p.start_loc,p.end_loc,2*p.num_slices),...
                          linspace(p.end_loc+slice_thk,p.end_loc+(slice_add*slice_thk),slice_add));
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

