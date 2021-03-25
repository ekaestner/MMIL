function mmil_aseg_analysis(fname,fspath,varargin)
%function mmil_aseg_analysis(fname,fspath,[options])
%
% Required Parameters:
%   fname: full path name of file containing image data
%   fspath: full path of freesurfer recon
%
% Optional Parameters:
%   'outdir': output directory
%     if empty, will attempt to use path of fname
%     {default = []}
%   'outstem': output file stem
%     if empty, will use file stem of fname
%     {default = []}
%   'fname_out': relative or full path name for output mat file
%     if supplied, outstem will be ignored
%     if full path, outdir will be ignored
%     if not, fname_out will be constructed from outstem plus other info
%     {default = []}
%   'csv_flag': [0|1] summarize results in csv file
%     name will have same stem as fname_out, except with .csv extension
%     {default = 1}
%   'fname_aseg': segmentation volume to specify ROIs
%     if empty, use fspath/mri/aseg.mgz
%       or aparc+aseg.mgz depending on aseg_aparc_flag
%     {default = []}
%   'aseg_aparc_flag': [0|1|2] whether to use cortical parcellation volume ROIs
%     0: aseg only
%     1: aparc only   NOTE: for best results with aparc, use erode_flag=0
%     2: aparc+aseg
%     {default = 0}
%   'fname_vals': full path name of file containing values to be used
%     to calculate weights based on ROI dispersion values in dispvec
%     {default = []}
%   'dispvec': vector of ROI dispersion values (MAD estimates)
%     required if fname_vals supplied
%     {default = []}
%   'disp_roicodes': vector of ROI codes corresponding to dispvec elements
%     required if fname_vals supplied
%     {default = []}
%   'disp_scalefact': scaling factor applied to fname_vals when
%     comparing values to dispersion values
%     {default = 1000}
%   'disp_suffix':  suffix attached to output file names
%     using dispersion weighting
%     {default = 'dwtd'}
%   'dispfact': multiplicative factor applied to dispersion values
%     {default = 4}
%   'erode_flag': [0|1] whether to "erode" ROIs by smoothing and thresholding
%     {default = 1}
%   'erode_nvoxels': number of voxels to erode (integer)
%     {default = 1}
%   'scalefact': scaling factor applied to extracted values
%     {default = 1}
%   'minval': minimum value; if NaN, include all voxels
%     {default = 1e-6}
%   'M_reg': 4x4 transformation matrix between fname and fspath/mri/orig.mgz
%     should be like M_T1_to_fname
%     If supplied, aseg will be resampled to space of fname
%     {default = []}
%   'res_outfix': output string added to resampled aseg fname
%     Ignored if M_reg is empty
%     {default = 'res'}
%   'fname_colorlut': full path of color look up table file (to get ROI names)
%     if empty, use FreeSurferColorLUT.txt in $FREESURFER_HOME
%    {default = []}
%   'verbose': [0|1] display status meassages
%     {default = 0}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  08/10/12 by Don Hagler
% Prev Mod: 06/03/12 by Don Hagler
% Last Mod: 08/15/17 by Don Hagler
%

%% todo: rename to mmil_volroi_analysis?
%% todo: accept aseg_roigroups
%% todo: accept additional mgh ROI files (flist_hippo_masks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = check_input(fname,fspath,varargin);

parms = set_output(parms);

parms = prep_aseg(parms);

extract_values(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,fspath,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    'fspath',fspath,[],...
  ...
    'outdir',[],[],...
    'outstem',[],[],...
    'fname_out',[],[],...
    'csv_flag',true,[false true],...
    'fname_aseg',[],[],...
    'aseg_aparc_flag',0,[0,1,2],...
    'fname_vals',[],[],...
    'dispvec',[],[],...
    'disp_roicodes',[],[],...
    'disp_scalefact',1000,[1e-10,1e10],...
    'disp_suffix','dwtd',[],...
    'dispfact',4,[1e-6,1e6],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'scalefact',1,[1e-10,1e10],...
    'minval',1e-6,[],...
    'M_reg',[],[],...
    'res_outfix','res',[],...
    'fname_colorlut',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden parameters
    'aseg_roilist',[1:28,40:60],[],...
    'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
    'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[],...
    'frames',1,[],...
  ...
    'code_tags',{'aseg_aparc_flag','aseg_roilist','aparc_roilist',...
                 'exclude_roilist'},[],...
    'roi_tags',{'fname_weights','fname_mask','roilist','roigroups',...
                'minval','mask_thresh','scalefact','frames',...
                'fname_colorlut','verbose'},[],...
    'erode_tags',{'nvoxels','sigma','thresh','verbose','forceflag'},[],...
    'weights_tags',{'fname_out','scalefact','dispfact','erode_flag',...
                    'erode_nvoxels','fname_aseg_erode','roicodes',...
                    'verbose','forceflag'},[],...
  });

  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;

  % check FreeSurfer recon exists if needed
  if isempty(parms.fname_aseg) || mmil_isrelative(parms.fname_aseg)
    if ~exist(parms.fspath,'dir')
      error('FreeSurfer recon dir %s not found',parms.fspath);
    end;
  end;

  % set aseg file name and check that it exists
  if isempty(parms.fname_aseg)
    if parms.aseg_aparc_flag
      parms.aseg_name = 'aparc+aseg';
    else
      parms.aseg_name = 'aseg';
    end;
    parms.fname_aseg = [parms.aseg_name '.mgz'];
  else
    [tmp,parms.aseg_name] = fileparts(parms.fname_aseg);
  end;
  if mmil_isrelative(parms.fname_aseg)
    parms.fname_aseg = [parms.fspath '/mri/' parms.fname_aseg];
  end;
  if ~exist(parms.fname_aseg,'file')
    error('aseg file %s not found\n',mfilename,parms.fname_aseg);
  end
  args = mmil_parms2args(parms,parms.code_tags);
  parms.roilist = mmil_aseg_roicodes(args{:});

  % set parameter name for fs_erode_aseg
  parms.nvoxels = parms.erode_nvoxels;

  % check weights-related parameters
  if ~isempty(parms.fname_vals)
    parms.weights_flag = 1;
    if ~exist(parms.fname_vals,'file')
      error('values for weights volume file %s not found',parms.fname_vals);
    end;
    if isempty(parms.dispvec)
      error('dispvec required if fname_vals supplied');
    end;
    if isempty(parms.disp_roicodes)
      error('disp_roicodes required if fname_vals supplied');
    end;
    if length(parms.dispvec) ~= length(parms.disp_roicodes)
      error('mismatch between dispvec and disp_roicodes');
    end;
  else
    parms.weights_flag = 0;
  end;
  
  % disable csv writing if num frames is > 1
  if length(parms.frames)>1 && parms.csv_flag
    if parms.verbose
      fprintf('%s: setting csv_flag = 0 because number of frames > 1\n',mfilename);
    end;
    parms.csv_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_output(parms)
  if isempty(parms.fname_out)  
    % set output file name
    if isempty(parms.outstem)
      [tmp,parms.outstem] = fileparts(parms.fname);
    end;
    switch parms.aseg_aparc_flag
      case {0,2}
        parms.outstem = [parms.outstem '_' parms.aseg_name];
      case 1
        parms.outstem = [parms.outstem '_aparc'];
    end;
    if parms.erode_flag
      parms.outstem = sprintf('%s_erode',parms.outstem);
      if parms.erode_nvoxels>1,
        parms.outstem = sprintf('%s%d',parms.outstem,parms.erode_nvoxels);
      end;
    end;
    if parms.weights_flag
      parms.outstem = sprintf('%s_%s%0.1f',...
        parms.outstem,parms.disp_suffix,parms.dispfact);
    end;
    parms.fname_out = sprintf('%s_roi_data.mat',parms.outstem);
  end;
  [tmp,tmp,ext] = fileparts(parms.fname_out);
  if ~strcmp(ext,'.mat')
    error('fname_out must have .mat extension');
  end;

  % set output directory and create if necessary
  if ~mmil_isrelative(parms.fname_out)
    parms.outdir = fileparts(parms.fname_out);
  else
    % set output directory
    if isempty(parms.outdir)
      parms.outdir = fileparts(parms.fname);
    end;
    parms.fname_out = [parms.outdir '/' parms.fname_out];
  end;
  mmil_mkdir(parms.outdir);

  parms.fname_csv = regexprep(parms.fname_out,'\.mat','.csv');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_aseg(parms)
  aseg_outstem = parms.aseg_name;
  % resample aseg to fname space
  fname_aseg = parms.fname_aseg;
  if ~isempty(parms.M_reg)
    aseg_outstem = [aseg_outstem '_' parms.res_outfix];
    fname_aseg = [parms.outdir '/' aseg_outstem '.mgz'];  
    if ~exist(fname_aseg,'file') || parms.forceflag
      [M,volsz] = fs_read_header(parms.fname);
      [vol_aseg,M_aseg] = fs_load_mgh(parms.fname_aseg);
      [vol_aseg,M_aseg] = mmil_resample_vol(vol_aseg,M_aseg,...
        'nvox_ref',volsz,'M_ref',M,'interpm',0,'M_reg',inv(parms.M_reg));
      fs_save_mgh(vol_aseg,fname_aseg,M_aseg);
    end
    parms.fname_aseg = fname_aseg;
  end;
  % erode aseg
  if parms.erode_flag
    aseg_outstem = [aseg_outstem '_erode'];
    if parms.erode_nvoxels>1
      aseg_outstem = sprintf('%s%d',aseg_outstem,parms.erode_nvoxels);
    end;
    fname_aseg_erode = [parms.outdir '/' aseg_outstem '.mgz'];  
    args = mmil_parms2args(parms,parms.erode_tags);
    fs_erode_aseg(parms.fname_aseg,fname_aseg_erode,args{:});
    parms.fname_aseg = fname_aseg_erode;
  end;
  % calculate weights
  if parms.weights_flag
    tmp_parms = parms;
    tmp_parms.fname_out = sprintf('%s/%s_weights_df%0.1f.mgz',...
      parms.outdir,aseg_outstem,parms.dispfact);
    if parms.erode_flag
      tmp_parms.fname_aseg_erode = fname_aseg_erode;
    end;
    tmp_parms.roicodes = parms.disp_roicodes;
    tmp_parms.scalefact = parms.disp_scalefact;
    args = mmil_parms2args(tmp_parms,parms.weights_tags);
    mmil_aseg_weights(fname_aseg,parms.fname_vals,parms.dispvec,args{:});
    parms.fname_weights = tmp_parms.fname_out;
    parms.fname_aseg = fname_aseg;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extract_values(parms)
  roi_data = [];
  if ~exist(parms.fname_out,'file') || parms.forceflag
    args = mmil_parms2args(parms,parms.roi_tags);
    roi_data = mmil_aseg_roi(parms.fname,parms.fname_aseg,args{:});    
    if isempty(roi_data), error('aseg ROI analysis failed'); end
    save(parms.fname_out,'roi_data');
  end
  if parms.csv_flag, write_csv(parms); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_csv(parms)
  roi_data = [];
  if ~exist(parms.fname_csv,'file') || parms.forceflag
    load(parms.fname_out);
    fid = fopen(parms.fname_csv,'wt');
    if fid==-1, error('failed to open %s for writing',parms.fname_csv); end;
    fprintf(fid,'"ROI","mean","median","stdv","nvals","nvals valid"\n');
    nrois = length(roi_data);
    for i=1:nrois
      fprintf(fid,'"%s",',roi_data(i).roiname);
      fprintf(fid,'%0.6f,%0.6f,%0.6f,%d,%d\n',...
        roi_data(i).avg,roi_data(i).median,roi_data(i).stdv,...
        roi_data(i).nvals,roi_data(i).nvals_valid);
    end
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

