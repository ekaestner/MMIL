function [vol,M] = mmil_read_dicom_4dvol(fnames,nrow,ncol,nslices,nreps,...
  mosaic_flag);
%function [vol,M] = mmil_read_dicom_4dvol(fnames,nrow,ncol,nslices,nreps,...
%  [mosaic_flag]);
%
% Purpose: read 4d volume (3d x nreps) from dicoms
%
% Required Input:
%   fnames: cell array of dicom file names (must be pre-sorted)
%           order should be with slices for inner loop, reps for outer loop
%   nrow: number of rows
%   ncol: number of columns
%   nslices: number of slices
%   nreps: number of repetitions
%
% Optional Input:
%   mosaic_flag: [0|1] convert mosaic images to volumes
%     {default = 0}
%
% Created:  10/31/06 Don Hagler
% Last Mod: 04/16/13 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,5), return; end;
vol = []; M=eye(4);

if ~exist('mosaic_flag','var') || isempty(mosaic_flag), mosaic_flag=0; end;

if ~iscell(fnames), fnames = {fnames}; end;
nfiles = length(fnames);

nii_flag = 1; %1 = convert to nifti using dcm2nii and read vol. 0 = read dicoms directly using dicomread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mosaic_flag
  if nfiles ~= nreps
    %% todo: error?
    fprintf('%s: ERROR: nfiles (%d) does not match nreps (%d)\n',...
      mfilename,nfiles,nreps);
    return;
  end;
  M = mmil_read_dicom_M(fnames(1));
  if isempty(M)
    vol=[];
    return;
  end;

  dcminfo = dicominfo(char(fnames{1}));

  if nii_flag
    vol = mmil_dcm2nii2vol(fnames);
    mnr = double(dcminfo.Rows);
    mnc = double(dcminfo.Columns);
  else
    vol = zeros(nrow,ncol,nslices,nreps, 'single');
    for f=1:nfiles
      fname = fnames{f};
      if ~exist(fname,'file')
        %% todo: error?
        fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
        vol=[];
        return;
      end;
      im = dicomread(fname)';
      if(isempty(im))
        %% todo: error?
        fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
        vol=[];
        return;
      end
      try
        [mnr,mnc] = size(im);
        mr_nimgs = mnr / nrow;
        mc_nimgs = mnc / ncol;
        mos_ns = mr_nimgs * mc_nimgs;
        z = 1;
        for i=1:mc_nimgs
          c1 = 1 + (i-1)*ncol;
          c2 = c1 + ncol - 1;
          for j=1:mr_nimgs
            r1 = 1 + (j-1)*nrow;
            r2 = r1 + nrow - 1;
            vol(:,:,z,f) = im(r1:r2,c1:c2);
            if z>=nslices
              break;
            else
              z = z + 1;
            end;
          end;
        end;
      catch
        fprintf('%s: ERROR: problem with pixel data from %s: %s\n',...
          mfilename,fname,lasterr);
        vol=[];
        return;
      end;
    end;
  end
  % adjust for incorrect mosaic slice position information
  % see: http://nipy.sourceforge.net/nibabel/dicom/dicom_mosaic.html
  a = [mnr-nrow;mnc-ncol]/2;
  M(1:3,4) = M(1:3,4) + M(1:3,1:2)*a;
  % make adjustments for Siemens scans
  %dcminfo = dicominfo(char(fnames{1}));
  if strfind(upper(dcminfo.Manufacturer),'SIEMENS')
    % check for reversed slice direction
    if isfield(dcminfo,'Private_0029_1020') % Siemens CSA
      [mdhStruct,txt] = mmil_read_mdhStruct(dcminfo);
      orient = fs_read_orient([],M);
      switch orient(3)
        case {'R','L'}
          tag = 'sSliceArrayducImageNumbSag';
        case {'A','P'}
          tag = 'sSliceArrayducImageNumbHor';
        case {'S','I'}
          tag = 'sSliceArrayducImageNumbTra';
      end;
      if isfield(mdhStruct,tag)
        slice_rev = get_mdh_tagval(mdhStruct,tag);
        if slice_rev~=0
          fprintf('%s: reversing slice order for Siemens DTI...\n',mfilename);
          if ~nii_flag, vol = vol(:,:,end:-1:1,:); end;
          % adjust offset
          M(1:3,4) = M(1:3,4) - M(1:3,1:3)*[0 0 nslices-1]';
        end;
      end;
      % see: http://www.nmr.mgh.harvard.edu/~greve/dicom-unpack
      % Siemens may reverse the slice order in order to make the images more readable by
      %   radiologists. This reversal is NOT accompanied by a change the slice direction cosine.
      %   Rather, it is indicated by the presence of sSliceArray.ucImageNumbAAA (any of the three)
      %   The FreeSurfer software reverses the slices upon read-in rather than chaning
      %   the direction cosine.
    end;
  end;
else
  if nfiles > nslices*nreps
    fprintf('%s: WARNING: nfiles (%d) is greater than nslices*nreps (%dx%d = %d) -- will use first %d reps only\n',...
      mfilename,nfiles,nslices,nreps,nslices*nreps,nreps);
  elseif nfiles ~= nslices*nreps
    %% todo: error?
    fprintf('%s: ERROR: nfiles (%d) does not match nslices*nreps (%dx%d)\n',...
      mfilename,nfiles,nslices,nreps);
    return;
  end;
  M = mmil_read_dicom_M({fnames{1:nslices}});
  if isempty(M)
    vol=[];
    return;
  end;
  if nii_flag
    vol = mmil_dcm2nii2vol(fnames);
  else
%  dcminfo = dicominfo(char(fnames{1}));
%  nrow = dcminfo.Rows;  %% should trust input?
%  ncol = dcminfo.Columns;
    vol = zeros(nrow,ncol,nslices,nreps, 'single');
    n=1;
    for r=1:nreps
      for z=1:nslices
        fname = char(fnames{n});
        if ~exist(fname,'file')
          %% todo: error?
          fprintf('%s: ERROR: file %s not found\n',mfilename,fname);
          vol=[];
          return;
        end;
        try
          im = dicomread(fname)';
        catch
          fprintf('%s: ERROR: dicom file %s may be corrupt\n',mfilename,fname);
          vol=[];
          return;
        end;
        if(isempty(im))
          %% todo: error?
          fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
          vol=[];
          return;
        end
        try
          vol(:,:,z,r) = im;
        catch
          fprintf('%s: ERROR: problem with pixel data from %s: %s\n',...
            mfilename,fname,lasterr);
          vol=[];
          return;
        end;
        n=n+1;
      end
    end;  
  end;
end;

