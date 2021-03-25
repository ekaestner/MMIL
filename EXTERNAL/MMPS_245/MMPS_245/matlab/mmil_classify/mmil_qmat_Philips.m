function [qmat,nb0,ndd,bval] = mmil_qmat_Philips(SeriesInfo) 
%function [qmat,nb0,ndd,bval] = mmil_qmat_Philips(SeriesInfo)
%
% Purpose: gather qmat information from SeriesInfo 
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   qmat: Qmat file
%   nb0: number of b0 images
%   ndd: number of diffusion directions
%   bval: b value
%
% Created:  03/01/12 by Vijay Venkatraman
% Prev Mod: 06/28/17 by Don Hagler
% Last Mod: 07/29/17 by Don Hagler
%

% NOTE: for Philips DTI scans, instance numbers are ordered by
%     slice (outer loop), then diffusion direction (inner loop)
%  on second occurence of b=0 image, instance number jumps (e.g. to 99)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

nb0 = 0; ndd = 0; qmat = []; bval = []; bvals = [];
smf = 1e-5;

% get position of first slice
tmpinfo = dicominfo(SeriesInfo.FileNames{1});
slice_pos_first = mmil_getfield(tmpinfo,'ImagePositionPatient');
if isempty(slice_pos_first)
  fprintf('%s: WARNING: failed to get qmat because of missing ImagePositionPatient\n',mfilename);
  return;
end;

for i=1:SeriesInfo.nimgs
  tmpinfo = dicominfo(SeriesInfo.FileNames{i});
	% stop if we've reached second slice
	slice_pos = tmpinfo.ImagePositionPatient;
	if any(slice_pos-slice_pos_first>smf)
		break;
	end;
  diffdirs = [];
  if tmpinfo.DiffusionBValue == 0
    nb0 = nb0 + 1;
    diffdirs = [0 0 0];
		tmp_bval = 0;
		if isempty(bval)
			bval = 0;
		end;
  else
    ndd = ndd + 1;
    tmp_bval = tmpinfo.DiffusionBValue;
		if tmp_bval > bval
			bval = tmp_bval;
		end;
    if isfield(tmpinfo,'DiffusionGradientOrientation') & ...
        isa(tmpinfo.DiffusionGradientOrientation,'double')
      diffdirs = tmpinfo.DiffusionGradientOrientation;
      diffdirs = reshape(diffdirs,1,3);
    end;
  end;
  qmat(i,:) = diffdirs(1,:);
  bvals(i,:) = tmp_bval;
end;

if bvals(end)>0 && all(qmat(end,:)==[0 0 0])
  % last scan is a synthesized image
	skip_last_flag = 1;
else
	skip_last_flag = 0;
end;

if ~isempty(qmat)
	if skip_last_flag
	  qmat = qmat(1:end-1,:);
  	ndd = ndd-1;
	end;
  % X and Y are flipped on DICOM qmat compared to Philips_difdirs_32.txt qmat
  qmat(:,1) = -qmat(:,1);
  qmat(:,2) = -qmat(:,2);
end;

% set qmat to empty, if no ndd (for/rev and other scans)    
if ndd == 0
  qmat = [];
end;

% if multi-shell, scale qmat, reset bval
if length(unique(bvals(bvals>0)))>1
  bval = max(bvals);
  qmat = bsxfun(@times,qmat,sqrt(bvals/bval));
end;

return;

