function [qmat,nb0,ndd,bval] = mmil_qmat_Siemens(SeriesInfo) 
%function [qmat,nb0,ndd,bval] = mmil_qmat_Siemens(SeriesInfo)
%
% Purpose: gather qmat information from SeriesInfo 
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   qmat: diffusion direction matrix
%   nb0: number of b0 images
%   ndd: number of diffusion directions
%   bval: b value
%
% Created:  08/16/11 by Vijay Venkatraman
% Prev Mod: 04/01/15 by Don Hagler
% Last Mod: 07/08/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

nb0 = 0; ndd=0; qmat = []; bval = []; bvals = [];
csa_flag = 1;
for i=1:SeriesInfo.nimgs
  tmpinfo = dicominfo(SeriesInfo.FileNames{i});
  tmp_bval = tmpinfo.Private_0019_100c;
  diffdirs = [];
  if isa(tmp_bval,'uint8'), tmp_bval = str2double(char(tmp_bval)'); end
  if tmp_bval == 0
    nb0 = nb0+1;
    diffdirs = [0 0 0];
    bval = 0;
  else
    ndd = ndd+1;
    if csa_flag
      try
        % read Siemens CSA header to find diffusion directions
        %   CSA header is more accurate than Private_0019_100e according to:
        %   http://www.cabiatl.com/mricro/mricron/bvec.doc
        if ndd==1
          fprintf('%s: reading diffusion directions from CSA header...\n',...
            mfilename);
        end;
        csaheader = spm_dicom_headers(SeriesInfo.FileNames{i});
        names = {csaheader{1}.CSAImageHeaderInfo.name};
        idx = find(strcmp(names,'DiffusionGradientDirection'));
        if isempty(idx)
          fprintf('%s: WARNING: missing DiffusionGradientDirection in CSA header from %s\n',...
            mfilename,SeriesInfo.FileNames{i});
          csa_flag = 0;
        end;
        diffdirs = [...
          csaheader{1}.CSAImageHeaderInfo(idx).item(1).val...
          csaheader{1}.CSAImageHeaderInfo(idx).item(2).val...
          csaheader{1}.CSAImageHeaderInfo(idx).item(3).val];
        diffdirs = str2num(diffdirs);    
      catch
        fprintf('%s: WARNING: failed to read CSA header from %s\n',...
          mfilename,SeriesInfo.FileNames{i});
        csa_flag = 0;
      end;
    end;
    if ~csa_flag && isfield(tmpinfo, 'Private_0019_100e') &&...
     isa(tmpinfo.Private_0019_100e,'double')
      if ndd==1
        fprintf('%s: reading diffusion directions from Private_0019_100e...\n',...
          mfilename);
      end;
      % if diffdirs is empty from csaheader than use Private_0019_100e 
      diffdirs = tmpinfo.Private_0019_100e;
      diffdirs = reshape(diffdirs,1,3);
    end
    if isempty(bval) || bval == 0
      bval = tmp_bval;
    end;
  end;
  if isempty(diffdirs)
    fprintf('%s: WARNING: failed to read diffusion directions from %s\n',...
      mfilename,SeriesInfo.FileNames{i});
    return;
  end;
  qmat(i,:) = diffdirs(1,:);
  bvals(i,:) = tmp_bval;
  if ndd+nb0>=SeriesInfo.nimgs
    break;
  end;
end;
% flipping directions x and y to make it RAS
qmat(:,1) = -qmat(:,1);
qmat(:,2) = -qmat(:,2);

% if multi-shell, scale qmat, reset bval
if length(unique(bvals(bvals>0)))>1
  bval = max(bvals);
  qmat = bsxfun(@times,qmat,sqrt(bvals/bval));
end;

return;

