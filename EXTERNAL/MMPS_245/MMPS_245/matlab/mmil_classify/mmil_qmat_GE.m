function [qmat,nb0,ndd,bval] = mmil_qmat_GE(SeriesInfo) 
%function [qmat,nb0,ndd,bval] = mmil_qmat_GE(SeriesInfo)
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
% Created:  03/05/12 by Vijay Venkatraman
% Last Mod: 02/19/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

nb0 = 0; ndd = 0; qmat= []; bval=[];
for i=1:SeriesInfo.nslices:SeriesInfo.nimgs
  tmpinfo = dicominfo(SeriesInfo.FileNames{i});
  tmp_bval =  tmpinfo.Private_0043_1039(1);
  difdir = [];
  if tmp_bval == 0
    nb0 = nb0 + 1;
    difdir = [0 0 0];
    bval = 0 ; % set zero for non-DTI scans/b0 DTI scan
  else 
    ndd = ndd +1;
    difdir = [tmpinfo.Private_0019_10bb tmpinfo.Private_0019_10bc tmpinfo.Private_0019_10bd];
    if isempty(bval) || bval == 0
      bval = tmp_bval;
    end;
  end;
  qmat = [qmat; difdir(1,:)];
end;    

if ismember(SeriesInfo.SeriesType,{'DTI_ipp','DTI_flex'}) &&...
   SeriesInfo.pepolar>1
  % for integrated scan, first b0 is rev/for scan followed by complete DTI scan
  nb0 = nb0 - 1;
  qmat = qmat(2:end,:);
end;

% set qmat to empty, if no ndd (for/rev and other scans)    
if ndd == 0 
  qmat = [];
end;
    
return;
