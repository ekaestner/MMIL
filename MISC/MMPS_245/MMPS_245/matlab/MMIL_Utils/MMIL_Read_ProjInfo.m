function ProjInfo = MMIL_Read_ProjInfo(fname,numvec_tags)
%function ProjInfo = MMIL_Read_ProjInfo(fname,[numvec_tags])
%
% Required Input:
%   fname: file name of csv file containing information about each project
%     each column will become a field in the ProjInfo output struct
%     each row is from a different subject or exam
%
% Optional Input:
%   numvec_tags: cell array of column headers specifying which columns
%     should be treated as numeric vectors
%
% Created:  08/15/10 by Don Hagler
% Last Mod: 05/02/12 by Don Hagler
%

default_tags = {...
...% MRI
  'DTI_trans'...
  'DTI_rot'...
  'DTI_resolution'...
  'DTI_nvoxels' ...
  'BOLD_resolution'...
  'BOLD_nvoxels' ...
  'kernelWidthMax_vec'...
  'lambda2_vec'...
  'DTI_CSD_seed_point_sampling'...
...% MEG
  'PROC_valid_event_codes'...
  'RF_pol_snums'...
  'RF_ecc_snums'...
  'RCSE_r_offset_range'...
  'RCSE_th_offset_range'...
};

ProjInfo = [];
if (~mmil_check_nargs(nargin,1)) return; end;
if ~exist('numvec_tags','var') || isempty(numvec_tags)
  numvec_tags = default_tags;
else
  if ~iscell(numvec_tags)
    numvec_tags = {numvec_tags};
  end;
  numvec_tags = reshape(numvec_tags,[1,numel(numvec_tags)]);
  numvec_tags = cat(2,numvec_tags,default_tags);
end;

if ~exist(fname,'file'), error('file %s not found',fname); end;

% creates struct array with field names from the column headers
ProjInfo = mmil_csv2struct(fname);

% exclude comments and blank lines
ID = 'ProjID';
if ~isfield(ProjInfo,ID)
  error('ProjInfo file %s is missing ProjID column');
end;
ProjIDs = ProjInfo.ProjID;
ind_keep = [];
for i=1:length(ProjInfo)
  if ~isempty(ProjInfo(i).(ID)) & ~strncmp(ProjInfo(i).(ID),'%',1)
    ind_keep = [ind_keep, i];
  end;
end;
ProjInfo = ProjInfo(ind_keep);
             
% convert certain fields from strings to numeric vectors
ProjInfo = mmil_structarr_str2num(ProjInfo,numvec_tags);
