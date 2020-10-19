function mmil_merge_csv(fname_in1,fname_in2,fname_out,...
    merge_field,merge_flag,forceflag)
%function mmil_merge_csv(fname_in1,fname_in2,fname_out,...
%   [merge_field],[merge_flag],[forceflag])
%
% purpose: merge two csv files
%
% required input:
%   fname_in1: file name of first input csv spreadsheet
%   fname_in2: file name of second input csv spreadsheet
%   fname_out: output file name
%
% optional input:
%   merge_field: string specifiying name of column on which to merge
%     if not supplied, will use first column header in common
%     {default = []}
%   merge_flag: how to handle missing rows in each cell array
%     0: intersection of rows in fname_in1 and fname_in2
%     1: all rows in fname_in1, with blanks for missing vals in fname_in2
%     2: all rows in fname_in2, with blanks for missing vals in fname_in1
%     3: union of rows in fname_in1 and fname_in2
%     {default = 0}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  12/06/12 by Don Hagler
% Last Mod: 07/19/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('merge_field','var'), merge_field = []; end;
if ~exist('merge_flag','var') || isempty(merge_flag), merge_flag = 0; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

if ~exist(fname_out,'file') || forceflag
  if ~exist(fname_in1,'file'), error('file %s not found',fname_in1); end;
  if ~exist(fname_in2,'file'), error('file %s not found',fname_in2); end;
  outdir = fileparts(fname_out);
  if ~isempty(outdir)
    mmil_mkdir(outdir);
  end;
  vals1 = mmil_readtext(fname_in1);
  vals2 = mmil_readtext(fname_in2);
  vals = mmil_merge_cells(vals1,vals2,merge_field,merge_flag);
  mmil_write_csv(fname_out,vals);
end;

