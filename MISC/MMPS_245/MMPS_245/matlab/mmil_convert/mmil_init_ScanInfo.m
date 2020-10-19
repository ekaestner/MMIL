function ScanInfo = mmil_init_ScanInfo(ScanTypes)
%function ScanInfo = mmil_init_ScanInfo(ScanTypes)
%
% Required Input:
%   ScanTypes: cell array of scan types to include as fields
%
% Output:
%   ScanInfo: struct containing empty fields for each of ScanTypes
%
% Created:  09/12/12 by Don Hagler
% Last Mod: 09/12/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

ScanInfo = [];

if ~iscell(ScanTypes), ScanTypes = {ScanTypes}; end;
ntypes = length(ScanTypes);
if ~ntypes, return; end;

tmpinfo = cell(2,ntypes);
tmpinfo(1,:) = ScanTypes;
tmpinfo = reshape(tmpinfo,[ntypes*2,1]);

ScanInfo = struct(tmpinfo{:});

