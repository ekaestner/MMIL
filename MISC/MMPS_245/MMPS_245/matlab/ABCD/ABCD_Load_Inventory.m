function inv_info = ABCD_Load_Inventory(ProjID,instem)
%function inv_info = ABCD_Load_Inventory(ProjID,[instem])
%
% Purpose: load inventory from csv, write to mat file for future access
%
% Required Input:
%   ProjID: Project ID string
%
% Optional Input:
%   instem: string included in inventory csv file name
%     NOTE: full name = '/home/{user}/MetaData/{ProjID}/{ProjID}_{instem}.csv
%     {default = 'merged_inventory_pcqcinfo'}
%
% Created:  01/06/16 by Don Hagler
% Last Mod: 01/12/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('instem','var') || isempty(instem)
  instem = 'merged_inventory_pcqcinfo';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indir = sprintf('%s/MetaData/%s',getenv('HOME'),ProjID);
fname_inventory = sprintf('%s/%s_%s.csv',indir,ProjID,instem);

inv_info = abcd_load_csv(fname_inventory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

