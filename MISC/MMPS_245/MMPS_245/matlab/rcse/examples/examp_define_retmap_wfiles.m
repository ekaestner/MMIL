function retmap = define_retmap_wfiles(smooth);
%function retmap = define_retmap_wfiles([smooth]);
% define retinotopic areas (and extra dipoles) for retinotopy inverse

if ~exist('smooth','var')
  smooth = 0;
end;

subjname = 'larap';
wfilestem = 'retpath';
hemilist = {'lh','rh'};

% dipole cluster "center" number for each theta
%   theta = stim loc numbered ccw 1-16 starting just above right horiz. merid.
%   if value = -1, no cluster (in that hemisphere) for that stim loc
lh_centers = [ 3  2  1  0 -1 -1 -1 -1 -1 -1 -1 -1  7  6  5  4];
rh_centers = [-1 -1 -1 -1  0  1  2  3  4  5  6  7 -1 -1 -1 -1];
% if desired, could allow ipsi activity, but only contra for now

% if near_nbr_weight=1 is passed to ret_inv, it will
%   assume that the locations are sequential, so theta=1:num_locs should
%   be sequential stimulus locations

num_locs = length(lh_centers);
if num_locs ~= length(rh_centers)
  fprintf('%s: error: lh_centers and rh_centers must be same length\n',mfilename);
  return;
end;

cond_order = [1:num_locs];
% cond_order is set so that a subset of locations can be used
%   as deteremined by use_stim_locs passed to ret_inv
%   construct_ret_mapping will remove the locations that are not used
%   and reset cond_order -- this is mainly for plotting purposes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define visual areas

areas(1).name = 'v1';
areas(2).name = 'v2';
areas(3).name = 'v3';
areas(4).name = 'v3a';

for a=1:length(areas)
  for theta=1:num_locs
    lh_c = lh_centers(theta);
    rh_c = rh_centers(theta);
    if lh_c==-1
      areas(a).wfiles(theta).lh_fname = [];
    else
      areas(a).wfiles(theta).lh_fname = ...
        sprintf('%s-%s-center%d-smooth%d-lh.w',...
        wfilestem,areas(a).name,lh_c,smooth);
    end;
    if rh_c==-1
      areas(a).wfiles(theta).rh_fname = [];
    else
      areas(a).wfiles(theta).rh_fname = ...
        sprintf('%s-%s-center%d-smooth%d-rh.w',...
        wfilestem,areas(a).name,rh_c,smooth);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define additional dipoles (optional)

% ret_dips:
ret_dips(1).name = 'ips1';
ret_dips(1).hemisphere = 'lh';
ret_dips(1).v = 15762;

ret_dips(2).name = 'ips1';
ret_dips(2).hemisphere = 'rh';
ret_dips(2).v = 15947;

ret_dips(3).name = 'ips2';
ret_dips(3).hemisphere = 'lh';
ret_dips(3).v = 12384;

ret_dips(4).name = 'ips2';
ret_dips(4).hemisphere = 'rh';
ret_dips(4).v = 22404;

ret_dips(5).name = 'ips3a';
ret_dips(5).hemisphere = 'lh';
ret_dips(5).v = 26807;

ret_dips(6).name = 'ips3b';
ret_dips(6).hemisphere = 'lh';
ret_dips(6).v = 23710;

ret_dips(7).name = 'ips3';
ret_dips(7).hemisphere = 'rh';
ret_dips(7).v = 22245;

ret_dips(8).name = 'lo';
ret_dips(8).hemisphere = 'lh';
ret_dips(8).v = 7226;

ret_dips(9).name = 'lo';
ret_dips(9).hemisphere = 'rh';
ret_dips(9).v = 9686;

ret_dips(10).name = 'mt';
ret_dips(10).hemisphere = 'lh';
ret_dips(10).v = 16132;

ret_dips(11).name = 'mt';
ret_dips(11).hemisphere = 'rh';
ret_dips(11).v = 24041;

ret_dips(12).name = 'v4v';
ret_dips(12).hemisphere = 'lh';
ret_dips(12).v = 19917;

ret_dips(13).name = 'v4v';
ret_dips(13).hemisphere = 'rh';
ret_dips(13).v = 10602;

ret_dips(14).name = 'v6';
ret_dips(14).hemisphere = 'rh';
ret_dips(14).v = 17847;

ret_dips(15).name = 'v7';
ret_dips(15).hemisphere = 'lh';
ret_dips(15).v = 12607;

ret_dips(16).name = 'v7';
ret_dips(16).hemisphere = 'rh';
ret_dips(16).v = 8803;

ret_dips(17).name = 'v8';
ret_dips(17).hemisphere = 'lh';
ret_dips(17).v = 21762;

ret_dips(18).name = 'v8';
ret_dips(18).hemisphere = 'rh';
ret_dips(18).v = 28519;

ret_dips(19).name = 'vip';
ret_dips(19).hemisphere = 'rh';
ret_dips(19).v = 36711;



% nonret_dips:
nonret_dips(1).name = 'ang-cing';
nonret_dips(1).hemisphere = 'rh';
nonret_dips(1).v = 108939;

nonret_dips(2).name = 'ant-LS';
nonret_dips(2).hemisphere = 'lh';
nonret_dips(2).v = 103972;

nonret_dips(3).name = 'ant-mfg';
nonret_dips(3).hemisphere = 'lh';
nonret_dips(3).v = 102482;

nonret_dips(4).name = 'ant-MFG';
nonret_dips(4).hemisphere = 'rh';
nonret_dips(4).v = 106420;

nonret_dips(5).name = 'ant-SFS';
nonret_dips(5).hemisphere = 'lh';
nonret_dips(5).v = 98949;

nonret_dips(6).name = 'ant-SFS';
nonret_dips(6).hemisphere = 'rh';
nonret_dips(6).v = 84691;

nonret_dips(7).name = 'ant-STG';
nonret_dips(7).hemisphere = 'rh';
nonret_dips(7).v = 125576;

nonret_dips(8).name = 'cun';
nonret_dips(8).hemisphere = 'rh';
nonret_dips(8).v = 11462;

nonret_dips(9).name = 'dormed-PPC';
nonret_dips(9).hemisphere = 'lh';
nonret_dips(9).v = 21850;

nonret_dips(10).name = 'dor-PPC';
nonret_dips(10).hemisphere = 'lh';
nonret_dips(10).v = 20400;

nonret_dips(11).name = 'dor-PPC';
nonret_dips(11).hemisphere = 'rh';
nonret_dips(11).v = 12359;

nonret_dips(12).name = 'frontpole';
nonret_dips(12).hemisphere = 'rh';
nonret_dips(12).v = 118014;

nonret_dips(13).name = 'lat-PPC';
nonret_dips(13).hemisphere = 'lh';
nonret_dips(13).v = 30003;

nonret_dips(14).name = 'lat-PPC';
nonret_dips(14).hemisphere = 'rh';
nonret_dips(14).v = 36807;

nonret_dips(15).name = 'pos-LS';
nonret_dips(15).hemisphere = 'rh';
nonret_dips(15).v = 62626;

nonret_dips(16).name = 'precun';
nonret_dips(16).hemisphere = 'rh';
nonret_dips(16).v = 38987;

nonret_dips(17).name = 'sup-frontpole';
nonret_dips(17).hemisphere = 'lh';
nonret_dips(17).v = 119079;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

retmap.cond_order = cond_order;
retmap.num_locs = num_locs;
retmap.areas = areas;
retmap.ret_dips = ret_dips;
retmap.nonret_dips = nonret_dips;

return;

