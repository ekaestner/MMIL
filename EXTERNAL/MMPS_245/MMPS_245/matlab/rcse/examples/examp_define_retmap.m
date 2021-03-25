function retmap = define_retmap();
%function retmap = define_retmaps();
% define retinotopic areas (and extra dipoles) for retinotopy inverse

subjname = 'larap';
wfilestem = 'retpath';
hemilist = {'lh','rh'};

% dipole cluster "center" number for each theta
%   theta = stim loc numbered ccw 1-16 starting just above right horiz. merid.
%   if value = -1, no cluster (in that hemisphere) for that stim loc
lh_locs = [ 3  2  1  0 -1 -1 -1 -1 -1 -1 -1 -1  7  6  5  4];
rh_locs = [-1 -1 -1 -1  0  1  2  3  4  5  6  7 -1 -1 -1 -1];
% if desired, could allow ipsi activity, but only contra for now

% if near_nbr_weight=1 is passed to ret_inv, it will
%   assume that the locations are sequential, so theta=1:num_locs should
%   be sequential stimulus locations

num_locs = length(lh_locs);
if num_locs ~= length(rh_locs)
  fprintf('%s: error: lh_locs and rh_locs must be same length\n',mfilename);
  return;
end;

cond_order = [1:num_locs];
% cond_order is set so that a subset of locations can be used
%   as deteremined by use_stim_locs passed to calc_ret_inverse
%   construct_ret_mapping will remove the locations that are not used
%   and reset cond_order -- this is mainly for plotting purposes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define visual areas

areas(1).name = 'v1';
areas(2).name = 'v2';
areas(3).name = 'v3';
areas(4).name = 'v3a';

% dipole "center" vertices, arranged from upper to lower field
%  note: this order is arbitrary and for convenience
centers(1).lh = [6354,5516,7233,7279,4067,2956,2514,2487];
centers(1).rh = [8384,9795,8381,5678,2837,3318,5007,4966];

centers(2).lh = [8055,11199,13750,16284,26,168,756,1150];
centers(2).rh = [9822,10564,13811,15598,1306,1734,2732,3260];

centers(3).lh = [20823,126447,18958,16207,206,191,740,1124];
centers(3).rh = [11301,12114,12953,14715,936,646,1279,1708];

centers(4).lh = [9338,7836,7069,6202,4665,4681,3978,3399];
centers(4).rh = [4825,4817,4805,3704,3179,2680,3193,3204];

% check for errors
tmp_numlocs = length(centers(1).lh);
for j=1:length(areas)
  if length(centers(j).lh) ~= tmp_numlocs | ...
     length(centers(j).rh) ~= tmp_numlocs
    fprintf('%s: error: areas should each have the same number of dipoles!\n',...
      mfilename);
    return;
  end;
end;

for a=1:length(areas)
  for theta=1:num_locs
    lh_loc = lh_locs(theta);
    rh_loc = rh_locs(theta);
    if lh_loc==-1
      areas(a).verts(theta).v_lh = [];
      areas(a).verts(theta).w_lh = [] ;
    else
      areas(a).verts(theta).v_lh = centers(a).lh(lh_loc+1);
      areas(a).verts(theta).w_lh = 1;
      % vertices could be differntially weighted (e.g. contra/ipsi)
      % and a vector of vertices can be specified
    end;
    if rh_loc==-1
      areas(a).verts(theta).v_rh = [];
      areas(a).verts(theta).w_rh = [] ;
    else
      areas(a).verts(theta).v_rh = centers(a).rh(rh_loc+1);
      areas(a).verts(theta).w_rh = 1;
    end;
  end;
end;

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

