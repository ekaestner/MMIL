procrootdir = '~/Dropbox/data/ABCD/proc';

subdirs = {'MRIPROC_G010_INV1W5TLB74_20170819_20170819.093848_1' ...
           'MRIPROC_G031_INV0PW5A9P1_20170816_20170816.180929_1' ...
           'MRIPROC_S076_INV31JMTLE9_20170817_20170817.121618.835000_1'};

% These do not have the occipital intensity issue -- check with RAs what the issue is
% 'MRIPROC_P043_INV0M9U5RX8_20161119_20161119.151150_1' % 

forceflag = true;

for subdiri = 1:length(subdirs)
  subdir = subdirs{subdiri};
  procdir = sprintf('%s/%s',procrootdir,subdir);
  fname_snap = sprintf('%s/snap.mat',procdir);
  if ~exist(fname_snap,'file') | forceflag
    fname_orig = sprintf('%s/MPR1_uw.mgz',procdir); 
    fname_wmbc = sprintf('%s/MPR1_uw_wmbc.mgz',procdir); 
    if ~exist(fname_orig,'file')
      fname_orig = sprintf('%s/MPR1.mgz',procdir); 
      fname_wmbc = sprintf('%s/MPR1_wmbc.mgz',procdir); 
    end
    vol_orig = ctx_load_mgh(fname_orig); 
    vol_wmbc = ctx_load_mgh(fname_wmbc); 
    [vol_wmbc.minI vol_wmbc.maxI] = deal(0,255);
    [vol_wm vol_bm vol_gm segStruct] = wmSeg(vol_orig);
    vol_corr = estimate_WM_biasfield_T1_ABCD(vol_orig,vol_wm,vol_bm,vol_gm,3); [vol_corr.minI vol_corr.maxI] = deal(0,255);
    [vol_wm2 vol_bm2 vol_gm2 segStruct2] = wmSeg(vol_corr);
    vol_corr2 = estimate_WM_biasfield_T1_ABCD(vol_orig,vol_wm2,vol_bm2,vol_gm2,3); [vol_corr2.minI vol_corr2.maxI] = deal(0,255);
    [vol_wm3 vol_bm3 vol_gm3 segStruct3] = wmSeg(vol_corr2);
    vol_corr3 = estimate_WM_biasfield_T1_ABCD(vol_orig,vol_wm3,vol_bm3,vol_gm3,3); [vol_corr3.minI vol_corr3.maxI] = deal(0,255);
    save(fname_snap,'vol_orig','vol_wmbc','vol_corr','vol_corr2','vol_corr3','segStruct','segStruct2','segStruct3');
%    showVol(vol_orig,vol_wmbc,vol_corr,vol_corr2,vol_corr3);
  else
    disp(fname_snap);
    load(fname_snap);
    showVol(vol_orig,vol_wmbc,vol_corr,vol_corr2,vol_corr3);
    pause
  end
end


% Look at single subject
ubdiri = 1;
subdir = subdirs{subdiri};
procdir = sprintf('%s/%s',procrootdir,subdir);
fname_snap = sprintf('%s/snap.mat',procdir);
disp(fname_snap);
load(fname_snap);
showVol(vol_orig,vol_wmbc,vol_corr,vol_corr2,vol_corr3,segStruct.topLabel,segStruct2.topLabel,segStruct3.topLabel);

[vol_wm vol_bm vol_gm segStruct] = wmSeg(vol_orig); % Should modify wmSeg to return vol_wm vol_bm vol_gm, if passed segStruct
save(fname_snap,'vol_orig','vol_wmbc','segStruct','vol_wm','vol_bm','vol_gm'); % Save progress so far
vol_corr_new = estimate_WM_biasfield_T1_ABCD(vol_orig,vol_wm,vol_bm,vol_gm,3); [vol_corr_new.minI vol_corr_new.maxI] = deal(0,255);
showVol(vol_orig,vol_wmbc,vol_corr,vol_corr_new);

% ToDo
%   Look into estimate_WM_biasfield, check wm mask volume, post-censoring
%     try sparse smoothing of sign of values relative to "target" value, increase weight of voxels with excessive signal intensity
%     also look at sparsely smoothed GM intensity?
%   Fix transpose -- make new version of ctx_load_mgh

