function [vol_wm vol_bm vol_gm segStruct] = ABCD_wmSeg(vol_T1)

if isfield(vol_T1,'imgs')
  vol_T1.imgs = double(vol_T1.imgs); % Should only do this if not already double
  segStruct = ABCD_corSeg(vol_T1,'reslice off'); 
else
  segStruct = vol_T1;
end
vol_bm = segStruct.brainmask;
vol_wm = vol_T1;
vol_wm.imgs = ismember(segStruct.topLabel.imgs,[1 9 17 19 21]); vol_wm.maxI = 1; % Include WM hypointensities
vol_gm = vol_T1;
vol_gm.imgs = ismember(segStruct.topLabel.imgs,[2 11 12 14 15]); vol_gm.maxI = 1;

%keyboard

% Should look into ~/matlab_new/DynamicAtlas/oldwork/corSeg_dynatl.m
% Look into /Users/dale/matlab_big/NQ_SegEdit/ctxatl_1_4_rbm_edited.mat for edited atlas
% Look at (prior) probability volumes of curently used atlas  
% Where are the intensity targets per ROI?

