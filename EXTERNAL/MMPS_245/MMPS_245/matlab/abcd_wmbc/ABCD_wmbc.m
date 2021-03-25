function [vol_corr vol_wm vol_bm vol_gm] = ABCD_wmbc(vol_orig)

[vol_wm vol_bm vol_gm segStruct] = ABCD_wmSeg(vol_orig);
vol_corr = estimate_WM_biasfield_T1_ABCD(vol_orig,vol_wm,vol_bm,vol_gm,3); [vol_corr.minI vol_corr.maxI] = deal(0,255);

