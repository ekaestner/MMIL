function ctx_vol = ctx_read_dicomvol(fnames)
% Last Mod: 12/21/12 by Don Hagler

[vol,M] = mmil_read_dicom_vol(fnames);
ctx_vol = ctx_mgh2ctx(vol,M);
ctx_vol = ctx_standardize_image_orientations(ctx_vol);

