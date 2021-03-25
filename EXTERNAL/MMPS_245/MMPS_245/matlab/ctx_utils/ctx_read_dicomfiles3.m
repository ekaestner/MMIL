function ctx_vol = ctx_read_dicomfiles2(dcminfo)

try
  [vol,M] = Convert_Dicom_Images3(dcminfo);
  ctx_vol = ctx_mgh2ctx(vol,M);
  ctx_vol = ctx_standardize_image_orientations(ctx_vol);
catch
  ctx_vol = [];
end

