function ctx_vol = ctx_read_dicomfiles(dcminfo)

[vol,M] = Convert_Dicom_Images(dcminfo);
ctx_vol = ctx_mgh2ctx(vol,M);
