rawdatadir = '/space/pogo1/6/dhagler/meg-raw';
session = '060325LP';
datafile = sprintf('%s/%s/%s_raw.fif',rawdatadir,session,session);

ref_EEG_coords=ts_extract_hpts_from_fiff(datafile)

