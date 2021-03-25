function ScanInfo = PING_dir_dcmodify(PINGdir,ScanInfo)

cd(PINGdir);

% dcmodify any DICOM files in the current directory
ScanInfo = PING_dcmodify_byScanInfo(PINGdir,ScanInfo);
if ~isempty(ScanInfo.ErrMsg), return; end

% Recursively go through all subdirectories, modifying any DICOM files found in them
D=dir;
for ii=1:length(D),
   if D(ii).name(1)=='.', continue; end
   if D(ii).isdir,
      subdirname = sprintf('%s/%s',PINGdir,D(ii).name);
      ScanInfo = PING_dir_dcmodify(subdirname,ScanInfo);
      if ~isempty(ScanInfo.ErrMsg), return; end
   end
end
