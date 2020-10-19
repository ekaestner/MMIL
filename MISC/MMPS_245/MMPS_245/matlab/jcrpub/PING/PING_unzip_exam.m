function [unzippath,ret_status,result] = PING_unzip_exam(indir,unzip_topdir,outdir)
% [unzippath,ret_status,result] = PING_unzip_exam(indir,unzip_topdir,outdir)

unzippath = '';
ret_status = 0;
result = '';

if ~exist(indir,'dir'),
   result = sprintf('The dir %s does not exist.',indir);
   ret_status = 1;
   return;
end

if ~exist(unzip_topdir,'dir')
   result = sprintf('The dir %s does not exist.',unzip_topdir);
   ret_status = 1;
   return;
end

% Check if all required files present
cd(indir);
Z = dir('*.zip');
if length(Z) ~= 1,
  result = sprintf('There must be a single .zip file in %s.\n',indir);
  ret_status = 1;
  return;
end
zipfile = Z.name;
if ~exist(zipfile,'file') | ~exist('ScanInfo.txt','file') | ~exist('UploadStatus.txt','file'),
   result = sprintf('Incomplete file set in %s.\n',indir);
   ret_status = 1;
   return;
end

% Check UploadStatus
statusfile = sprintf('%s/UploadStatus.txt',indir);
fid=fopen(statusfile,'r');
tline = fgetl(fid);
fclose(fid);
if tline ~= '1',
   result = sprintf('%s: UploadStatus failure.\n',statusfile);
   ret_status = 1;
   return;
end

[status,ret_str] = unix(sprintf('mktemp -d -p %s %s.XXXX',unzip_topdir,outdir));
if status,
   result = sprintf('%s: Cannot create scratch dir for %s',mfilename,outdir);
   ret_status = 1;
   return;
end
unzippath = deblank(ret_str);

% unzip zipfile contents to unzippath
zipfilepath = sprintf('%s/%s',indir,zipfile);
unzip_cmd = sprintf('unzip -q -u "%s" -d %s',zipfilepath,unzippath);
[status,res] = system(unzip_cmd);
if status, 
   result = res;
   ret_status = 1;
   return;
end

% Copy ScanInfo.txt
[status,res] = system(sprintf('cp %s/ScanInfo.txt %s',indir,unzippath));
if status,
   result = res;
   ret_status = 1;
   return;
end
