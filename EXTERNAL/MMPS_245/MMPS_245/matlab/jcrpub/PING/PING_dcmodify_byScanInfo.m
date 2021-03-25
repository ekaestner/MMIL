function [ScanInfo] = PING_dcmodify_byScanInfo(dicom_dir,ScanInfo)

errmsg = nargchk(2,2,nargin);
if ~isempty(errmsg)
   ScanInfo.ErrMsg = errmsg;
   return;
end

%dcmtk_dir = '/usr/pubsw/packages/dcmtk/3.5.3/bin';
dcmtk_dir = '/usr/pubsw/packages/dcmtk/3.6.0/bin';

dcmodify_cmd = sprintf('%s/dcmodify',dcmtk_dir);
%dcmdump_cmd = sprintf('%s/dcmdump',dcmtk_dir);

anon_dcm_cmd = '/home/cooper/bin/PING_anonymize_dcm_DIR.sh';

if ~isdir(dicom_dir),
   ScanInfo.ErrStr = sprintf('%s: ERROR: %s is not a directory.',mfilename,dicom_dir);
   return;
end

ind1 = regexp(ScanInfo.ScheduledScanDate,'\d\d\d\d-\d\d-\d\d');
ind2 = regexp(ScanInfo.DOB,'\d\d/\d\d\d\d');

if isempty(ind1) || ind1~=1 || isempty(ind2) || ind2~=1,
   age = NaN;
else
   try age = (datenum(ScanInfo.ScheduledScanDate,'yyyy-mm-dd') - datenum(ScanInfo.DOB,'mm/yyyy') - 15)/365.25;
   catch age = NaN;  end
end
if isnan(age) || age < 0 || age > 200,
  age_str = '';
else
  age_str = sprintf('%.2fY',age);
end

if ~isempty(ind2) && ind2==1,
  DOB = datestr(datenum(ScanInfo.DOB,'mm/yyyy'),'yyyymm');
else
  DOB = '';
end

dcm_mod_cmd = sprintf('%s -nb -i ''(0008,0080)=%s -i (0008,0081)= -i (0010,0010)=%s -i (0010,0020)=%s -i (0010,0030)= -i (0010,1010)= -i (0010,0040)=''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID);

%dcm_mod_cmd = sprintf('%s -nb -i ''"(0008,0080)=%s" -i "(0008,0081)=" -i "(0010,0010)=%s" -i "(0010,0020)=%s" -i "(0010,0030)=" -i "(0010,1010)=" -i "(0010,0040)="''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID);

% Build dcmodify cmd to modify Patient ID, scan site, DOB, age, and sex.
%dcm_mod_cmd = sprintf('%s -i ''"(0008,0080)=%s" -i "(0008,0081)=" -i "(0010,0010)=%s" -i "(0010,0020)=%s" -i "(0010,0030)=%s" -i "(0010,1010)=%s" -i "(0010,0040)=%s"''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID,DOB,age_str,ScanInfo.Sex);

% The -nb flag specifies "no backup files created"
%dcm_mod_cmd = sprintf('%s -nb -i ''(0008,0080)=%s -i (0008,0081)= -i (0010,0010)=%s -i (0010,0020)=%s -i (0010,0030)=%s -i (0010,1010)=%s -i (0010,0040)=%s''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID,DOB,age_str,ScanInfo.Sex);

%dcm_mod_cmd = sprintf('%s -ma ''"(0008,0080)=%s" -ma "(0008,0081)=" -ma "(0010,0010)=%s" -ma "(0010,0020)=%s" -ma "(0010,0030)=%s" -i "(0010,1010)=%s" -i "(0010,0040)=%s"''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID,DOB,age_str,ScanInfo.Sex);

% Use the following version if passing this string as an arg to a shell script
%dcm_mod_cmd = sprintf('%s -ma ''''InstitutionName=%s -ma InstitutionAddress=''''',dcmodify_cmd,ScanInfo.SiteID);
%dcm_mod_cmd = sprintf('%s -ma ''''InstitutionName=%s -ma InstitutionAddress= -ma PatientsName=%s -ma PatientID=%s -ma PatientsBirthDate=%s -i PatientsAge=%s -i PatientsSex=%s''''',dcmodify_cmd,ScanInfo.SiteID,ScanInfo.SubjectID,ScanInfo.SubjectID,DOB,age_str,ScanInfo.Sex);

[status,result] = unix(sprintf('%s "%s" %s',anon_dcm_cmd,dicom_dir,dcm_mod_cmd));

if (status == 0), 
   [mod_dcm_count,cnt]=sscanf(result,'anon_count=%d');
   if cnt==1 && mod_dcm_count>=0,  
      ScanInfo.AnonFileCount = ScanInfo.AnonFileCount + mod_dcm_count; 
   end
else
   result = deblank(result);
   ind = strfind(result,sprintf('\n'));
   result(ind) = '|';
   ScanInfo.ErrMsg = sprintf('Error anonymizing files in %s: %s\n',dicom_dir,result); 
   return;
end
