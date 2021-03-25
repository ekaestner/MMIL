function PING_unzip_anonymize(PING_rsync_dir,PING_Subj_pattern,Setup_script)

% THIS FUNCTION IS OBSOLETE!!!!

% Setup_script - full pathname of our m-file setup script
if ~exist(Setup_script,'file'),
   error(sprintf('The script %s cannot be found.',Setup_script));
end

% Execute Setup_script which specifies some file and dir names
run(Setup_script);

% Create subject exclusion list
ExcludeSubj_List = {''};
if exist(PING_Exclude_file,'file'),
   fid=fopen(PING_Exclude_file,'r');
   ii=0;
   while 1
      subj_id=fgetl(fid);
      if ~ischar(subj_id), break; end
      ii=ii+1;
      ExcludeSubj_List(ii) = {subj_id};
   end
   fclose(fid);
end

% Get our list of subject dirs
Subj_pattern = sprintf('%s/%s',PING_rsync_dir,PING_Subj_pattern);
if isempty(strfind(Subj_pattern,'*'))
   P.isdir = isdir(Subj_pattern);
   P.name = PING_Subj_pattern;
else
   P = dir(Subj_pattern);
end

% Process each subject
for ii=1:length(P),
   if ~P(ii).isdir, continue; end

   SubjID = char(regexp(P(ii).name,'^P\d+$','match'));
   if isempty(SubjID), continue; end

   if ~isempty(strmatch(SubjID,ExcludeSubj_List,'exact')), 
      fprintf(1,'Subject %s found in exclude list. Skipping...\n',SubjID);
      continue; 
   end

   fprintf(1,sprintf('Unzipping and anonymizing %s\n',SubjID));
   PING_Subject_unzip(PING_rsync_dir,SubjID,PING_scratch_dir);

   % Anonymize DICOMs for all visits and attempts for this subject. Each visit_attempt
   % will be in a different dir (e.g. P032300001 vs. P032300002).
   D = dir(sprintf('%s/%s*',PING_scratch_dir,SubjID));
   for jj=1:length(D),
      subj_scratch_dir = sprintf('%s/%s',PING_scratch_dir,D(jj).name);
      if ~isdir(subj_scratch_dir), continue; end
       %  error(sprintf('%s is not a directory.',subj_scratch_dir));

      PING_anonymize(subj_scratch_dir);

      subj_tmp_dir = sprintf('%s/%s',PING_anon_dir,D(jj).name);
      if exist(subj_tmp_dir,'dir'),
         [status,result] = system(sprintf('rm -rf %s',subj_tmp_dir));
         if status, error(result); end
      end

      [status,result] = system(sprintf('cp -rp %s %s/.',subj_scratch_dir,PING_anon_dir));
      if status, error(result); end

      % Remove this scratch dir
      cd(PING_scratch_dir);
      [status,result] = system(sprintf('rm -rf %s',subj_scratch_dir));
      if status, error(result); end
   end
end
