function [Pheno] = extract_alex_pheno(asegstats,SubjectIDs,outfile)

write_flag=false;
if nargin==3
   write_flag=true;
end

leftamyg_ind = strmatch('Left-Amygdala',asegstats.colnames);
rightamyg_ind = strmatch('Right-Amygdala',asegstats.colnames);
lefthippo_ind =  strmatch('Left-Hippocampus',asegstats.colnames);
righthippo_ind =  strmatch('Right-Hippocampus',asegstats.colnames);

leftcerewhite_ind = strmatch('Left-Cerebellum-White-Matter',asegstats.colnames);
rightcerewhite_ind = strmatch('Right-Cerebellum-White-Matter',asegstats.colnames);
leftcerectx_ind =  strmatch('Left-Cerebellum-Cortex',asegstats.colnames);
rightcerectx_ind =  strmatch('Right-Cerebellum-Cortex',asegstats.colnames);

leftcaud_ind = strmatch('Left-Caudate',asegstats.colnames);
rightcaud_ind = strmatch('Right-Caudate',asegstats.colnames);
leftputa_ind = strmatch('Left-Putamen',asegstats.colnames);
rightputa_ind = strmatch('Right-Putamen',asegstats.colnames);
leftpall_ind = strmatch('Left-Pallidum',asegstats.colnames);
rightpall_ind = strmatch('Right-Pallidum',asegstats.colnames);

numsubj = length(SubjectIDs);

amygdata = asegstats.table(:,[leftamyg_ind rightamyg_ind]);
hippodata = asegstats.table(:,[lefthippo_ind righthippo_ind]);
ceredata = asegstats.table(:,[leftcerewhite_ind rightcerewhite_ind leftcerectx_ind rightcerectx_ind]);
bgdata = asegstats.table(:,[leftcaud_ind rightcaud_ind leftputa_ind rightputa_ind leftpall_ind rightpall_ind]);

amygvol = sum(amygdata,2);
hippovol = sum(hippodata,2);
cerevol = sum(ceredata,2);
bgvol = sum(bgdata,2);

Pheno.dat = [amygvol, hippovol, cerevol, bgvol];
Pheno.names = [{'AmygdalaVol'},{'HippocampusVol'},{'CerebellumVol'},{'BasalGangliaVol'}];
Pheno.subjID = SubjectIDs;

if write_flag,
   fid=fopen(outfile,'w');
   if fid==-1,
      error(sprintf('Cannot open %s for writing', outfile));
   end
   fprintf(fid,'SubjectID'); 
   fprintf(fid,',%s',Pheno.names{:});  fprintf('\n');
   for ii=1:numsubj,   
      fprintf(fid,'%d,%.1f,%.1f,%.1f,%.1f\n', str2double(SubjectIDs(ii)), amygvol(ii), hippovol(ii), cerevol(ii), bgvol(ii));
   end
   fclose(fid);
end

