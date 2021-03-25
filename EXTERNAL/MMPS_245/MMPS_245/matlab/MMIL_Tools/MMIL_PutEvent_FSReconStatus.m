function errcode = MMIL_PutEvent_FSReconStatus(ContainerPath,ReconStatus,ReconMessage);
%function errcode = MMIL_PutEvent_FSReconStatus(ContainerPath,ReconStatus,ReconMessage);
%
% created:  01/11/07 Don Hagler
% last mod: 05/25/07 Don Hagler
%

if nargin<3, help(mfilename); return; end;

errcode = 0;
quiet_flag = 0;

if ~ReconStatus, return; end;
fname = sprintf('%s/scripts/EVENT.xml',ContainerPath);

% write to xml file
fid = fopen(fname,'wt');
if fid==-1
  errcode = 1;
  fprintf('%s: ERROR: unable to write to file %s\n',...
    mfilename,fname);
  return;
end;
fprintf(fid,'<ReconInfo>\n');
fprintf(fid,'  <ReconStatus>%d</ReconStatus>\n',ReconStatus);
fprintf(fid,'  <ReconMessage>%s</ReconMessage>\n',ReconMessage);
fprintf(fid,'</ReconInfo>\n');
fclose(fid);

if ~quiet_flag
  fprintf('%s: %s: ReconStatus = %s\n',...
    mfilename,ContainerPath,ReconMessage);
end;

