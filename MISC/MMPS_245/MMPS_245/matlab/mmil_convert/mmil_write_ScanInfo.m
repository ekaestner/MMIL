function mmil_write_ScanInfo(ScanInfo,outdir,forceflag)
%function mmil_write_ScanInfo(ScanInfo,[outdir],[forceflag])
%
% Required Input:
%   ScanInfo: scan info struct containing struct arrays for multiple scan types
%
% Optional Input:
%   outdir: output directory
%     {default = pwd}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 03/20/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = sprintf('%s/ScanInfo.csv',outdir);
if ~exist(fname,'file') || forceflag
  fid = fopen(fname,'wt');
  if fid<0, error('failed to open %s for writing',fname); end;
  fprintf(fid,'"ScanType","ScanIndex","SeriesIndex"\n');
  ScanTypes = fieldnames(ScanInfo);
  for t=1:length(ScanTypes)
    type = ScanTypes{t};
    tmpInfo = ScanInfo.(type);
    for s=1:length(tmpInfo)
      fprintf(fid,'"%s",%d,%d\n',type,s,tmpInfo(s).SeriesIndex);
    end;
  end;
  fclose(fid);
end;

