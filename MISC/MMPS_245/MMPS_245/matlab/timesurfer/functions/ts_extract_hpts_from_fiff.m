function ref_EEG_coords = ts_extract_hpts_from_fiff(datafile,hptsfile)
%function ref_EEG_coords = ts_extract_hpts_from_fiff(datafile,[hptsfile])
%
% Required Input:
%   datafile: input fif format MEG/EEG data file
%
% Optional Input:
%   hptsfile: output hpts file
%     If empty, will replace datafile fif extension with hpts
%     {default = []}
%
% uses Uutela's meg_pd fiff access toolbox hpipoints function
%
% Early Mod: 08/05/09 by Don Hagler
% Last Mod:  04/04/14 by Don Hagler
%

ref_EEG_coords=[];
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('hptsfile','var'), hptsfile = []; end;

if ~exist(datafile,'file')
  error('file %s not found',datafile);
end;
if isempty(regexp(datafile,'\.fif$'))
  error('datafile %s is missing .fif extension',datafile);
end;
if isempty(hptsfile)
  hptsfile = regexprep(datafile,'.fif','.hpts');
end;

[coords,kind,num] = hpipoints(datafile);

points = [];
for p=1:length(num)
  switch kind(p)
    case 1
      points(p).name = 'cardinal';
    case 2
      points(p).name = 'hpi';
    case 3
      points(p).name = 'eeg';
    case 4
      points(p).name = 'extra';
    otherwise
      points(p).name = [];
  end;
  points(p).num = num(p);
  points(p).x = coords(1,p);
  points(p).y = coords(2,p);
  points(p).z = coords(3,p);
end

% get coordinates of EEG reference
eeg_i = find(strcmp('eeg',{points.name}));
eeg0_i = find(cell2mat({points(eeg_i).num})==0);
if ~isempty(eeg0_i)
  p = eeg_i(eeg0_i);
  ref_EEG_coords=[points(p).x,points(p).y,points(p).z];
  points(p).name = 'eeg_ref';
end;

npoints = length(points);

fid=fopen(hptsfile,'wt');
for p=1:npoints
  switch points(p).name
    case {'cardinal','extra','eeg','eeg_ref'}
      fprintf(fid,'%s %03d %f %f %f\n',...
        points(p).name,...
        points(p).num,...
        points(p).x,...
        points(p).y,...
        points(p).z);
  end;
end;
fclose(fid);

