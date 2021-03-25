function [roinums,roilabels,ctab] = fs_read_annotation(fname)
%function [roinums,roilabels,ctab] = fs_read_annotation(fname)
%
% Purpose: reads annotation file
%
% Input:
%   fname: full or relative path of annotation file
%
% Output:
%   roinums: vector of roi numbers (one for each vertex)
%   roilabels: cell array of roi labels (one for every unique roinum)
%   ctab: color table struct
%     ctab.numEntries = number of Entries
%     ctab.orig_tab = name of original ct
%     ctab.struct_names = list of structure names (e.g. central sulcus and so on)
%     ctab.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
%       is b, 4th column is flag, 5th column is resultant integer values
%     calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0
%
% Created:  02/01/07 by Don Hagler
% Last Mod: 11/29/12 by Don Hagler
%

% based on read_annotation.m distributed with freesurfer

roinums = []; roilabels = []; ctab = [];

if ~mmil_check_nargs(nargin,1), return; end;

fp = fopen(fname, 'r', 'b');

if(fp < 0)
  error('failed to open file %s',fname);
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
codevec = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
  fclose(fp);
  error('no colortable found in %s',fname);
end

if(bool)
  %Read colortable
  numEntries = fread(fp, 1, 'int');
  if(numEntries > 0)
    ctab.numEntries = numEntries;
    len = fread(fp, 1, 'int');
    ctab.orig_tab = fread(fp, len, '*char')';
    ctab.orig_tab = ctab.orig_tab(1:end-1);
    ctab.struct_names = cell(numEntries,1);
    ctab.table = zeros(numEntries,5);
    for i = 1:numEntries
      len = fread(fp, 1, 'int');
      ctab.struct_names{i} = fread(fp, len, '*char')';
      ctab.struct_names{i} = ctab.struct_names{i}(1:end-1);
      ctab.table(i,1) = fread(fp, 1, 'int');
      ctab.table(i,2) = fread(fp, 1, 'int');
      ctab.table(i,3) = fread(fp, 1, 'int');
      ctab.table(i,4) = fread(fp, 1, 'int');
      ctab.table(i,5) = ctab.table(i,1) + ctab.table(i,2)*2^8 + ctab.table(i,3)*2^16 + ctab.table(i,4)*2^24;
    end
  else
    version = -numEntries;
    if(version~=2)    
      error('version %d not supported',version);
    else
    end
    numEntries = fread(fp, 1, 'int');
    ctab.numEntries = numEntries;
    len = fread(fp, 1, 'int');
    ctab.orig_tab = fread(fp, len, '*char')';
    ctab.orig_tab = ctab.orig_tab(1:end-1);
    ctab.struct_names = cell(numEntries,1);
    ctab.table = zeros(numEntries,5);
    numEntriesToRead = fread(fp, 1, 'int');
    for i = 1:numEntriesToRead
      structure = fread(fp, 1, 'int')+1;
      if (structure < 0)
        error('bad structure %d',structure);
      end
      if(~isempty(ctab.struct_names{structure}))
        error('duplicated structure %d',structure);
      end
      len = fread(fp, 1, 'int');
      ctab.struct_names{structure} = fread(fp, len, '*char')';
      ctab.struct_names{structure} = ctab.struct_names{structure}(1:end-1);
      ctab.table(structure,1) = fread(fp, 1, 'int');
      ctab.table(structure,2) = fread(fp, 1, 'int');
      ctab.table(structure,3) = fread(fp, 1, 'int');
      ctab.table(structure,4) = fread(fp, 1, 'int');
      ctab.table(structure,5) = ctab.table(structure,1) + ctab.table(structure,2)*2^8 + ctab.table(structure,3)*2^16 + ctab.table(structure,4)*2^24;       
    end
  end    
else
  error('bad value in annotation file');
end

fclose(fp);

codes = ctab.table(:,5);
nrois = length(codes);
roilabels = ctab.struct_names;
roinums = zeros(size(codevec));
for i=1:nrois
  roinums(codevec==codes(i))=i;
end;

