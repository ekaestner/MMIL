function unwarpflag = mmil_detect_gradwarp(dcmfname)
%function unwarpflag = mmil_detect_gradwarp(dcmfname)
%
% loads dicom file, looks for more than 10% zeros on edges
%   (indicates gradwarp applied)
%   edges with all zeros will be discarded before testing fraction of zeros
%
%   unwarpflag=1  (gradwarp has been applied)
%   unwarpflag=0  (gradwarp has not been applied)
%
% Created:  12/19/08 by Don Hagler
% Last Mod: 01/27/17 by Don Hagler
%

unwarpflag = 0;

thresh = 0.1;

if ~exist(dcmfname,'file')
  error('file %s not found',dcmfname);
end;
try
  im = dicomread(dcmfname);
catch
  error('unable to read file %s as dicom',dcmfname);
end;

[nr,nc] = size(im);
nr2 = round(nr/2);
nc2 = round(nc/2);

edge1 = nan;
for i=1:nr2
  edge1 = mmil_rowvec(im(i,:));
  if ~all(edge1==0), break; end;
end;

edge2 = nan;
for i=nr:-1:nr2
  edge2 = mmil_rowvec(im(i,:));
  if ~all(edge2==0), break; end;
end;

edge3 = nan;
for i=1:nc2
  edge3 = mmil_rowvec(im(:,i));
  if ~all(edge3==0), break; end;
end;

edge4 = nan;
for i=nc:-1:nc2
  edge4 = mmil_rowvec(im(:,i));
  if ~all(edge4==0), break; end;
end;

border = cat(2,edge1,edge2,edge3,edge4);

if any(isnan(border)), return; end;

if length(find(border==0))/length(border) > thresh
  unwarpflag = 1;
else
  unwarpflag = 0;
end;

