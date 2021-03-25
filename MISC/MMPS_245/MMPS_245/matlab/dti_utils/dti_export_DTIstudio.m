function dti_export_DTIstudio(vol,fname,SliceThickness,readoutFOV,phaseFOV,...
  bval,qmat,sliceorient,sliceseq)
%function dti_export_DTIstudio(vol,fname,SliceThickness,readoutFOV,phaseFOV,...
%  bval,[qmat],[sliceorient],[sliceseq])
%
% Created:  12/12/06 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,6), return; end;

if ~exist('qmat','var'), qmat=[]; end;
if ~exist('sliceorient','var'), sliceorient=1; end;
if ~exist('sliceseq','var'), sliceseq=0; end;

k=findstr(fname,'.rec');
if ~isempty(k)
  fstem = fname(1:k(end)-1);
else
  fstem = fname;
end;

% write out diffusion data in Philip's rec format (readable by DTI studio)
vol(find(isnan(vol)))=0;
vol = permute(vol,[2,1,3,4]);;
imgfilename=sprintf('%s.rec',fstem);
fimg=fopen(imgfilename,'w');
count=fwrite(fimg,vol,'uint16');
fclose(fimg);
k=findstr(imgfilename,'/');
if isempty(k)
  localfilename = imgfilename;
else
  localfilename = imgfilename(k(end)+1:end);
end;

% write header file
hdrname=sprintf('%s.dpf',fstem);
fhdr=fopen(hdrname,'w');
fprintf(fhdr,'Begin:\n');
fprintf(fhdr,'        ImageWidth: %d\n', size(vol,2));
fprintf(fhdr,'        ImageHeight: %d\n', size(vol,1));
fprintf(fhdr,'        ImageSlices: %d\n\n', size(vol,3));
fprintf(fhdr,'        ProcSliceStart: 0\n');
fprintf(fhdr,'        ProcSliceEnd: %d\n\n', size(vol,3)-1);
fprintf(fhdr,'        FieldOfView(X): %.1f\n', readoutFOV);
fprintf(fhdr,'        FieldOfView(Y): %.1f\n', phaseFOV);
fprintf(fhdr,'        SliceThickness: %.1f\n', SliceThickness);
fprintf(fhdr,'        SliceOrientation: %d\n', sliceorient);
fprintf(fhdr,'        SliceSequencing: %d\n', sliceseq);
fprintf(fhdr,'        B-Value: %d\n\n', bval);
if ~isempty(qmat)
  if size(qmat,2)~=3
    fprintf(fhdr,'        2nd dim of qmat should be 3 (it is %d)\n',size(qmat,2));
  else
    fprintf(fhdr,'        Begin_Of_Gradient_Table:\n');
    for q=1:size(qmat,1)
      fprintf(fhdr,'                 %d: %0.6f, %0.6f, %0.6f\n',...
        q-1,qmat(q,1),qmat(q,2),qmat(q,3));
    end;
    fprintf(fhdr,'        End_Of_Gradient_Table:\n\n');
  end;
end;
fprintf(fhdr,'        InputImgFile: %s\n', localfilename);
fprintf(fhdr,'        InputImgOrder:  Gradient_By_Gradient\n');
fprintf(fhdr,'End:\n\n');
fprintf(fhdr,'*** End of File ***\n');
fclose(fhdr);

