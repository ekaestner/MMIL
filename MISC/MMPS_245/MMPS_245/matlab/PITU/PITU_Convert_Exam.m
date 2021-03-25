function PITU_Convert_Exam(dirname,forceflag)
%function PITU_Convert_Exam(dirname,forceflag)
%
% Created:  06/28/07 by Anders Dale
% Last Mod: 12/21/12 by Don Hagler
%

if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end

fprintf(1,'PITU_Convert_Exam(''%s'',%d)\n',dirname,forceflag);

donefname = sprintf('%s/done.log',dirname);
if exist(donefname,'file') & ~forceflag
  fprintf(1,'%s: %s already exists\n',mfilename,donefname);
  return
end

fieldnamelist = {'SequenceName','SeriesNumber','FlipAngle','InstanceNumber','EchoNumber','EchoTime'};

cd(dirname)
filelist = dir(sprintf('%s/*.IMA',dirname));
if isempty(filelist)
  filelist = dir(sprintf('%s/*.dcm',dirname));
end
fieldvals = {};
for i = 1:length(filelist)
  fname = char(filelist(i).name);
  fnames{i} = fname;
%  fprintf(1,'reading file %s: %d of %d\n',fname,i,length(filelist));
  hdr = dicominfo(fname);
  for j = 1:length(fieldnamelist)
    fieldvals{i,j} = mmil_getfield(hdr,fieldnamelist{j},0);
  end
end
indvec = strcmp('fl3d8_ns',{fieldvals{:,strmatch('SequenceName',fieldnamelist)}}); % Multi-echo images
SeriesNumbers = unique([fieldvals{indvec,strmatch('SeriesNumber',fieldnamelist)}]);
for serindx = 1:length(SeriesNumbers)
  indvec2 = ([fieldvals{:,strmatch('SeriesNumber',fieldnamelist)}]==SeriesNumbers(serindx));
  EchoNumbers = unique([fieldvals{indvec2,strmatch('EchoNumber',fieldnamelist)}]);
  EchoTimes = unique([fieldvals{indvec2,strmatch('EchoTime',fieldnamelist)}]);
  if length(find(indvec2))==1024
    volsum = 0;
    for echoindx = 1:length(EchoNumbers)
      indvec3 = indvec2 & ([fieldvals{:,strmatch('EchoNumber',fieldnamelist)}]==EchoNumbers(echoindx));
      indlist = find(indvec3);
      [sortvals,sortorder] = sort([fieldvals{indlist,strmatch('InstanceNumber',fieldnamelist)}]); 
      indlist = indlist(sortorder);
      [vol,M] = mmil_read_dicom_vol({fnames{indlist}});
      volsum = volsum + vol.^2;
    end
    vol = sqrt(volsum/length(EchoNumbers));
    ctx_vol = ctx_mgh2ctx(vol,M);
    [gradwarpinfo,errmsg] = ctx_get_gradwarpinfo(hdr);
    ctx_voluw = ctx_unwarp_grad(ctx_vol, gradwarpinfo.gwtype, gradwarpinfo.unwarpflag, gradwarpinfo.isoctrflag);
    fname = sprintf('%s/MEDIC%d.mgh',dirname,serindx);
    ctx_write2_mgh(ctx_voluw,fname,'MRI_FLOAT');
    fprintf(1,'file %s written (%d of %d)\n',fname,serindx,length(SeriesNumbers));
  else
    fprintf(1,'Series %d (%d of %d): incorrect number of images in series (%d ~= 1024))\n',SeriesNumbers(serindx),serindx,length(SeriesNumbers),length(find(indvec2)));
  end
end

fname_mask = sprintf('%s/MASK.mgh',dirname);
volmask = mmil_dct_brainmask(ctx_vol,'atlasname','PITU/pitu_afm','smooth',50);
ctx_save_mgh(volmask,fname_mask);

system(sprintf('echo "%s completed at %s" >> %s\n',mfilename,datestr(now),donefname));
