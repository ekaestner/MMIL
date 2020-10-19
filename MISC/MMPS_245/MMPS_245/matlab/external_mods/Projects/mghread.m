function [mghhr,mghimg]=mghread(file)
% mghread: Read mgh files
% [mghhr,mghimg]=mghread(file)
% file: input mgh file name
% mghhr: header structure
% mghimg: 3-D images
% Gen-Nan Chen, CorTechs Labs.

%[pathstr, fname, fext, versn] = fileparts(file);
[pathstr, fname, fext] = fileparts(file); % DH 07/20/2017

bzip =false;
if strcmp(fext,'.mgz')
    str = sprintf('!gzip -S .mgz -d %s -c > /tmp/tmp.mgh', file);
    eval(str);
    file = '/tmp/tmp.mgh';
    bzip = true;
end

fid = fopen(file,'r', 'b');
mghhr.ver=fread(fid,1,'int32');
mghhr.width=fread(fid,1,'int32');
mghhr.height=fread(fid,1,'int32');
mghhr.depth=fread(fid,1,'int32');
mghhr.nframes=fread(fid,1,'int32');
mghhr.type=fread(fid,1,'int32');
mghhr.dof=fread(fid,1,'int32');
mghhr.ras_good_flag=fread(fid,1,'int16');
mghhr.xsize=fread(fid,1,'float32');
mghhr.ysize=fread(fid,1,'float32');
mghhr.zsize=fread(fid,1,'float32');
mghhr.x_r=fread(fid,1,'float32');
mghhr.x_a=fread(fid,1,'float32');
mghhr.x_s=fread(fid,1,'float32');
mghhr.y_r=fread(fid,1,'float32');
mghhr.y_a=fread(fid,1,'float32');
mghhr.y_s=fread(fid,1,'float32');
mghhr.z_r=fread(fid,1,'float32');
mghhr.z_a=fread(fid,1,'float32');
mghhr.z_s=fread(fid,1,'float32');
mghhr.c_r=fread(fid,1,'float32');
mghhr.c_a=fread(fid,1,'float32');
mghhr.c_s=fread(fid,1,'float32');
unused=194;
emp=fread(fid,unused,'uchar');
unknown=0;
switch mghhr.type
    case 0
        dt='uchar';
    case 1
        dt='int';
    case 2
        dt='int64';
    case 3
        dt='float';
    case 4
        dt='uint16';
    otherwise
        unknown=1;
        disp('Unknown datatype.')
end

if (~unknown)
  [data count]=fread(fid,mghhr.width*mghhr.height*mghhr.depth,dt);
  mghhr.maxI = max(data);
  mghhr.minI = min(data);
  imgt=reshape(data,mghhr.width,mghhr.height,mghhr.depth);
  %mghimg=reshape(data,mghhr.height,mghhr.width,mghhr.depth);
  mghimg=zeros(mghhr.height,mghhr.width,mghhr.depth);
  for i=1:mghhr.depth
    mghimg(:,:,i)=imgt(:,:,i)';
  end
else
  mghimg=0;
end
fclose(fid);

if (bzip)
    !rm -f /tmp/tmp.mgh
end
