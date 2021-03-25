function ADNI_MMIL_Build_MRI_hist()
%function ADNI_MMIL_Build_MRI_hist()
%
% Last Mod: 09/14/12 by Don Hagler
%

ContainerRootDir = '/space/invaders/1/data/MMILDB/ADNI/Containers';

h = {};
dirlist = dir(sprintf('%s/MRIPROC*',ContainerRootDir));
for i = 1:length(dirlist)
  ContainerDir = char(dirlist(i).name);
  fprintf(1,'i=%d of %d: %s\n',i,length(dirlist),ContainerDir);
  MPRfname = sprintf('%s/%s/MPR_res.mgz',ContainerRootDir,ContainerDir);
  if exist(MPRfname,'file')
    vol = ctx_load_mgh(MPRfname);
    [hc,bv] = hist(vol.imgs(:),1000);
    h{length(h)+1} = struct('hc',hc,'bv',bv);
  end
  figure(1); plot(bv,cumsum(hc)/sum(hc)); axis tight;
  drawnow; pause(0.1);
end

