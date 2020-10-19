function ADNI_MMIL_Build_PET_hist()
%function ADNI_MMIL_Build_PET_hist()
%
% NOTE: this function is out of date and will not work as is
%
% Last Mod:  03/24/10 by Don Hagler
%

ContainerRootDir = '/space/invaders/1/data/MMILDB/ADNI/Containers';

h = {};
dirlist = dir(sprintf('%s/PETRAW*',ContainerRootDir));
for i = 1:length(dirlist)
  ContainerDir = char(dirlist(i).name);
  fprintf(1,'i=%d of %d: %s\n',i,length(dirlist),ContainerDir);
  fname = sprintf('%s/%s/vols.mat',ContainerRootDir,ContainerDir);
  if exist(fname,'file')
    load(fname);
    if length(vols) == 6
      vol = 0;
      for j = 1:length(vols)
        vol = vol+max(0,double(vols{j}.vol));
      end
      vol = vol/length(vols); 
      [hc,bv] = hist(vol(:),1000);
      h{length(h)+1} = struct('hc',hc,'bv',bv);
    end
  end
  figure(1); clf; plot(bv,cumsum(hc)/sum(hc)); axis tight;
  crange = [0 max(vol(:))];
  dims = size(vol);
  figure(2); clf;
  subplot(2,2,1); imagesc(squeeze(vol(floor(dims(1)/2),:,:))',crange); colormap(hot); axis image;
  subplot(2,2,2); imagesc(squeeze(vol(:,floor(dims(2)/2),:))',crange); colormap(hot); axis image;
  subplot(2,2,3); imagesc(squeeze(vol(:,:,floor(dims(3)/2))),crange); colormap(hot); axis image;
  drawnow; pause(0.1);
end

