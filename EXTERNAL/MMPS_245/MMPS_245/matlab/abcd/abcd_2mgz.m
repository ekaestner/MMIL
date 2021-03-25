function abcd_2mgz(imgpath,mode,orientation,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function abcd_2mgz(imgpath,mode,varargin)
% Purpose
%  1. convert dcm to nii using dcm2niix
%  2. convert nii to mgh using fs_mri_convert
%  3. make sure orientation is correct
%  4. save as mgz
% Required:
%  'imgpath': where dcm/nii/mgh are saved
%     {default = []}
%  'mode': [1|2|3] convertion mode
%     1 = starting from dcm
%     2 = starting from nii
%     3 = starting from mgh
%     4 = starting from mgz
%     {default = []}
%
% Created:  07/18/17 by Feng Xue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes:
%for GE:
%  T1/T2(tested on G031/G032/G087_INV1JAFHY62/G087_INVZYFB43H7/G075_INVZFYGCHKB/G010_INVYTH97JCP/G010_INVZZZ2ALR6/G054_INVWFTA7NUU/G054_INVYHMMV4AC/G075_INVZMWM1Y06): 
%    dcm2niix
%    nii=load_nifti
%    reorient if necessary
%    fs_save_mgh(nii.vol,'T1.mgz',nii.sform) (for M matrix, can use one from nii or dcm even though they are not equal, use M_new = M_dcm; M_new(1,3)=M_nii(1,3))
%   DTI(tested on G031 and G032/G087_INV1JAFHY62/G087_INVZYFB43H7/G075_INVZFYGCHKB/G010_INVYTH97JCP/G010_INVZZZ2ALR6/G054_INVWFTA7NUU/G054_INVYHMMV4AC/G075_INVZMWM1Y06)/fm(tested on G031 and G032/G087_INV1JAFHY62/G087_INVZYFB43H7/G075_INVZFYGCHKB/G010_INV0F78WV5U/G054_INVWFTA7NUU/G054_INVYHMMV4AC/G075_INVZMWM1Y06)/fmri(tested on G032/G031/G087_INV1JAFHY62/G087_INVZYFB43H7/G075_INVZFYGCHKB/G010_INVYTH97JCP/G010_INVZZZ2ALR6/G054_INVWFTA7NUU/G054_INVYHMMV4AC/G075_INVZMWM1Y06): 
%      # images in a frame = size(vol,3)
%      M_test=mmil_read_dicom_M({a(3:size(vol,3)-1).name})
%    fmri_fm(tested on G031/G032/G087_INV1JAFHY62/G087_INVZYFB43H7/G075_INVZFYGCHKB/G010_INVYTH97JCP/G054_INVWFTA7NUU/failed on G010_INVZZZ2ALR6/G054_INVYHMMV4AC/G075_INVZMWM1Y06):
%      # images in a frame = size(vol,3)
%      M_test=mmil_read_dicom_M({a(3:size(vol,3)+3-1).name})
%       vol(:,:,:,2) = for
%       flip(vol(:,:,:,1),2) = rev 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  path_curr = pwd;
  imgname = [];
  switch exist(imgpath)
    %regular file
    case 2 
      [imgpath,imgname,imgext] = fileparts(imgpath);
      if isempty(imgpath)
        imgpath = pwd;
      else
        cd(imgpath);
        imgpath = pwd;
      end

      switch lower(imgext)
        case 'nii'
          mode = 2;
        case 'gz'
          mode = 2;
        case 'dcm'
          mode = 1;
        case 'mgh'
          mode = 3;
        case 'mgz'
          mode = 4;
      end
    %directory
    case 7
      cd(imgpath);
      imgpath = pwd;
    otherwise
      error(sprintf('imgpath: %s does not exist',imgpath));
  end

  switch mode
    case 1
      %Could place a flag here to ignore derived images.
      if isempty(imgname)
        cmd=sprintf('dcm2niix -z n -p n -f %s -o %s %s','%n_%k_%j_%t_%s',imgpath,imgpath);
        unix(cmd);
        imglist = dir('*.nii');
      else
        cmd=sprintf('dcm2niix -s y -z n -p n -f %s -o %s %s %s%s',imgname,imgpath,imgpath,imgname,imgext);
        unix(cmd);
        imglist = dir(sprintf('%s.nii',filename));
      end
    case 2
      if isempty(imgname)
        imglist = dir('*.nii');
      else
        imglist = dir(sprintf('%s%s',imgname,imgext));
      end
    case 3
      if isempty(imgname)
        imglist = dir('*.mgh');
      else
        imglist = dir(sprintf('%s%s',imgname,imgext));
      end
    case 4
      if isempty(imgname)
        imglist = dir('*.mgz');
      else
        imglist = dir(sprintf('%s%s',imgname,imgext));
      end
    otherwise
      error('Wrong mode, only 1-4 are accepted');
  end
  if ~isempty(imglist)
    if isempty(orientation), orientation = 'LPS'; end;
    for idx=1:length(imglist)
      [~,filestem] = fileparts(imglist(idx).name);
      if mode<3, fs_mri_convert(sprintf('%s/%s',imgpath,imglist(idx).name),sprintf('%s/%s.mgh',imgpath,filestem)); end
      %%%%%%%%%%%%%%%%%%%%%%%
      %Consider use fs_load_nii etc.
      %%%%%%%%%%%%%%%%%%%%%%%
      if mode == 4
        [vol,M] = fs_load_mgh(sprintf('%s/%s.mgz',imgpath,filestem));
      else
        [vol,M] = fs_load_mgh(sprintf('%s/%s.mgh',imgpath,filestem));
      end
      orientation_curr = fs_read_orient('',M);
      try
        %image must be full image (not partial slices) or fs_reorient will fail.
        if ~strcmp(orientation_curr,orientation), [vol,M] = fs_reorient(vol,M,orientation); end;
        fs_save_mgh(vol,sprintf('%s/%s.mgz',imgpath,filestem),M);
      catch ME
        warning(ME.message);
      end
      if mode < 4
        %Could place a flag here to decide whether to delete mgh file
        delete(sprintf('%s/%s.mgh',imgpath,filestem));
      end
    end
  else
    error('No image exist');
  end
  
  %cd back 
  cd(path_curr);
end
