function PITU_Modify_Brainmask(dirname,forceflag)
%function PITU_Modify_Brainmask(dirname,forceflag)
%
% Rcnt mod: 03/05/10 by Don Hagler
% Last Mod: 05/21/12 by Don Hagler
%

atlas_dir = [getenv('MMPS_DIR') '/atlases'];

load [atlas_dir '/T1_Atlas/T1_atlas'];
[maskedvol, volmask] = getmaskvol(volm, bm, eye(4));
showVol(volm,volmask);

volmask.imgs = volmask.imgs*128;
ctx_save_mgh(volm,'/home/dale/volm.mgh');
ctx_save_mgh(volstd,'/home/dale/volstd.mgh');
ctx_save_mgh(volmask,'/home/dale/volmask.mgh');

volmask = ctx_load_mgh('/home/dale/matlab/atlases/volmask.mgh'); % Read in edited volmask
volmask.imgs = 1.0*(volmask.imgs>5); volmask.maxI = 1.0;
%volmask.imgs = mmil_smooth3d(volmask.imgs,10,10,10);

save([atlas_dir '/PITU/pitu_afm.mat'],'volm','volstd','volmask');

if 0
FV = isosurface(volmask.imgs,0.5);
FV = reducepatch(FV, 5000);
%icoords = [FV.vertices ones(size(FV.vertices,1),1)];
icoords = [FV.vertices(:,[2 1 3]) ones(size(FV.vertices,1),1)];
rcoords = (volmask.Mvxl2lph*icoords')';
bm2 = FV;
bm2.vertices = rcoords(:,1:3);
bm2 = preprocessQ(bm2);
[maskedvol, volmask2] = getmaskvol(volm, bm2, eye(4)); % Coredumps for some reason
end

