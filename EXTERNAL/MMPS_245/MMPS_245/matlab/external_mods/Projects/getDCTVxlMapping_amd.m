function vxlmap = getDCTVxlMapping_amd(vol, regStruct)

% Create voxel mapping from atlas to subject
minI = regStruct.range(1,1);
maxI = regStruct.range(1,2);
minJ = regStruct.range(2,1);
maxJ = regStruct.range(2,2);
minK = regStruct.range(3,1);
maxK = regStruct.range(3,2);
sizeatl = size(regStruct.VL.imgs);
atl_Mvxl2lph = regStruct.VL.Mvxl2lph;

Mvxl2vxl_atlaf2subj = inv(vol.Mvxl2lph) * (regStruct.M_atl_to_vol_af)*atl_Mvxl2lph;
Tr = regStruct.Tr;

count=0;
dim=[maxI-minI+1  maxJ-minJ+1 maxK-minK+1];
stablize=1;
[iK jK kK n]=size(Tr);
bi = stablize*dctb(dim(1), iK, (0:dim(1)-1)',0);
bj = stablize*dctb(dim(2), jK, (0:dim(2)-1)',0);
bk = stablize*dctb(dim(3), kK, (0:dim(3)-1)',0);
atlK=minK:maxK;
[atlI atlJ]=ndgrid(minI:maxI, minJ:maxJ);
VIimgs = zeros(sizeatl, 'single');
VJimgs = zeros(sizeatl, 'single');
VKimgs = zeros(sizeatl, 'single');
for l=1:dim(3)   % Cycle over planes
           % Nonlinear deformations
           %----------------------------------------------------------------------------
           ti = get_2Dtrans(Tr(:,:,:,1),bk,l);
           tj = get_2Dtrans(Tr(:,:,:,2),bk,l);
           tk = get_2Dtrans(Tr(:,:,:,3),bk,l);
           atlI_af = atlI    + bi*ti*bj';
           atlJ_af = atlJ    + bi*tj*bj';
           atlK_af = atlK(l)    + bi*tk*bj';

           [vi vj vk] = mmult(atlI_af, atlJ_af, atlK_af, Mvxl2vxl_atlaf2subj);
           VIimgs(minI:maxI,minJ:maxJ,atlK(l))=vi;
           VJimgs(minI:maxI,minJ:maxJ,atlK(l))=vj;
           VKimgs(minI:maxI,minJ:maxJ,atlK(l))=vk;
end       

ind  = find( abs(VIimgs)>eps | abs(VJimgs)>eps | abs(VKimgs)>eps);
vxlmap.atl2subj.indVol.imgs = zeros(sizeatl, 'uint32');
vxlmap.atl2subj.indVol.Mvxl2lph = atl_Mvxl2lph;
vxlmap.atl2subj.indVol.imgs(ind) = uint32(1:length(ind));
vxlmap.atl2subj.I = single(VIimgs(ind));
vxlmap.atl2subj.J = single(VJimgs(ind));
vxlmap.atl2subj.K = single(VKimgs(ind));

clear VIimgs VJimgs Vkimgs ti tj tk atlJ_af atlI_af atlK_af atlI atlJ atlK vi vj vk;

% Create voxel mapping from subject to atlas
%ind = find(vol.imgs> sqrt(eps));
ind = find(vol.imgs>=0); % Include all voxels (AMD added 8/28/06 -- old version masks out good parts of volume)
[I J K] = ind2sub(size(vol.imgs), ind);
vxlmap.subj2atl.indVol.imgs = zeros(size(vol.imgs), 'uint32');
vxlmap.subj2atl.indVol.Mvxl2lph = vol.Mvxl2lph;
vxlmap.subj2atl.indVol.imgs(ind) = uint32(1:length(ind));

np =length(ind);
vxl_subj=ones(4, np);
vxl_subj(1,:) = I';
vxl_subj(2,:) = J';
vxl_subj(3,:) = K';

vxl_atl_af = inv(Mvxl2vxl_atlaf2subj)* vxl_subj;
nK=[iK jK kK];
% Need to shift around the range box
[iI iJ iK] = getInvDCTMapMEX(np, vxl_atl_af(1,:)-minI+1, vxl_atl_af(2,:)-minJ+1, vxl_atl_af(3,:)-minK+1, ...
             squeeze(Tr(:,:,:,1)), squeeze(Tr(:,:,:,2)), squeeze(Tr(:,:,:,3)), dim, nK);
 
vxlmap.subj2atl.I = single(iI+minI-1);
vxlmap.subj2atl.J = single(iJ+minJ-1);
vxlmap.subj2atl.K = single(iK+minK-1);
