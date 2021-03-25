function [Mvox2ras] = compute_vox2ras(mdhStruct)

% Assemble vox2ras matrix straight from header (not relying on Rudolph's pre-computed)

nslices = get_mdh_tagval(mdhStruct,'sSliceArray.lSize');
npartitions = get_mdh_tagval(mdhStruct,'sKSpace.lPartitions');
baseres = get_mdh_tagval(mdhStruct,'sKSpace.lBaseResolution');
PhaseFOV = get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].dPhaseFOV');
ReadoutFOV = get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].dReadoutFOV');
thickness = get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].dThickness');
slicepos1 = [get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].sPosition.dSag');...
             get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].sPosition.dCor');...
             get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[0].sPosition.dTra')];

if isfield(mdhStruct,'m_tMRAcquisitionType')
  is3dflag = (strcmp(mdhStruct.m_tMRAcquisitionType,'3D'));
else
  is3dflag = (strcmp(mdhStruct.tMRAcquisitionType,'3D'));
end

if is3dflag
  voldim = [baseres baseres npartitions];
  voxsize = [PhaseFOV/voldim(1) ReadoutFOV/voldim(2) thickness/voldim(3)]; 
else
  voldim = [baseres baseres nslices];
  slicepos2 = [get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[1].sPosition.dSag');...
               get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[1].sPosition.dCor');...
               get_mdh_tagval(mdhStruct,'sSliceArray.asSlice[1].sPosition.dTra')];
  
  voxsize = [PhaseFOV/voldim(1) ReadoutFOV/voldim(2) norm(slicepos2-slicepos1)]; 
end

M_voxScale 	= eye(3);
M_voxScale(1,1)	= voxsize(1);
M_voxScale(2,2)	= voxsize(2);
M_voxScale(3,3)	= voxsize(3);

M_Wn = mdhStruct.adRM';
[r c]	= size(M_Wn);
if r~=0
  M_mask1	= [1 0 0; 0 1 0; 0 0 -1];
  M_mask2	= [-1 0 0; 0 1 0; 0 0 -1];
  M = M_mask1*M_Wn*M_voxScale*M_mask2;
  slicepos1ras = [-slicepos1(1); -slicepos1(2); slicepos1(3)];
  if is3dflag
    transvec = slicepos1ras-M*[(voldim(1)+1)/2; (voldim(2)+1)/2; (voldim(3)+1)/2];
  else
    transvec = slicepos1ras-M*[(voldim(1)+1)/2; (voldim(2)+1)/2; 1.5];
  end
  Mvox2ras = eye(4);
  Mvox2ras(1:3,1:3) = M;
  Mvox2ras(1:3,4) = transvec;
end
