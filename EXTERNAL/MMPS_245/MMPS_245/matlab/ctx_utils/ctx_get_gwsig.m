function gwsig = ctx_get_gwsig(vol)
%function gwsig = ctx_get_gwsig(vol)
%
% Purpose: determine gradwarp signature from an image
%
% Required Input:
%   vol: input volume (ctx format)
%
% Output:
%   gwsig: struct with these fields
%     valmat ?
%     z1 ?
%     z2 ?
% 
% Created:  12/11/10 by Don Hagler
% Last Mod: 12/11/10 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
gwsig = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mvxl2lph = vol.Mvxl2lph;
dims = size(vol.imgs);
ctr_ijk = (1+dims')/2;
ctr_lph = Mvxl2lph(1:3,:)*[ctr_ijk; 1];
ctr_lph(3) = 0; % Adjust for iso-center scanning
Mvxl2lph(1:3,4) = Mvxl2lph(1:3,4)+(ctr_lph-Mvxl2lph(1:3,:)*[ctr_ijk; 1]);
Mlph2vxl = inv(Mvxl2lph);
ctr_ijk = Mlph2vxl(1:3,:)*[ctr_lph; 1];
v1 = Mvxl2lph(1:3,1); v1 = v1/norm(v1);
v2 = Mvxl2lph(1:3,2); v2 = v2/norm(v2);
v3 = Mvxl2lph(1:3,3); v3 = v3/norm(v3);
vz = v1(3)*v1 + v2(3)*v2;
vt = cross(vz,v3);
zvals = [-50:0.5:50];
tvals = [-200:0.5:200];
svals = [0];
[zmat,tmat,smat] = ndgrid(zvals,tvals,svals);
lphmat = (ones(3,1)*zmat(:)').*(vz*ones(1,length(zmat(:))))+(ones(3,1)*tmat(:)').*(vt*ones(1,length(tmat(:))))+(ones(3,1)*smat(:)').*(v3*ones(1,length(smat(:))));
lphmat = [lphmat; ones(1,size(lphmat,2))];
vxlmat = Mlph2vxl(1:3,:)*lphmat;
valvec = interp3(vol.imgs,vxlmat(2,:),vxlmat(1,:),vxlmat(3,:),'nearest');
valmat = reshape(valvec,size(zmat));
z1 = zeros(size(valmat,1),1);
z2 = zeros(size(valmat,1),1);
z1(:) = NaN;
z2(:) = NaN;
for j = 1:size(valmat,1)
  indx = find(isfinite(valmat(j,:)));
  if ~isempty(indx)
    indx2 = find(valmat(j,:)>0 & isfinite(valmat(j,:)));
    if indx2(1)>indx(1) z1(j) = (indx2(1)-indx(1)); end
    if indx2(end)<indx(end) z2(j) = (indx(end)-indx2(end)); end
  end
end

gwsig.valmat = valmat;
gwsig.z1 = z1;
gwsig.z2 = z2;
