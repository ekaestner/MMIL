function coords_subj = Transform_Coords_Atlas_to_Subject(coords_atl,regStruct)

tmp = zeros(size(coords_atl,1),4);
tmp(:,1:3) = coords_atl;
tmp(:,4)=1;
% transform RAS to LPH
tmp = (M_RAS_TO_LPH*tmp')';
dl = vol_getvxlsval(tmp, regStruct.VL, eye(4,4), 2);
dp = vol_getvxlsval(tmp, regStruct.VP, eye(4,4), 2);
dh = vol_getvxlsval(tmp, regStruct.VH, eye(4,4), 2);
tmp(:,1) = tmp(:,1)+dl;
tmp(:,2) = tmp(:,2)+dp;
tmp(:,3) = tmp(:,3)+dh;
% transform LPH to RAS
tmp = (M_LPH_TO_RAS*tmp')';
coords_subj = tmp(:,1:3);
