function ocode = ctx_get_image_orientation(ctx_vol_in)

slice_normal = ctx_vol_in.Mvxl2lph(1:3,3);
slice_normal = slice_normal/norm(slice_normal);
sag_corr = dot(slice_normal,[1; 0; 0]);
cor_corr = dot(slice_normal,[0; -1; 0]);
hor_corr = dot(slice_normal,[0; 0; 1]);

if abs(sag_corr)>=0.5
  ocode = 2;
elseif abs(cor_corr)>=0.5
  ocode = 1;
else
  ocode = 0;
end
