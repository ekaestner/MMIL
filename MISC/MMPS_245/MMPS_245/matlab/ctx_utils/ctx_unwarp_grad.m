function ctx_voluw = ctx_unwarp_grad(ctx_vol,gwtype,unwarpflag,isocntrflag,interpm)
 
Mvxl2lph = ctx_vol.Mvxl2lph;
if isocntrflag
  firstLPH = Mvxl2lph*[1; 1; 1; 1];
  lastLPH = Mvxl2lph*[ctx_vol.dimr; ctx_vol.dimc; ctx_vol.dimd; 1];
  ctx_vol.Mvxl2lph(3,4) = ctx_vol.Mvxl2lph(3,4)-(firstLPH(3)+lastLPH(3))/2;
end
if exist('interpm','var')
  ctx_voluw = vol_unwarp_grad(ctx_vol, gwtype, 1, 1, unwarpflag,interpm);
else
  ctx_voluw = vol_unwarp_grad(ctx_vol, gwtype, 1, 1, unwarpflag);
end
if isocntrflag
  ctx_vol.Mvxl2lph = Mvxl2lph;
  ctx_voluw.Mvxl2lph = Mvxl2lph;
end

