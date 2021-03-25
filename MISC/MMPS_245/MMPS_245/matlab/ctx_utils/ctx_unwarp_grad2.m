function ctx_voluw = ctx_unwarp_grad(ctx_vol,gwtype,unwarpflag,isocntrflag,jacobianflag)

if ~exist('jacobianflag','var')
  jacobianflag = 1;
end
 
Mvxl2lph = ctx_vol.Mvxl2lph;
if isocntrflag
  firstLPH = Mvxl2lph*[1; 1; 1; 1];
  lastLPH = Mvxl2lph*[ctx_vol.dimr; ctx_vol.dimc; ctx_vol.dimd; 1];
  ctx_vol.Mvxl2lph(3,4) = ctx_vol.Mvxl2lph(3,4)-(firstLPH(3)+lastLPH(3))/2;
end
ctx_voluw = vol_unwarp_grad(ctx_vol, gwtype, jacobianflag, 1, unwarpflag);
if isocntrflag
  ctx_vol.Mvxl2lph = Mvxl2lph;
  ctx_voluw.Mvxl2lph = Mvxl2lph;
end

