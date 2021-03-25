function [vol_sze,M] = ejk_read_vol_sze(fname,parms)

  vol_sze = []; M = [];
  [~,~,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [~,M,vol_sze] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [M,vol_sze] = mmil_load_mgh_info(fname,0,[]);
  end
  vol_sze = vol_sze(1:3);
  
return;