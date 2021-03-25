function [vol,M,volsz] = ejk_load_vol(fname)

  vol = []; M = []; volsz = [];
  [~,~,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [vol,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [vol,M,~,volsz] = fs_load_mgh(fname);
  end
  
return;