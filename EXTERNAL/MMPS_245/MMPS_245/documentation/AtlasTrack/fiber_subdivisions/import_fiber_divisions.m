indir = '/space/md18/2/data/mmildev/cjpung/fiber_segmentation';
outdir = '/space/md5/2/data/dhagler/extra/papers/RSI_tracto/fiber_subdivisions';
fname_legend = [outdir '/DTI_Fiber_Subdivisions_Legend.csv'];
suffix = 'mask';
thresh = eps;
forceflag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FiberLegend = DTI_MMIL_Read_Fiber_Legend(fname_legend);

for f=1:length(FiberLegend)
  fnum = FiberLegend(f).FiberNumber;
  fname_in = [indir '/' FiberLegend(f).FileName];
  fname_out = sprintf('%s/fiber_%04d_%s.mat',outdir,fnum,suffix);
  if ~exist(fname_out,'file') || forceflag
    [vol,M] = fs_load_mgh(fname_in);
    vol = 1.0*(vol > thresh);
    mmil_save_sparse(vol,fname_out,M);
  end;
end;
  
