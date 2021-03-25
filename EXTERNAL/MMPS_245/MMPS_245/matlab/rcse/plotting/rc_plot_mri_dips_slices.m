function rc_plot_mri_dips_slices(subject,prefix,slicerange,plane,inplaneval,outdir)
%function rc_plot_mri_dips_slices(subject,prefix,[slicerange],[plane],[inplaneval],[outdir])
%
% Purpose: plot MRI and RCSE dipoles for multiple slices
%
% Required Input:
%   subject: FreeSurfer subject name
%   prefix: RCSE prefix
%
% Optional Input:
%   slicerange: vector with 2 numbers between 1 and 256 (slicenums)
%   plane: 'cor', 'sag', or 'hor'
%     {default = 'cor'}
%   inplaneval: 0, 1,or greater
%     if 0, plot all dips on all slices in slicerange
%     if 1, only plot dips in the mri plane for each slice in slicerange
%     if >=2, plot dips on closest of inplaneval # of slices in slicerange
%    {default: 0}
%  outdir: output directory name
%    {default = []}
%
% Early Mod: 12/18/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('slicerange','var') | isempty(slicerange), slicerange = [1 256]; end;
if ~exist('plane','var') | isempty(plane), plane = 'cor'; end;
if ~ismember(plane,{'cor','sag','hor'})
  error('plane must be cor, sag, or hor');
end;
if ~exist('inplaneval','var'), inplaneval = []; end;
if isempty(inplaneval), inplaneval = 0; end;
if ~exist('outdir','var'), outdir = []; end;

label_flag = 1;

if length(slicerange)==1
  slicerange = [slicerange slicerange];
end;
for s=1:2
  slicenum = slicerange(s);
  if slicenum<1 | slicenum>256
    error('slicenum must be between 1 and 256');
  end;
end;

if isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir')
  [success,msg] = mkdir(outdir);
  if ~success
    error('failed to create dir %s:\n%s',outdir,msg);
  end;
end;

if inplaneval < 2
  for slicenum=slicerange(1):1:slicerange(2)
    clf;
    rc_plot_mri(subject,slicenum,plane);
    if inplaneval
      rc_plot_dips(prefix,slicenum,plane,[],label_flag);
      if label_flag
        outname = sprintf('%s/mri_dips_%s_slice%d_wlabel.tif',outdir,plane,slicenum);
      else
        outname = sprintf('%s/mri_dips_%s_slice%d.tif',outdir,plane,slicenum);
      end;
    else
      rc_plot_dips(prefix,0,plane,[],label_flag);
      if label_flag
        outname = sprintf('%s/mri_dips_%s_slice%d_alldips_wlabel.tif',outdir,plane,slicenum);
      else
        outname = sprintf('%s/mri_dips_%s_slice%d_alldips.tif',outdir,plane,slicenum);
      end;
    end;  
    print('-dtiff',outname);
  end;
else
  slicelist = [slicerange(1):slicerange(2)];
%fprintf('%s: slicelist:\n',mfilename);     
%disp(slicelist);
  nslices = slicerange(2)-slicerange(1)+1;
  nslabs = inplaneval;
  nslabslices = nslices/nslabs;
  plotslicelist = [];
  for i=1:nslabs
    tmpslicelist = round((i-1)*nslabslices) + [1:ceil(nslabslices)];
    midslice = floor(mean(tmpslicelist));
    plotslicelist = [plotslicelist midslice];
  end;
  plotslicelist = slicelist(plotslicelist);
  for i=1:length(plotslicelist)
    slicenum = plotslicelist(i);
    dipslicelist = round([slicenum-nslabslices/2:1:slicenum+nslabslices/2]);
%  fprintf('%s: plotslice = %d, dipslicelist:\n',mfilename,slicenum);     
%  disp(dipslicelist);
    clf;
    rc_plot_mri(subject,slicenum,plane);
    rc_plot_dips(prefix,dipslicelist,plane,[],label_flag);
    if label_flag
      outname = sprintf('%s/mri_dips_%s_slice%d_inplane%d_wlabel.tif',...
        outdir,plane,slicenum,inplaneval);
    else
      outname = sprintf('%s/mri_dips_%s_slice%d_inplane%d.tif',...
        outdir,plane,slicenum,inplaneval);
    end;
    print('-dtiff',outname);
  end;
end;
