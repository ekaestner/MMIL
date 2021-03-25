function rc_create_ring_weights(subjname,areas,centers,fwhm_inner,fwhm_outer,...
                             fwhm_post,outdir,outstem,surfname);
%function rc_create_ring_weights(subjname,areas,centers,fwhm_inner,fwhm_outer,...
%                             [fwhm_post],[outdir],[outstem],[surfname]);
%
% create w files from center vertices
%
%
% Created: 11/16/06 Don Hagler
% Lst Mod: 02/19/11 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargin(nargin,3), return; end;
if ~exist('fwhm_post','var') || isempty(fwhm_post), fwhm_post = 2; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('outstem','var') || isempty(outstem), outstem = 'retweights'; end;
if ~exist('surfname','var') || isempty(surfname), surfname = 'white'; end;

initval = 10;

sm_inner = round((fwhm_inner/1.25).^2); % FWHM ~ 1.25*sqrt(N)
sm_outer = round((fwhm_outer/1.25).^2);
sm_post =  round((fwhm_post/1.25).^2);

hemilist = {'lh','rh'};
nhemi = length(hemilist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for errors
ncenters = length(centers(1).lh);
for j=1:length(areas)
  if length(centers(j).lh) ~= ncenters | ...
     length(centers(j).rh) ~= ncenters
    fprintf('%s: error: areas should each have the same number of dipoles!\n',...
      mfilename);
    return;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!mkdir -p matfiles
if ~exist(outdir,'dir')
  [success,message]=mkdir(outdir);
  if ~success
    fprintf('%s: unable to create output dir %s... quitting\n',mfilename,outdir);
    return;
  end;
end;

for h=1:nhemi
  hemi = hemilist{h};

  % load subject's surface
  matname=sprintf('matfiles/%s-surf-%s-%s.mat',subjname,surfname,hemi);
  if exist(matname,'file')
    load(matname);
  else
    subjdir = getenv('SUBJECTS_DIR');
    if isempty(subjdir)
      fprintf('%s: SUBJECTS_DIR not defined as an environment variable... quitting\n',mfilename);
      return;
    end;
    surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subjname,hemi,surfname);
    if ~exist(surffile,'file')
      fprintf('%s: surface file %s not found... quitting\n',mfilename,surffile);
      return;
    end
    surf = fs_read_surf(surffile);
    surf = fs_find_neighbors(surf);

    % save to mat file
    save(matname,'surf');
  end;

  % for each center vertex, smooth, make ring, smooth, then save as w file
  fprintf('%s: creating w files\n',mfilename);
  for a=1:length(areas)
    for c=1:ncenters
      v = getfield(centers,{a},hemi,{c});
      w = zeros(surf.nverts,1);
      w(v)=initval;
      w_inner=fs_smooth(surf,w,sm_inner);
      v_inner=find(w_inner);
      w=w_inner(v_inner);
      fname = sprintf('%s/%s-%s-center%d-fwhm_in%d-fwhm_out%d-fwhm_post%d-%s.w',...
        outdir,outstem,areas(a).name,c-1,fwhm_inner,fwhm_outer,fwhm_post,hemi);
      [status]=fs_write_wfile(fname,w,v_inner); if ~status, return; end;

      v = getfield(centers,{a},hemi,{c});
      w = zeros(surf.nverts,1);
      w(v)=initval;
      w_outer=fs_smooth(surf,w,sm_outer);
      v_outer=find(w_outer);
      w_inner=w_inner(v_inner);
      w_outer=w_outer(v_outer);
      w = zeros(surf.nverts,1);
      w_inner=rc_norm_weights(w_inner,1);  % norm to max
      [w_inner,v_inner]=rc_thresh_weights(w_inner,v_inner,0.5); % thresh relative to max
      w(v_outer) = w_outer;
      w(v_inner) = 0;
      w=fs_smooth(surf,w,sm_post);
      v=find(w);
      w=w(v);
      fname = sprintf('%s/%s-%s-ring%d-fwhm_in%d-fwhm_out%d-fwhm_post%d-%s.w',...
        outdir,outstem,areas(a).name,c-1,fwhm_inner,fwhm_outer,fwhm_post,hemi);
      [status]=fs_write_wfile(fname,w,v); if ~status, return; end;

    end;
  end;
end;

return;
