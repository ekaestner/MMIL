function rc_get_offsets_from_corr(varargin)
%function rc_get_offsets_from_corr[options])
%
% 'prefix_stem': stem of RCSE output prefix
%   {default = 'RCSE'}
% 'areas': vector of area numbers
%   If empty, assume idependent offset fits
%     for different visual areas were not done
%   {default = []}
% 'hemis':  vector of hemifield indices (1=right hemifield, 2=left hemifield)
%   If empty, assume idependent offset fits
%     for different hemifields were not done
%   {default = []}
% 'uplows':  vector of upper/lower field indices (1=upper field, 2=lower field)
%   If empty, assume idependent offset fits
%     for different fields were not done
%   {default = []}
% 'eccs': vector of eccentricity indices
%   If empty, assume idependent offset fits
%     for different eccentricities were not done
%   {default = []}
% 'contvec': vector of contrast indices used simultaneously for fits
%   If empty, assume all contrast levels were used
%   {default = []}
% 'fname_out': name of output file
%   If empty, write to stdout
%   {default = []}
% 'hemi_corr_flag': calculate correlation between hemispheres
%   and use as weighting for selection of dipoles
%   {default = 1}
% 'self_corr_flag': calculate correlation between offset grid locations
%    within a quarterfield
%   {default = 1}
% 'ref_corr_flag': calculate correlation with reference waveform
%   and use as weighting or mask for selection of dipoles
%   reference waveform assumes zero baseline, 1 between fit_time0 and fit_time1
%   {default = 1}
% 'ref_corr_pol': [-1,0,1] polarity of reference waveform
%   if -1, assumes negative polarity
%   if 0, allows either negative or positive polarity
%   if 1, asssumes positive polarity
%   {default = 0}
% 'ref_corr_thresh': correlation with reference threshold
%   0 keeps values, but must be positive
%   otherise, apply threshold and binarize
%   {default = 0}
% 'expvar_thresh': explained variance threshold
%   0 keeps values, but must be positive
%   otherise, apply threshold and binarize
%   {default = 0}
% 'expvar_fact': exponential weighting applied to explained variance in
%   determining which offsets to choose
%   {default = 0}
% 'plotflag': [0|1] whether to make colorful plots
%   {default = 1}
% 'plot_outdir': output directory for plot tif files
%   {default = 'plots'}
% 'forceflag': [0|1] whether to overwrite output file
%   {default = 1}
%
% Created:  04/10/09 by Don Hagler
% Last Mod: 02/19/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin, { ...
  'prefix_stem','RCSE',[],...
  'areas',[],[],...
  'hemis',[],[],...
  'uplows',[],[],...
  'eccs',[],[],...
  'contvec',[],[],...
  'fname_out',[],[],...
  'expvar_thresh',0,[0,1],...
  'ref_corr_thresh',0,[0,1],...
  'expvar_fact',5,[0,1000],...
  'expvar_clim',[0,0.4],[],...
  'forceflag',true,[false true],...
  'hemi_corr_flag',true,[false true],...
  'self_corr_flag',true,[false true],...
  'ref_corr_flag',false,[false true],...
  'ref_corr_pol',0,[-1 0 1],...
  'plotflag',true,[false true],...
  'plot_outdir','plots',[],...
...
  'color_order',{'b' 'g' 'r' 'c' 'm' 'y' 'k'},[],...
  'hemistrs',{'right','left'},{'right','left'},...
  'uplowstrs',{'upper','lower'},{'upper','lower'},...
  'areanames',{'V1','V2','V3'},{'V1','V2','V3'},...
  'wf_ylim',[-1.5,1.5],[],...
  'norm_wf_flag',true,[false true],...
  'cond_offsets_flag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.areas), parms.areas = 0; end;
if isempty(parms.hemis), parms.hemis = 0; end;
if isempty(parms.uplows), parms.uplows = 0; end;
if isempty(parms.eccs), parms.eccs = 0; end;

if ~isempty(parms.fname_out) && exist(parms.fname_out,'file') && ~parms.forceflag
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create reference waveform

% set full prefix
a=parms.areas(1);
h=parms.hemis(1);
u=parms.uplows(1);
e=parms.eccs(1);
prefix = parms.prefix_stem;
if h~=0
  prefix = sprintf('%s_hemi%d',prefix,h);
end;
if u~=0
  prefix = sprintf('%s_uplow%d',prefix,u);
end;
if e~=0
  prefix = sprintf('%s_ecc%d',prefix,e);
end;
if length(parms.contvec)==1
  prefix = sprintf('%s_cont%d',...
    prefix,parms.contvec);
elseif ~isempty(parms.contvec)
  prefix = sprintf('%s_cont%s',...
    prefix,sprintf('_%d',parms.contvec));
end;
if a~=0
  prefix = sprintf('%s_areas_%d',prefix,a);
end;
if parms.cond_offsets_flag
  prefix = sprintf('%s_condoffsets',prefix);
end;
matfile = sprintf('matfiles/%s_parms.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;

% load parameters from RCSE
tmp_parms = parms;
load(matfile);
ri_parms = parms;
parms = tmp_parms;
baseline_start = ri_parms.baseline_start_samp;
baseline_end = ri_parms.baseline_end_samp;
fit_start = ri_parms.fit_t0;
fit_end = ri_parms.fit_t1;

if parms.plotflag
  matfile = sprintf('matfiles/%s_avg_data.mat',prefix);
  load(matfile);
  time = avg_data.averages(1).time;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load offset waveforms

wforms = [];
err = [];
if max(parms.areas) > length(parms.areas)
  nareas = max(parms.areas);
else
  nareas = length(parms.areas);
end;
nhemis = length(parms.hemis);
nuplows = length(parms.uplows);
neccs = length(parms.eccs);
for a=parms.areas
  for h=parms.hemis
    for u=parms.uplows
      for e=parms.eccs
        if a==0, tmp_a=1; else tmp_a=a; end;
        if h==0, tmp_h=1; else tmp_h=h; end;
        if u==0, tmp_u=1; else tmp_u=u; end;
        if e==0, tmp_e=1; else tmp_e=e; end;
        prefix = parms.prefix_stem;
        if h~=0
          prefix = sprintf('%s_hemi%d',prefix,h);
        end;
        if u~=0
          prefix = sprintf('%s_uplow%d',prefix,u);
        end;
        if e~=0
          prefix = sprintf('%s_ecc%d',prefix,e);
        end;
        if length(parms.contvec)==1
          prefix = sprintf('%s_cont%d',...
            prefix,parms.contvec);
        elseif ~isempty(parms.contvec)
          prefix = sprintf('%s_cont%s',...
            prefix,sprintf('_%d',parms.contvec));
        end;
        if a~=0
          prefix = sprintf('%s_areas_%d',prefix,a);
        end;
        if parms.cond_offsets_flag
          prefix = sprintf('%s_condoffsets',prefix);
        end;
        matfile = sprintf('matfiles/%s_ret_forward.mat',prefix);
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,matfile);
          continue;
        end;
        clear retforward;
        load(matfile);
        tmp_wforms = retforward.retmap.offset_wforms;
        tmp_err = retforward.retmap.offset_errors;
        if isempty(wforms)
          matsize = size(tmp_wforms);
          gridsize = matsize(1:2);
          ngrid = prod(gridsize);
          if nareas==1 & parms.areas==0
            nwfareas = matsize(3);
            ntpoints = matsize(4);
            wforms = zeros([nwfareas,nhemis,nuplows,neccs,matsize([1,2,4])]);
            err = zeros([nwfareas,nhemis,nuplows,neccs,gridsize]);
          else
            ntpoints = matsize(3);
            wforms = zeros([nareas,nhemis,nuplows,neccs,matsize]);
            err = zeros([nareas,nhemis,nuplows,neccs,gridsize]);
          end;
        end;
        if nareas==1 & parms.areas==0
          for wa=1:nwfareas
            wforms(wa,tmp_h,tmp_u,tmp_e,:,:,:) = squeeze(tmp_wforms(:,:,wa,:));
            err(wa,tmp_h,tmp_u,tmp_e,:,:) = tmp_err;
          end;
        else
          wforms(tmp_a,tmp_h,tmp_u,tmp_e,:,:,:) = tmp_wforms;
          err(tmp_a,tmp_h,tmp_u,tmp_e,:,:) = tmp_err;
        end;
      end;
    end;
  end;
end;

if nareas==1 & parms.areas==0
  parms.areas = [1:nwfareas];
  nareas = nwfareas;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations with reference

if parms.ref_corr_flag
  fprintf('%s: calculating correlations to reference...\n',mfilename);

  % create reference waveform
  baseline_range = [baseline_start:baseline_end];
  fit_range = [fit_start:fit_end];
  nbaseline = length(baseline_range);
  nfit = length(fit_range);
  wform_ref = zeros(nbaseline+nfit,1);
  if parms.ref_corr_pol==0
    wform_ref(nbaseline+1:end) = 1;
  else
    wform_ref(nbaseline+1:end) = parms.ref_corr_pol;
  end;

  allR_ref = zeros([nareas,nhemis,nuplows,neccs,ngrid]);
  for a=parms.areas
    for h=parms.hemis
      for u=parms.uplows
        for e=parms.eccs
          if a==0, tmp_a=1; else tmp_a=a; end;
          if h==0, tmp_h=1; else tmp_h=h; end;
          if u==0, tmp_u=1; else tmp_u=u; end;
          if e==0, tmp_e=1; else tmp_e=e; end;
          tmp_wforms = squeeze(wforms(tmp_a,tmp_h,tmp_u,tmp_e,:,:,:));
          tmp_matsize = size(tmp_wforms);
          tmp_wforms = reshape(tmp_wforms,[ngrid,tmp_matsize(3)])';
          tmp_wforms = tmp_wforms([baseline_range,fit_range],:);
          R = corr(tmp_wforms,wform_ref);
          allR_ref(tmp_a,tmp_h,tmp_u,tmp_e,:) = R;
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create probability maps based on explained variance and ref corr
allmasks = zeros([nareas,nhemis,nuplows,neccs,gridsize]);
for a=parms.areas
  for h=parms.hemis
    for u=parms.uplows
      for e=parms.eccs
        if a==0, tmp_a=1; else tmp_a=a; end;
        if h==0, tmp_h=1; else tmp_h=h; end;
        if u==0, tmp_u=1; else tmp_u=u; end;
        if e==0, tmp_e=1; else tmp_e=e; end;
        tmp_expvar = 1-squeeze(err(tmp_a,tmp_h,tmp_u,tmp_e,:,:));
        if parms.expvar_thresh
          tmp_expvar(tmp_expvar<parms.expvar_thresh)=0;
          tmp_expvar(tmp_expvar>=parms.expvar_thresh)=1;
        else
          tmp_expvar(tmp_expvar<0)=0;
        end;
        tmp_expvar = tmp_expvar.^parms.expvar_fact;
        if parms.ref_corr_flag
          tmp_R_ref = reshape(allR_ref(tmp_a,tmp_h,tmp_u,tmp_e,:),gridsize);
          if parms.ref_corr_pol==0
            tmp_R_ref = abs(tmp_R_ref);
          end;
          if parms.ref_corr_thresh
            tmp_R_ref(tmp_R_ref<parms.ref_corr_thresh) = 0;
            tmp_R_ref(tmp_R_ref>=parms.ref_corr_thresh) = 1;
          else
            tmp_R_ref(tmp_R_ref<0) = 0;
          end;
          allmasks(tmp_a,tmp_h,tmp_u,tmp_e,:,:) = tmp_expvar.*tmp_R_ref;
        else
          allmasks(tmp_a,tmp_h,tmp_u,tmp_e,:,:) = tmp_expvar;
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlation to self

allR_self = zeros([nareas,nhemis,nuplows,neccs,ngrid,gridsize]);
if parms.self_corr_flag && parms.plotflag
  fprintf('%s: calculating correlations to self...\n',mfilename);
  tic
  for a=parms.areas
    for h=parms.hemis
      for u=parms.uplows
        for e=parms.eccs
          if a==0, tmp_a=1; else tmp_a=a; end;
          if h==0, tmp_h=1; else tmp_h=h; end;
          if u==0, tmp_u=1; else tmp_u=u; end;
          if e==0, tmp_e=1; else tmp_e=e; end;
          tmp_wforms = squeeze(wforms(tmp_a,tmp_h,tmp_u,tmp_e,:,:,:));
          tmp_matsize = size(tmp_wforms);
          tmp_wforms = reshape(tmp_wforms,[ngrid,tmp_matsize(3)])';
          R_self = corr(tmp_wforms,tmp_wforms);
          R_self = permute(reshape(R_self,[gridsize,ngrid]),[3,1,2]);
          allR_self(tmp_a,tmp_h,tmp_u,tmp_e,:,:,:) = R_self;
        end;
      end;
    end;
  end;
  toc
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations between hemispheres

if nhemis>1 & parms.hemi_corr_flag
  fprintf('%s: calculating correlations between hemispheres...\n',mfilename);
  tic
  allR_lvsr = zeros([nareas,nuplows,neccs,ngrid,ngrid]);
  allR = zeros([nareas,nhemis,nuplows,neccs,ngrid,gridsize]);
  for a=parms.areas
    for u=parms.uplows
      for e=parms.eccs
        if a==0, tmp_a=1; else tmp_a=a; end;
        if u==0, tmp_u=1; else tmp_u=u; end;
        if e==0, tmp_e=1; else tmp_e=e; end;
        wforms_right = squeeze(wforms(tmp_a,1,tmp_u,tmp_e,:,:,:));
        wforms_left = squeeze(wforms(tmp_a,2,tmp_u,tmp_e,:,:,:));
        tmp_matsize = size(wforms_right);
        wforms_right = reshape(wforms_right,[ngrid,tmp_matsize(3)])';
        wforms_left = reshape(wforms_left,[ngrid,tmp_matsize(3)])';
        R = corr(wforms_right,wforms_left);
        % right hemi 1D grid points vs. left hemi 2D offsets
        R_left = reshape(R,[ngrid,gridsize]);
        % left hemi 1D grid points vs. right hemi 2D offsets
        R_right = permute(reshape(R,[gridsize,ngrid]),[3,1,2]);
        allR_lvsr(tmp_a,tmp_u,tmp_e,:,:) = R;
        allR(tmp_a,1,tmp_u,tmp_e,:,:,:) = R_right;
        allR(tmp_a,2,tmp_u,tmp_e,:,:,:) = R_left;
      end;
    end;
  end;
  toc;
else
  % cannot calculate correlation between hemispheres with nhemis<2
  allR_lvsr = ones([nareas,nuplows,neccs,ngrid,ngrid]);
  allR = ones([nareas,nhemis,nuplows,neccs,ngrid,gridsize]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find offsets

alloffsets_1D = zeros([nareas,nhemis,nuplows,neccs]);
alloffsets_2D = zeros([nareas,nhemis,nuplows,neccs,2]);
for a=parms.areas
  for u=parms.uplows
    for e=parms.eccs
      if a==0, tmp_a=1; else tmp_a=a; end;
      if u==0, tmp_u=1; else tmp_u=u; end;
      if e==0, tmp_e=1; else tmp_e=e; end;
      if nhemis>1 & parms.hemi_corr_flag
        maskvec_right = reshape(allmasks(tmp_a,1,tmp_u,tmp_e,:,:),[ngrid,1]);
        maskvec_left = reshape(allmasks(tmp_a,2,tmp_u,tmp_e,:,:),[1,ngrid]);
        maskgrid = maskvec_right*maskvec_left;
        R = squeeze(allR_lvsr(tmp_a,tmp_u,tmp_e,:,:));
        R_masked = maskgrid.*R;
        [max_r,ind_max] = max(R_masked(:));
        [ind_right,ind_left] = ind2sub([ngrid,ngrid],ind_max);
        [ind_right_r,ind_right_th] = ind2sub(gridsize,ind_right);
        [ind_left_r,ind_left_th] = ind2sub(gridsize,ind_left);
        alloffsets_1D(tmp_a,1,tmp_u,tmp_e) = ind_right;
        alloffsets_1D(tmp_a,2,tmp_u,tmp_e) = ind_left;
        alloffsets_2D(tmp_a,1,tmp_u,tmp_e,1) = ind_right_r;
        alloffsets_2D(tmp_a,1,tmp_u,tmp_e,2) = ind_right_th;
        alloffsets_2D(tmp_a,2,tmp_u,tmp_e,1) = ind_left_r;
        alloffsets_2D(tmp_a,2,tmp_u,tmp_e,2) = ind_left_th;
      else
        for h=parms.hemis
          if h==0, tmp_h=1; else tmp_h=h; end;
          maskvec = reshape(allmasks(tmp_a,tmp_h,tmp_u,tmp_e,:,:),[ngrid,1]);
          if all(maskvec==1)
            ind = round(length(maskvec)/2);
          else
            [max_r,ind] = max(maskvec);
          end;
          [ind_r,ind_th] = ind2sub(gridsize,ind);
          alloffsets_1D(tmp_a,tmp_h,tmp_u,tmp_e) = ind;
          alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,1) = ind_r;
          alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,2) = ind_th;
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to csv file

% open csv file for output
if isempty(parms.fname_out)
  fid = 1; % write to stdout instead
else
  fid = fopen(parms.fname_out,'wt');
  if fid==-1
    error('failed to open %s for writing',parms.fname_out);
  end;
end;
fprintf(fid,'"area", "hemi", "uplow", "ecc", "r_offset", "th_offset", "min_error"\n');

for a=parms.areas
  for h=parms.hemis
    for u=parms.uplows
      for e=parms.eccs
        if a==0, tmp_a=1; else tmp_a=a; end;
        if h==0, tmp_h=1; else tmp_h=h; end;
        if u==0, tmp_u=1; else tmp_u=u; end;
        if e==0, tmp_e=1; else tmp_e=e; end;
        ind_r = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,1);
        ind_th = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,2);
        r_offset = retforward.retmap.r_offsets(ind_r);
        th_offset = retforward.retmap.th_offsets(ind_th);
        min_err = err(tmp_a,tmp_h,tmp_u,tmp_e,ind_r,ind_th);
        fprintf(fid,'%d, %d, %d, %d, %0.2f, %0.2f, %0.4f\n',...
          a,h,u,e,r_offset,th_offset,min_err);
      end;
    end;
  end;
end;
if fid~=1, fclose(fid); end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots

if parms.plotflag
  fprintf('%s: plotting results...\n',mfilename);
  f = 1;
  if ~strncmp(parms.plot_outdir,'/',1) % relative
    parms.plot_outdir = [pwd '/' parms.plot_outdir];
  end;
  if ~exist(parms.plot_outdir,'dir')
    [s,m]=mkdir(parms.plot_outdir);
    if ~s
      error('failed to create plot_outdir %s:\n%s',parms.plot_outdir,m);
    end;
  end;
  for e=parms.eccs
    if e==0, tmp_e=1; else tmp_e=e; end;
    for a=parms.areas
      if a==0, tmp_a=1; else tmp_a=a; end;
      figure(f); clf; p=1;
      plotname = 'explained variance';
      if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
      plotname = [parms.areanames{a} ' ' plotname];
      for u=parms.uplows
        for h=parms.hemis
          if u==0, tmp_u=1; else tmp_u=u; end;
          if h==0, tmp_h=1; else tmp_h=h; end;
          subplot(nuplows,nhemis,p); p=p+1;
          im = 1-squeeze(err(tmp_a,tmp_h,tmp_u,tmp_e,:,:));
          imagesc(im,parms.expvar_clim); colorbar;
          if h==0 & u~=0
            title(sprintf('%s visual field',parms.uplowstrs{u}));
          elseif h~=0 & u==0
            title(sprintf('%s visual field',parms.hemistrs{h}));
          elseif h~=0 & u~=0
            title(sprintf('%s %s visual field',parms.hemistrs{h},parms.uplowstrs{u}));
          end;
        end;
      end;
      suptitle(plotname);
      fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
      set(gcf,'Visible','off'); print('-dtiff',fname);

      if parms.ref_corr_flag
        figure(f); clf; p=1;
        plotname = 'correlation with reference';
        if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
        plotname = [parms.areanames{a} ' ' plotname];
        for u=parms.uplows
          for h=parms.hemis
            if u==0, tmp_u=1; else tmp_u=u; end;
            if h==0, tmp_h=1; else tmp_h=h; end;
            subplot(nuplows,nhemis,p); p=p+1;
            im = reshape(allR_ref(tmp_a,tmp_h,tmp_u,tmp_e,:),gridsize);
            imagesc(im,[-1,1]); colorbar;
            if h==0 & u~=0
              title(sprintf('%s visual field',parms.uplowstrs{u}));
            elseif h~=0 & u==0
              title(sprintf('%s visual field',parms.hemistrs{h}));
            elseif h~=0 & u~=0
              title(sprintf('%s %s visual field',parms.hemistrs{h},parms.uplowstrs{u}));
            end;
          end;
        end;
        suptitle(plotname);
        fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
        set(gcf,'Visible','off'); print('-dtiff',fname);
      end;

      figure(f); clf; p=1;
      plotname = 'dipole location probability';
      if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
      plotname = [parms.areanames{a} ' ' plotname];
      for u=parms.uplows
        for h=parms.hemis
          if u==0, tmp_u=1; else tmp_u=u; end;
          if h==0, tmp_h=1; else tmp_h=h; end;
          subplot(nuplows,nhemis,p); p=p+1;
          im = squeeze(allmasks(tmp_a,tmp_h,tmp_u,tmp_e,:,:));
          im = im/(eps+max(im(:)));
          imagesc(im); colorbar;
          if h==0 & u~=0
            title(sprintf('%s visual field',parms.uplowstrs{u}));
          elseif h~=0 & u==0
            title(sprintf('%s visual field',parms.hemistrs{h}));
          elseif h~=0 & u~=0
            title(sprintf('%s %s visual field',parms.hemistrs{h},parms.uplowstrs{u}));
          end;
        end;
      end;
      suptitle(plotname);
      fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
      set(gcf,'Visible','off'); print('-dtiff',fname);

      if nhemis>1 & parms.hemi_corr_flag
        % plot between-hemifield- and self-correlations with dipole offsets
        ylim = [-1,1];
        for h=parms.hemis
          figure(f); clf; p=1;
          plotname = sprintf('%s hemifield max cross-correlation',parms.hemistrs{h});
          if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
          plotname = [parms.areanames{a} ' ' plotname];
          for u=parms.uplows
            if u==0, tmp_u=1; else tmp_u=u; end;
            subplot(nuplows,3,p); p=p+1;
            im = zeros(gridsize);
            ind = alloffsets_1D(tmp_a,h,tmp_u,tmp_e);
            im(ind) = 1;
            imagesc(im,ylim); colorbar;
            title(sprintf('%s %s grid loc',parms.hemistrs{h},parms.uplowstrs{u}));

            for h2=parms.hemis
              subplot(nuplows,3,p); p=p+1;
              if h2==h
                im = squeeze(allR_self(tmp_a,h2,tmp_u,tmp_e,ind,:,:));
              else
                im = squeeze(allR(tmp_a,h2,tmp_u,tmp_e,ind,:,:));
              end;
              imagesc(im,ylim); colorbar;
              title(sprintf('%s %s corr',parms.hemistrs{h2},parms.uplowstrs{u}));
            end;
          end;  
          suptitle(plotname);
          fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
          set(gcf,'Visible','off'); print('-dtiff',fname);
        end;
      else
        % plot self-correlations with dipole offsets
        ylim = [-1,1];
        for h=parms.hemis
          if h==0, tmp_h=1; else tmp_h=h; end;
          figure(f); clf; p=1;
          if h~=0
            plotname = sprintf('%s hemifield dipole with self-correlation',...
              parms.hemistrs{h});
          else
            plotname = 'dipole with self-correlation';
          end;
          if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
          plotname = [parms.areanames{a} ' ' plotname];
          for u=parms.uplows
            if u==0, tmp_u=1; else tmp_u=u; end;
            subplot(nuplows,2,p); p=p+1;
            im = zeros(gridsize);
            ind = alloffsets_1D(tmp_a,tmp_h,tmp_u,tmp_e);
            im(ind) = 1;
            imagesc(im,ylim); colorbar;
            if h==0 & u~=0
              title(sprintf('%s grid loc',parms.uplowstrs{u}));
            elseif h~=0 & u==0
              title(sprintf('%s grid loc',parms.hemistrs{h}));
            elseif h~=0 & u~=0
              title(sprintf('%s %s grid loc',parms.hemistrs{h},parms.uplowstrs{u}));
            else
              title('grid loc');
            end;
            subplot(nuplows,2,p); p=p+1;
            im = squeeze(allR_self(tmp_a,tmp_h,tmp_u,tmp_e,ind,:,:));
            imagesc(im,ylim); colorbar;
            if h==0 & u~=0
              title(sprintf('%s corr',parms.uplowstrs{u}));
            elseif h~=0 & u==0
              title(sprintf('%s corr',parms.hemistrs{h}));
            elseif h~=0 & u~=0
              title(sprintf('%s %s corr',parms.hemistrs{h},parms.uplowstrs{u}));
            else
              title('corr');
            end;
          end;  
          suptitle(plotname);
          fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
          set(gcf,'Visible','off'); print('-dtiff',fname);
        end;
      end;
      
      figure(f); clf; p=1;
      plotname = 'source waveforms';
      if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
      plotname = [parms.areanames{a} ' ' plotname];
      ylim = parms.wf_ylim;
      for u=parms.uplows
        for h=parms.hemis
          if u==0, tmp_u=1; else tmp_u=u; end;
          if h==0, tmp_h=1; else tmp_h=h; end;
          subplot(nuplows,nhemis,p); p=p+1;
          ind_r = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,1);
          ind_th = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,2);
          wform = squeeze(wforms(tmp_a,tmp_h,tmp_u,tmp_e,ind_r,ind_th,:));
          if parms.norm_wf_flag
            wform = wform/abs(mean(wform(fit_start:fit_end)));
          end;
          plot(1000*time,wform);
          xlabel('time (ms)')
          axis tight;
          set(gca,'YLim',ylim);
          if h==0 & u~=0
            title(sprintf('%s visual field',parms.uplowstrs{u}));
          elseif h~=0 & u==0
            title(sprintf('%s visual field',parms.hemistrs{h}));
          elseif h~=0 & u~=0
            title(sprintf('%s %s visual field',parms.hemistrs{h},parms.uplowstrs{u}));
          end;
        end;
      end;
      suptitle(plotname);
      fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
      set(gcf,'Visible','off'); print('-dtiff',fname);
    end;

    % overlay area waveforms
    figure(f); clf; p=1;
    plotname = 'source waveforms';
    if e~=0, plotname = sprintf('ecc %d %s',plotname,e); end;
    ylim = parms.wf_ylim;
    for u=parms.uplows
      for h=parms.hemis
        if u==0, tmp_u=1; else tmp_u=u; end;
        if h==0, tmp_h=1; else tmp_h=h; end;
        subplot(nuplows,nhemis,p); p=p+1;
        hold on;
        for a=parms.areas
          if a==0, tmp_a=1; else tmp_a=a; end;
          ind_r = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,1);
          ind_th = alloffsets_2D(tmp_a,tmp_h,tmp_u,tmp_e,2);
          wform = squeeze(wforms(tmp_a,tmp_h,tmp_u,tmp_e,ind_r,ind_th,:));
          if parms.norm_wf_flag
            wform = wform/abs(mean(wform(fit_start:fit_end)));
          end;
          plot(1000*time,wform,parms.color_order{a});
        end;
        legend(parms.areanames(parms.areas),'Location','EastOutside');
        xlabel('time (ms)')
        axis tight;
        set(gca,'YLim',ylim);
        if h==0 & u~=0
          title(sprintf('%s visual field',parms.uplowstrs{u}));
        elseif h~=0 & u==0
          title(sprintf('%s visual field',parms.hemistrs{h}));
        elseif h~=0 & u~=0
          title(sprintf('%s %s visual field',parms.hemistrs{h},parms.uplowstrs{u}));
        end;
      end;
    end;
    suptitle(plotname);
    fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
    set(gcf,'Visible','off'); print('-dtiff',fname);
  end;
  close(gcf);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

