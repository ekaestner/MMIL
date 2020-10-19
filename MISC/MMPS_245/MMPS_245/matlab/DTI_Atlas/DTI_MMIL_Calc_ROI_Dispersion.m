function results = DTI_MMIL_Calc_ROI_Dispersion(ProjID,varargin)
%function results = DTI_MMIL_Calc_ROI_Dispersion(ProjID,[options])
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     if empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      proc_dti, fsurf
%     if both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%
% Other Optional Parameters:
%   'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'DTI_ROI_Dispersion'}
%   'outstem': output file stem
%     if full path, outdir will be ignored
%     {default = 'DTI_ROI'}
%   'minval': minimum value
%     {default = 1e-6}
%   'hist_flag': [0|1] create multi-subject histogram images for each ROI
%     {default = 0}
%   'hist_wtd_flag': [0|1] use ROI weights to calculate weighted histograms
%     ignored if hist_flag = 0
%     {default = 1}
%   'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Optional Parameters that specify which analysis results to compile:
%   'measlist': list of DT "measures" to extract for each fiber
%      e.g. 'FA', 'MD', 'LD', 'TD', 'b0', 'b0N', 'T1w'
%       'T1w' is the T1-weighted FreeSurfer nu.mgz
%        use 'none' to compile none of these measures
%      if more than one element in measlist, analysis should have been done
%        with minval = 'NaN' to ensure identical voxels used for each measure
%     { default = {'MD'} }
%   'scalefacts': scaling factors applied to each measure in measlist
%     if empty, will use 1 for each
%     if not empty, length must match length of measlist
%     {default = []}
%   'fiber_flag': [0|1] compile white matter fiber ROI results
%     {default = 1}
%   'aseg_flag': [0|1] compile aseg ROI results
%     {default = 0}
%   'wmparc_flag': [0|1] extract wmparc ROI results
%     {default = 0}
%   'cortsurf_flag': [0|1] compile cortical surface ROI results
%     {default = 0}
%
% Optional Parameters that determine input diffusion data:
%  'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 1}
%  'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet (must have run DT
%     calculations separately for each scan)
%     If empty, use DT measures calcd from all DTIScanNums
%     {default = []}
%  'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix: 'corr_resDTI'
%     {default = []}
%  'auto_infix_flag': [0|1] set infix automatically based on typical
%     processing and settings in ProjInfo
%     ignored if infix is not empty
%     {default = 1}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%  'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in DT measures analysis output
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     if 0, options above (snums_flag, snum_index, infix, etc.) are ignored
%     {default = 0}
%
% Optional Parameters specific to fiber analysis:
%   'fiber_indir': name of DTI fiber analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%   'fibers': vector of fiber numbers
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'atlas_flag': whether to use atlas fibers and if so,
%      what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'weighted_avg_flag': [0|1] use weighted average fiber ROI measures
%     {default = 1}
%   'thresh_FA': FA threshold used to calc fiber ROI measures
%     {default = 0}
%   'thresh_prob': fiber probability threshold
%     {default = 0}
%   'fiber_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'fiber_dispfact': multiplicative factor applied to fiber dispersion values
%     {default = 4}
%   'xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%     {default = 0}
%   'masksf_flag': [0|1] voxels with multiple fibers excluded
%     {default = 0}
%   'fiber_atlasname': name of the atlas used in AtlasTrack
%     and attached to output analysis mat files (empty = default atlas)
%     {default = []}
%   'fname_fiber_legend': name of fiber legend csv file
%     if empty, use default ($MMPS_PARMS/DTI_Fiber/DTI_Fiber_Legend.csv)
%     {default = []}
%
% Optional Parameters specific to aseg analysis:
%   'aseg_indir': name of DTI aseg analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%   'aseg_disp_flag': [0|1] whether to calculate weighted averages
%     based on MD values and dispvec
%     {default = 0}
%   'aseg_dispfact': multiplicative factor applied to aseg dispersion values
%     {default = 4}
%   'erode_flag': [0|1] whether aseg ROIs were "eroded" by
%     smoothing and thresholding to reduce partial voluming
%     {default = 1}
%   'erode_nvoxels': number of voxels eroded
%     {default = 1}
%   'aseg_aparc_flag': [0|1|2] whether cortical parcellation volume ROIs were used
%     0: aseg only
%     1: aparc only
%     2: aparc+aseg
%     {default = 0}
%
% Optional Parameters specific to cortical surface analysis:
%   'cortsurf_indir': name of DTI cortsurf analysis folder
%     relative to ContainerPath
%     {default = 'DTanalysis'}
%   'projdist': distance (mm) along normal vector
%     negative = white matter, positive = gray matter
%     {default = 1}
%
% Output:
%   results: struct containing fields:
%     dispersion: matrix of ROI dispersion values with size = [nrois,nmeas]
%       "dispersion" defined as median of robust estimates of standard deviation
%       calculated as median absolute deviation (mad) / 0.6745
%       averaged across hemispheres
%     centers: matrix of ROI center values with size = [nrois,nmeas]
%       "center" defined as median of medians, averaged across hemispheres
%     roicodes: vector of ROI codes
%     roinames: cell array of ROI names
%     roitypes: cell array of ROI types (i.e. fiber, aseg, cortsurf)
%     measlist: cell array of measure names
%     nrois: number of ROIs
%     nmeas: number of measures
%
% Created:  10/22/12 by Don Hagler
% Last Mod: 06/03/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
results = [];

parms = check_input(ProjID,varargin);

fname_out = [parms.outstem '_dispersion.mat'];
roivals = [];
if ~exist(fname_out,'file') || parms.forceflag
  % load ROI values for each subject
  roivals = compile_vals(parms);
  % retain useful values
  results = init_results(roivals);
  % calc ROI dispersion values
  [results.dispersion,results.centers] = calc_dispersion(roivals,parms);
  % average values across hemispheres
  [results.dispersion,results.centers] = average_hemispheres(results,parms);
  % save results
  save(fname_out,'results');
elseif nargout>0
  load(fname_out);
end;

if parms.hist_flag
  % load ROI values for each subject
  if isempty(roivals)
    roivals = compile_vals(parms);
  end;
  if isempty(results)
    load(fname_out);
  end;
  create_histograms(roivals,results,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
...
    'outdir','DTI_ROI_Dispersion',[],...
    'outstem','DTI_ROI',[],...
    'minval',1e-6,[],...
    'hist_flag',false,[false true],...
    'hist_wtd_flag',true,[false true],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
... % specify which results to compile
    'measlist',{'MD'},[],...
    'scalefacts',[],[],...
    'inputlist',[],[],...
    'fiber_flag',true,[false true],...
    'aseg_flag',false,[false true],...
    'wmparc_flag',false,[false true],...
    'cortsurf_flag',false,[false true],...
... % specify DTI data used for tensor calculations and analysis
    'snums_flag',1,[0:3],...
    'snum_index',[],[],...
    'infix',[],[],...
    'auto_infix_flag',true,[false true],...
    'revflag',0,[0:2],...
    'nob0_flag',false,[false true],...
    'min_bval',1,[],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'min_ndirs',6,[],...
    'full_fstem_flag',false,[false true],...
... % fiber analysis
    'fiber_indir','DTanalysis',[],...
    'fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'atlas_flag',2,[0:4],...
    'weighted_avg_flag',true,[false true],...
    'thresh_FA',0,[],...
    'thresh_prob',0,[],...
    'fiber_disp_flag',false,[false true],...
    'fiber_dispfact',4,[1e-6,1e6],...
    'xcg_flag',false,[false true],...
    'masksf_flag',false,[false true],...
    'fiber_atlasname',[],[],...
    'fname_fiber_legend',[],[],...
    'ICVnorm_flag',false,[false true],...
... % aseg analysis
    'aseg_indir','DTanalysis',[],...
    'aseg_disp_flag',false,[false true],...
    'aseg_dispfact',4,[1e-6,1e6],...
    'erode_flag',true,[false true],...
    'erode_nvoxels',1,[1:100],...
    'aseg_aparc_flag',0,[0,1,2],...
... % cortsurf analysis
    'cortsurf_indir','DTanalysis',[],...
    'projdist',1,[-5,5],...
... % hidden
    'aseg_disp_suffix','dwtd',[],...
    'fiber_disp_suffix','dwtd',[],...
    'xcg_suffix','xcg',[],...
    'masksf_suffix','masksf',[],...
    'hist_range',[],[],...
    'hist_nbins',100,[],...
    'min_nvals',3,[],...
    'fiber_pat','(?<hemi>[LR])[ _](?<name>[ \w]+)',[],...
    'aseg_pat','(?<hemi>[LeftRigh]+)-(?<name>[-\w]+)',[],...
    'cortsurf_pat','ctx-(?<hemi>[lr]h)-(?<name>[\w]+)',[],...
...
    'compile_tags',{'StudyInfo','RootDirs','measlist','scalefacts',...
                    'inputlist','fiber_flag','aseg_flag','wmparc_flag',...
                    'cortsurf_flag',...
                    'snums_flag','snum_index','infix','auto_infix_flag',...
                    'revflag','nob0_flag','min_bval','flex_flag',...
                    'min_nb0','min_ndirs',...
                    'full_fstem_flag','fiber_indir','fibers','atlas_flag',...
                    'weighted_avg_flag','thresh_FA','thresh_prob',...
                    'fiber_disp_flag','fiber_dispfact','xcg_flag',...
                    'masksf_flag','fiber_atlasname','fname_fiber_legend',...
                    'aseg_indir','aseg_disp_flag','aseg_dispfact',...
                    'erode_flag','erode_nvoxels','aseg_aparc_flag',...
                    'wmparc_aparc_flag',...
                    'cortsurf_indir','projdist_list','verbose',...
                    'aseg_disp_suffix','fiber_disp_suffix',...
                    'xcg_suffix','masksf_suffix'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);

  parms.projdist_list = parms.projdist;

  if ~iscell(parms.measlist)
    parms.measlist = {parms.measlist};
  end;

  if mmil_isrelative(parms.outstem)
    if mmil_isrelative(parms.outdir)
      parms.outdir = [getenv('HOME') '/MetaData/' parms.ProjID '/' parms.outdir];
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  else
    parms.outdir = fileparts(parms.outstem);
  end;

  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roivals = compile_vals(parms)
  roivals = [];
  fname_out = [parms.outstem '_vals.mat'];
  if ~exist(fname_out,'file') || parms.forceflag
    args = mmil_parms2args(parms,parms.compile_tags);
    roivals = MMIL_Compile_DTI_Analysis_Vals(parms.ProjID,args{:});
    save(fname_out,'roivals','-v7.3');
  else
    load(fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = init_results(roivals)
  results = [];
  tags = {'roicodes','roinames','roitypes','measlist','nrois','nmeas'};  
  for t=1:length(tags)
    tag = tags{t};
    results.(tag) = roivals.(tag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dispersion,centers] = calc_dispersion(roivals,parms);
  dispersion = nan(roivals.nrois,roivals.nmeas);
  centers = nan(roivals.nrois,roivals.nmeas);
  for r=1:roivals.nrois
    mad_mat = nan(roivals.nsubs,roivals.nmeas);
    med_mat = nan(roivals.nsubs,roivals.nmeas);
    for s=1:roivals.nsubs
      % calc median absolute deviation
      vals = roivals.StudyData(s).roi_data(r).vals;
      nvals = roivals.StudyData(s).roi_data(r).nvals;
      nvals_valid = roivals.StudyData(s).roi_data(r).nvals_valid;
      if nvals_valid<parms.min_nvals, continue; end;
      weights = mmil_getfield(roivals.StudyData(s).roi_data(r),'weights',[]);
      if ~isempty(weights) && ~all(weights==0)
        weights = roivals.StudyData(s).roi_data(r).weights;
        weights = repmat(weights,[1,roivals.nmeas]);
      else
        weights = ones(nvals,roivals.nmeas);
      end;
      weights(abs(vals)<parms.minval) = 0;
      [mad_mat(s,:),med_mat(s,:)] = mmil_wtd_mad(vals,weights);
    end;
    % exclude NaNs from median
    if any(isnan(mad_mat))
      for m=1:roivals.nmeas
        mad_vec = mad_mat(:,m);
        med_vec = med_mat(:,m);
        ind_valid = find(~isnan(mad_vec));
        if ~isempty(ind_valid)
          dispersion(r,m) = median(mad_vec(ind_valid));
          centers(r,m) = median(med_vec(ind_valid));
        end;
      end;
    else
      dispersion(r,:) = median(mad_mat,1);
      centers(r,:) = median(med_mat,1);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dispersion,centers] = average_hemispheres(results,parms)
  dispersion = results.dispersion;
  centers = results.centers;
  hemis = cell(results.nrois,1);
  roinames = results.roinames;
  for r=1:results.nrois
    roiname = results.roinames{r};
    roitype = results.roitypes{r};
    switch roitype
      case 'fiber'
        n = regexp(roiname,parms.fiber_pat,'names');
        if ~isempty(n)
          switch n.hemi
            case 'L'
              hemis{r} = 'lh';
            case 'R'
              hemis{r} = 'rh';
          end;
          roinames{r} = n.name;
        end;
      case 'aseg'
        n = regexp(roiname,parms.aseg_pat,'names');
        if ~isempty(n)
          switch n.hemi
            case 'Left'
              hemis{r} = 'lh';
            case 'Right'
              hemis{r} = 'rh';
            otherwise
              continue;
          end;
          roinames{r} = n.name;
        end;
      case 'cortsurf'
        n = regexp(roiname,parms.cortsurf_pat,'names');
        if ~isempty(n)
          hemis{r} = n.hemi;
          roinames{r} = n.name;
        end;
    end;
  end;
  matched_flags = zeros(results.nrois,1);
  for r=1:length(roinames)
    if matched_flags(r), continue; end;
    roiname = roinames{r};
    roitype = results.roitypes{r};
    ind = find(strcmp(roiname,roinames) & strcmp(roitype,results.roitypes));
    if length(ind)==2
      dispersion(ind,:) = repmat(mean(dispersion(ind,:),1),[2,1]);
      centers(ind,:) = repmat(mean(centers(ind,:),1),[2,1]);
      matched_flags(ind) = 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_histograms(roivals,results,parms);
  SubjIDs = {roivals.StudyInfo.SubjID};
  VisitIDs = {roivals.StudyInfo.VisitID};
  for m=1:roivals.nmeas
    if isempty(parms.hist_range)
      hist_range = [0,max(results.centers(:,m) + 3*results.dispersion(:,m))];
    else
      hist_range = parms.hist_range;
    end;
    hist_edges = ...
      linspace(hist_range(1),hist_range(2),parms.hist_nbins);
    nbins = length(hist_edges);
    for r=1:roivals.nrois
      fname_tif = sprintf('%s_%s_%s_hist.tif',...
        parms.outstem,roivals.measlist{m},roivals.roinames{r});
      fname_mat = sprintf('%s_%s_%s_hist.mat',...
        parms.outstem,roivals.measlist{m},roivals.roinames{r});
      if ~exist(fname_tif,'file') ||...
         ~exist(fname_mat,'file') || parms.forceflag
        hcmat = zeros(roivals.nsubs,nbins);
        for s=1:roivals.nsubs
          if roivals.StudyData(s).roi_data(r).nvals
            vals = roivals.StudyData(s).roi_data(r).vals(:,m);
            weights = roivals.StudyData(s).roi_data(r).weights;
            if all(weights==0), weights = ones(size(weights)); end;
            ind_valid = find(vals>=parms.minval);
            vals = vals(ind_valid);
            weights = weights(ind_valid);
            [h,bin] = histc(vals,hist_edges);
            if parms.hist_wtd_flag
              % replace bin numbers for out of range values
              weights(bin==0) = 0;
              bin(bin==0) = 1;
              % make sure output of accumarray includes all bins
              weights(end+1) = 0;
              bin(end+1) = 1;
              weights(end+1) = 0;
              bin(end+1) = parms.hist_nbins;
              bin = mmil_colvec(bin);
              weights = mmil_colvec(weights);
              % weighted histogram
              h = accumarray(bin,weights);
            end;
            hcmat(s,:) = h;
          end;
        end;
        save(fname_mat,'hcmat','SubjIDs','VisitIDs',...
                       'hist_range','hist_edges');
        figure;
        imagesc(hcmat);
        colormap hot;
        axis off;
        title = sprintf('%s %s',...
          roivals.measlist{m},regexprep(roivals.roinames{r},'_',' '));
        set(gcf,'visible','off');
        print(gcf,'-dtiff',fname_tif);
        close(gcf);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
