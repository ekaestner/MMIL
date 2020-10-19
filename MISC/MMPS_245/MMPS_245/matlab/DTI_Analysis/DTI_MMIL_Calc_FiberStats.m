function DTI_MMIL_Calc_FiberStats(ProjID,varargin);
%function DTI_MMIL_Calc_FiberStats(ProjID,varargin);
%
% Usage:
%  DTI_MMIL_Calc_FiberStats(ProjID,'key1', value1,...);
%
% Required Parameters:
%   ProjID: project ID to run (e.g. 'REC_TEST')
%
% Optional Parameters:
%   'StudyInfo': struct array of study information
%               If empty, will use ProjID to get StudyInfo  {default = []}
%   'RootDirs': root directory struct -- should have proc_dti field
%   'FiberLegend': struct array containing information about fiber numbers
%     (see DTI_MMIL_Read_Fiber_Legend)
%   'outdir': output directory
%     {default = []}
%   'plot_flag': [0|1|2|3] whether to plot results
%     and generate image files
%     0: no plots
%     1: some plots
%     2: lots of plots
%     3: lots of plots and scatter plots
%     {default = 1}
%   'plot_eps_flag': [0|1] whether to save plots as eps files
%     otherwise tiffs (ignored if plot_flag = 0)
%     {default = 0}
%   'save_csv_flag': [0|1] whether to save results as csv file
%     {default = 0}
%   'fibernums': list of fiber numbers to include in plots
%     {default = [101:110,115:123,133:138,141:150,2000:2004]}
%   'forceflag': overwrite existing output files
%     {default = 0}
%
% Optional Parameters that specify how DT ROI measures were generated:
%   'measlist': cell array of DT-derived measures
%     {default = {'FA','MD','LD','TD','b0N'}}
%   'scalefacts': vector of scaling factors for each measure in measlist
%     if a single value, will be applied to all measures
%     {default = 1}
%   'atlas_flag': whether atlas fibers were used and if so,
%      what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%     {default = 2}
%   'snums_flag': [0|1|2|3] which set of "snums" to use
%      0: use all available scans
%      1: use scan numbers in StudyInfo.DTIScanNums
%      2: use scan numbers in StudyInfo.DTIScanNums2
%      3: use scan numbers in StudyInfo.DTIScanNums3
%     {default = 0}
%   'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2 or DTIScanNums3) were used in spreadsheet
%     (if empty, then all)
%     {default = []}
%   'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%   'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1000}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'weighted_avg_flag': [0|1] whether weighted averages were calculated for
%      each ROI using fiber counts or fiber probabilities to weight
%      contribution from each voxel
%     {default = 1}
%   'thresh_FA': mask fiber ROIs by applying threshold to FA image
%     {default = 0}
%   'thresh_prob': mask atlas-derived fiber ROIs by applying
%     probability threshold
%     {default = 0}
%   'xcg_flag': [0|1] CSF and gray-matter voxels excluded
%     using FreeSurfer-derived aseg segmentation
%     {default = 1}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     {default = 0}
%   'ICVnorm_flag': [0|1] whether fiber volumes were normalized by
%       intracranial volume as calculated by freesurfer
%       (applies only to 'meas' = 'vol')
%     {default = 1}
%   'fiber_atlasname': name of atlas used in AtlasTrack 
%     if empty, uses the default atlas
%     {default = []}
%
% Created:  04/04/11 by Vijay Venkatraman
% Last Mod: 02/04/13 by Don Hagler
%

%% todo: create subfunctions to improve readability
%% todo: create functions that can be shared with separate RSI/FOD function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
  'outdir',[],[],...
  'FiberLegend',[],[],...
  'plot_flag',1,[0 1 2 3],...
  'plot_eps_flag',false,[false true],...
  'save_csv_flag',true,[false true],...
  'forceflag',false,[false true],...
  'measlist',{'DT_FA','DT_MD','DT_LD','DT_TD','DT_b0N'},[],...
  'scalefacts',1,[-Inf,Inf],...
  'atlas_flag',2,[0:4],...
  'snums_flag',0,[0:3],...
  'snum_index',[],[],...
  'nob0_flag',false,[false true],...
  'min_bval',1000,[],...
  'flex_flag',false,[false true],...
  'weighted_avg_flag',true,[false true],...
  'thresh_FA',0,[0,1],...
  'thresh_prob',0,[0,1],...
  'xcg_flag',true,[false true],...
  'masksf_flag',false,[false true],...
  'ICVnorm_flag',false,[false true],...
  'fibernums',[101:110,115:123,133:138,141:150,2000:2004],[],...
  'required_rootdirs',{'proc_dti'},[],...
  'qcflag',true,[false true],...
  'modality',{'MRI'},{'MRI','MEG','PET'},...
... % hidden parameters:
  'age_flag',true,[false true],... % whether to do ANCOVA with age as regressor
  'xlim_FA',[0.2,0.6],[],...
  'xlim_MD',[0.5e-3,1.5e-3],[],...
  'infix','corr_regT1',[],...
  'revflag',0,[0,1,2],...
  'outlier_flag',false,[false true],...
  'outlier_nSD',2,[1,100],... % number of standard deviations away from mean for rejection criteria
  'min_ndirs',6,[],...
  'legend_flag',true,[false true],...
  'annot_flag',true,[false true],...
  'paper_width',5.5,[1,10],...
  'paper_height',4,[1,10],...
  'dpi',250,[100,600],...
  'linewidth',1,[0.1,4],...
  'fontsize',6,[1,100],...
  'barwidth',1,[0.01,100],...
  'outfix',[],[],...
  'title_flag',true,[false true],...
  'xcg_suffix','xcg',[],...
  'masksf_suffix','masksf',[],...
  'fiber_atlasname',[],[],...
...
  'compile_tags',{'StudyInfo','infix','snums_flag','snum_index','revflag',...
                  'nob0_flag','min_bval','flex_flag','atlas_flag','measlist',...
                  'weighted_avg_flag','thresh_FA','thresh_prob','xcg_flag',...
                  'masksf_flag','xcg_suffix','masksf_suffix','outdir',...
                  'min_ndirs','fiber_atlasname','forceflag'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args= MMIL_Args(parms,'MMIL_Check_ProjID');
if isempty(parms.StudyInfo) || isempty(parms.RootDirs)
  [ProjInfo,parms.StudyInfo,parms.RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
end;
DEFAULT_MIN_BVAL = 1;
DEFAULT_MIN_NDIRS = 6;
DEFAULT_GROUPA='ace'; %lower case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.outdir)
  parms.outdir=[parms.RootDirs.home '/FiberStats'];
  if ~exist(parms.outdir,'dir'), mkdir(parms.outdir); end;
end;

if isempty(parms.FiberLegend)
  parms.FiberLegend= DTI_MMIL_Read_Fiber_Legend;
end;

if isempty(parms.measlist)
  error('measlist is empty')
end;

if ~iscell(parms.measlist), parms.measlist = {parms.measlist}; end;
if length(parms.scalefacts)==1
  parms.scalefacts = parms.scalefacts*ones(length(parms.measlist),1);
elseif length(parms.scalefacts)~=length(parms.measlist)
  error('number of scalefacts (%d) does not match number of measures (%d)',...
    length(parms.scalefacts),length(parms.measlist));
end;

% exclude specified fibers, reorder
legend_fibernums = cell2mat({parms.FiberLegend.FiberNumber});
legend_fibernames = {parms.FiberLegend.FiberName};
legend_rnums = cell2mat({parms.FiberLegend.ReorderNumber});
[legend_rnums,ind_reorder] = sort(legend_rnums);
legend_fibernums = legend_fibernums(ind_reorder);
legend_fibernames = legend_fibernames(ind_reorder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set file name infixes (must match those set by DTI_MMIL_Analyze_Fibers_Exam)
switch parms.atlas_flag
  case 0 % manual
    outstr = 'fibers';
    figtitle = 'fibers';
  case 1 % loc only, count atlas
    outstr = 'fibers_from_loc_countatlas';
    figtitle = 'fibers from loc countatlas';
  case 2 % loc+dir, count atlas
    outstr = 'fibers_from_countatlas';
    figtitle = 'fibers from countatlas';
  case 3 % loc only, mask atlas
    outstr = 'fibers_from_loc_atlas';
    figtitle = 'fibers from loc mask atlas';
  case 4 % loc+dir, mask atlas
    outstr = 'fibers_from_atlas';
    figtitle = 'fibers from mask atlas';
end;    

if parms.weighted_avg_flag
  outstr = [outstr '_wtd'];
  figtitle = [figtitle ' weighted'];
end;

if parms.thresh_prob>0 && parms.atlas_flag>0
  outstr = sprintf('%s_pthresh%0.2f',outstr,parms.thresh_prob);
  figtitle = sprintf('%s pthresh%0.2f',figtitle,parms.thresh_prob);
end;
if parms.thresh_FA>0
  outstr = sprintf('%s_FAthresh%0.2f',outstr,parms.thresh_FA);
  figtitle = sprintf('%s FAthresh%0.2f',figtitle,parms.thresh_FA);
end;
if parms.xcg_flag
  outstr = [outstr '_' parms.xcg_suffix];
  figtitle = [figtitle ' ' parms.xcg_suffix];
end;
if parms.masksf_flag
  outstr = [outstr '_' parms.masksf_suffix];
  figtitle = [figtitle ' ' parms.masksf_suffix];
end;
if parms.snums_flag
  outstr = sprintf('%s_snums%d',outstr,parms.snums_flag);
  figtitle = sprintf('%s snums%d',figtitle,parms.snums_flag);
end;
if ~isempty(parms.snum_index)
  outstr = sprintf('%s_snumidx%d',outstr,parms.snum_index);
  figtitle = sprintf('%s snumidx%d',figtitle,parms.snum_index);
end;
if ~isempty(parms.fiber_atlasname)
  figtitle= sprintf('%s Atlas: %s',figtitle,parms.fiber_atlasname);
  outstr= sprintf('%s_%s',outstr,parms.fiber_atlasname);
end;
if parms.nob0_flag
  outstr = [outstr '_nob0'];
  figtitle = [figtitle ' nob0'];
end;
if parms.min_bval ~= DEFAULT_MIN_BVAL
  outstr = sprintf('%s_minb%d',outstr,parms.min_bval);
  figtitle = sprintf('%s minb%d',figtitle,parms.min_bval);
end
if parms.flex_flag
  outstr = sprintf('%s_flex',outstr);
  figtitle = sprintf('%s flex',figtitle);
end;
if parms.min_ndirs ~= DEFAULT_MIN_NDIRS
  outstr = sprintf('%s_mind%d',outstr,parms.min_ndirs);
end
if parms.outlier_flag
  outstr = sprintf('%s_outlier%0.1fSD',outstr,parms.outlier_nSD);
end;

if ~parms.title_flag, figtitle=[]; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s: running DTIAtlas_MMIL_Calc_FiberStats...\n',mfilename);
fprintf('      atlas=%d, weighted=%d, thresh_prob=%0.2f, thresh_FA=%0.2f, AtlasName= %s\n',...
    parms.atlas_flag,parms.weighted_avg_flag,parms.thresh_prob,parms.thresh_FA,parms.fiber_atlasname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compile data for all subjects
args = mmil_parms2args(parms,parms.compile_tags);
matfiles = DTI_MMIL_Compile_ROIs(parms.RootDirs,args{:});
if isempty(matfiles), error('data compilation failed'); end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:length(parms.measlist)
  meas = parms.measlist{m};
  scalefact = parms.scalefacts(m);
  if parms.ICVnorm_flag & strcmp(meas,'vol')
    outstr = ['ICVnorm_' outstr];
  end;
  fname_out = sprintf('%s/%s_%s_groupstats.mat',parms.outdir,...
    meas,outstr);
  if ~exist(fname_out,'file') || parms.forceflag
    % load compiled ROI data
    fname_in = matfiles{m};
    clear group_data;
    load(fname_in);
    group_data = select_subjects(group_data,parms.age_flag);
    nsubs = length(group_data);
    if nsubs==0
      error('no valid subjects found');
    end;

    fprintf('%s: calculating %s group stats for %d subjects...\n',...
      mfilename,meas,nsubs);
    fiberstats = [];
    subject_inclusion = []; within_group = []; between_group = [];
    nrois = length(group_data(1).fiber_data);
    fibernums = cell2mat({group_data(1).fiber_data.fibernum});
    groups = unique({group_data.Group});
    for f=1:nrois
      TTEST = []; ANOVA = []; ANCOVA = [];
      fibernum = fibernums(f);
      ind = find(fibernum==legend_fibernums);
      if isempty(ind)
        error('fiber number %d not found in legend',fibernum);
      end;
      fibername = legend_fibernames{ind};
      nsubs = length(group_data);
      ind_good_subs = [1:nsubs];
      % exclude NaNs
      ind_bad_subs = [];
      for s=1:nsubs
        SubjID = group_data(s).SubjID;
        if isempty(group_data(s).fiber_data)
          if f==1
            fprintf('%s: WARNING: empty fiber_data for subject %s\n',...
              mfilename,SubjID);
          end;
          ind_bad_subs = [ind_bad_subs s];
          continue;
        end;
        val = group_data(s).fiber_data(f).avg*scalefact;
        % exclude NaNs
        if any(isnan(val))
          fprintf('%s: WARNING: NaN for subject %s, fiber %d (%s)\n',...
            mfilename,SubjID,fibernum,fibername);
          ind_bad_subs = [ind_bad_subs s];
          continue;
        end;
      end;
      ind_good_subs = setdiff(ind_good_subs,ind_bad_subs);
      subject_inclusion(f).fibernum = fibernum;
      subject_inclusion(f).fibername = fibername;
      subject_inclusion(f).NaNs = {};
      for s=ind_bad_subs
        subject_inclusion(f).NaNs{end+1} = group_data(s).SubjID;
      end;
      % exclude outliers
      if parms.outlier_flag
        % calculate mean and stdv across all subjects
        tmp_mean = 0;
        tmp_stdv = 0;
        tmp_N = 0;
        for s=ind_good_subs
          SubjID = group_data(s).SubjID;
          val = group_data(s).fiber_data(f).avg*scalefact;
          tmp_mean = tmp_mean + val;
          tmp_stdv = tmp_stdv + val.*val;
          tmp_N = tmp_N + 1;
        end;
        if tmp_N <= 1
          error('must have more than 1 subject to do outlier rejection');
        end;
        tmp_stdv = sqrt((tmp_N*tmp_stdv - tmp_mean.^2)./(eps+tmp_N*(tmp_N-1)));
        tmp_mean = tmp_mean/tmp_N;
        % find outliers
        ind_bad_subs = [];
        for s=ind_good_subs
          val = group_data(s).fiber_data(f).avg*scalefact;
          zscore = (val-tmp_mean)./tmp_stdv;
          if any(abs(zscore)>parms.outlier_nSD)
            ind_bad_subs = [ind_bad_subs s];
          end;
        end;
        ind_good_subs = setdiff(ind_good_subs,ind_bad_subs);
        subject_inclusion(f).outliers = {};
        for s=ind_bad_subs
          subject_inclusion(f).outliers{end+1} = group_data(s).SubjID;
        end;
      end;
      subject_inclusion(f).included = {};
      for s=ind_good_subs
        subject_inclusion(f).included{end+1} = group_data(s).SubjID;
      end;
      nsubs = length(ind_good_subs);

      for g=1:length(groups)
        within_group = setfield(within_group,{f,g},'fibernum',fibernum);
        within_group = setfield(within_group,{f,g},'fibername',fibername);
        within_group = setfield(within_group,{f,g},'group',groups{g});
        within_group = setfield(within_group,{f,g},'mean',0);
        within_group = setfield(within_group,{f,g},'stdv',0);
        within_group = setfield(within_group,{f,g},'N',0);
      end;
      for s=ind_good_subs
        g = find(strcmp(group_data(s).Group,groups));
        tmp_mean = getfield(within_group,{f,g},'mean');
        tmp_stdv = getfield(within_group,{f,g},'stdv');
        tmp_N = getfield(within_group,{f,g},'N');
        val = group_data(s).fiber_data(f).avg*scalefact;
        if isnan(val), continue; end;
        tmp_mean = tmp_mean + val;
        tmp_stdv = tmp_stdv + val.*val;
        tmp_N = tmp_N + 1;
        
        within_group = setfield(within_group,{f,g},'mean',tmp_mean);
        within_group = setfield(within_group,{f,g},'stdv',tmp_stdv);
        within_group = setfield(within_group,{f,g},'N',tmp_N);
      end;

      % calculate within group stats (mean, stdv) from sum and sum of squares
      for g=1:length(groups)
        tmp_mean = getfield(within_group,{f,g},'mean');
        tmp_stdv = getfield(within_group,{f,g},'stdv');
        tmp_N = getfield(within_group,{f,g},'N');
        if tmp_N > 1
          tmp_stdv = sqrt((tmp_N*tmp_stdv - tmp_mean.^2)./(eps+tmp_N*(tmp_N-1)));
        else
          tmp_stdv = [];
        end;
        if tmp_N > 0
          tmp_mean = tmp_mean/tmp_N;
        else
          tmp_mean = [];
        end;
        tmp_stderr = tmp_stdv/sqrt(tmp_N);
        tmp_cv = tmp_stdv/(tmp_mean+eps);
        within_group = setfield(within_group,{f,g},'mean',tmp_mean);
        within_group = setfield(within_group,{f,g},'stdv',tmp_stdv);
        within_group = setfield(within_group,{f,g},'sterr',tmp_stderr);
        within_group = setfield(within_group,{f,g},'cv',tmp_cv);
        within_group = setfield(within_group,{f,g},'N',tmp_N);
      end;

      % calculate between-group stats (t-test, ANOVA, ANCOVA)
      if length(groups)>1
        % pairwise t-tests
        subj_groups = {group_data.Group};
        k = 1;
        for g=1:length(groups)
          for h=g+1:length(groups)
            groupA = groups{g};
            groupB = groups{h};
            indA = find(strcmp(groupA,subj_groups));
            indB = find(strcmp(groupB,subj_groups));
            valsA = [];
            valsB = [];
            for i=1:length(indA)
              s = indA(i);
              if ~ismember(s,ind_good_subs), continue; end;
              val = group_data(s).fiber_data(f).avg*scalefact;
              valsA = [valsA val];
            end;
            for i=1:length(indB)
              s = indB(i);
              if ~ismember(s,ind_good_subs), continue; end;
              val = group_data(s).fiber_data(f).avg*scalefact;
              valsB = [valsB val];
            end;
            [h,p,ci,stats] = ttest2(valsA,valsB);
            plog = -log10(p);
            d = cohensd(valsA,valsB);
%           [amd_d,ci] = andersd(valsA,valsB);
            TTEST(k).groupA = groupA;
            TTEST(k).groupB = groupB;
            TTEST(k).p = p;
            TTEST(k).plog = plog;
            TTEST(k).d = d;
            TTEST(k).t = stats.tstat;
            k = k + 1;
          end;
        end;

        % ANOVA
        vals = zeros(1,nsubs);
        grouplist = {};
        for i=1:nsubs
          s = ind_good_subs(i);
          vals(i) = group_data(s).fiber_data(f).avg*scalefact;
          grouplist{i} = group_data(s).Group;
        end;
        [p,anovatab,stats] = anova1(vals,grouplist,'off');
        ANOVA.p = p;
        ANOVA.plog = -log10(p);
        ANOVA.F = anovatab{2,5};
        SSeffect = anovatab{2,2};
        SStotal = anovatab{4,2};
        ANOVA.eta2 = SSeffect/SStotal; % effect size eta squared
        ANOVA.anovatab = anovatab;
        ANOVA.stats = stats;

        % ANCOVA with age as regressor
        if parms.age_flag
          ages = zeros(1,nsubs);
          for i=1:nsubs
            s = ind_good_subs(i);
            ages(i) = group_data(s).Age;
          end;
          [h,anovatab,ctab,stats] = aoctool(ages,vals,grouplist,...
            0.05,[],[],[],'off','parallel lines');
          p = anovatab{2,6};
          p_age = anovatab{3,6};
          ANCOVA.p = p;
          ANCOVA.plog = -log10(p);
          ANCOVA.p_age = p_age;
          ANCOVA.plog_age = -log10(p_age);
          ANCOVA.F = anovatab{2,5};
          SSeffect = anovatab{2,3};
          SStotal = anovatab{2,3} + anovatab{3,3} + anovatab{4,3};
          ANCOVA.eta2 = SSeffect/SStotal; % effect size eta squared
          ANCOVA.ctab = ctab;
          ANCOVA.anovatab = anovatab;
          ANCOVA.stats = stats;
        end;
      end;
      between_group(f).fibernum = fibernum;
      between_group(f).fibername = fibername;
      between_group(f).TTEST = TTEST;
      between_group(f).ANOVA = ANOVA;
      between_group(f).ANCOVA = ANCOVA;
    end;

    % include fibers in common between legend and stats, reorder
    fibernums = cell2mat({within_group(:,1).fibernum});
    [tmp_fibernums,ind_leg,ind_stats] = intersect(legend_fibernums,fibernums);    
    [tmp,ind] = sort(ind_leg);
    ind_reorder = ind_stats(ind);
    fiberstats.within_group = within_group(ind_reorder,:);
    fiberstats.between_group = between_group(ind_reorder);
    fiberstats.subject_inclusion = subject_inclusion(ind_reorder);
    nfibers = length(ind_reorder);
    if ~nfibers
      error('no fibers with valid fiber numbers');
    end;

    fprintf('%s: saving results to %s...\n',mfilename,fname_out);
    save(fname_out,'fiberstats');
  end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save fiberstats as csv file (comma separated value spreadsheet)
if parms.save_csv_flag
  fprintf('%s: saving fiberstats to csv...\n',mfilename);
  for m=1:length(parms.measlist)
    meas = parms.measlist{m};
    scalefact = parms.scalefacts(m);
    if parms.ICVnorm_flag & strcmp(meas,'vol')
      outstr = ['ICVnorm_' outstr];
      if ~isempty(figtitle)
        tmp_figtitle = [figtitle ' ICVnorm'];
      end;
    else
      tmp_figtitle = figtitle;
    end;
    fname_in = sprintf('%s/%s_%s_groupstats.mat',parms.outdir,...
      meas,outstr);
    fname_out = sprintf('%s/%s_%s_groupstats.csv',parms.outdir,...
      meas,outstr);
    if ~exist(fname_out,'file') || parms.forceflag
      clear fiberstats;
      load(fname_in);
      groups = unique({fiberstats.within_group.group});
      fid = fopen(fname_out,'wt');
      if fid==-1
        error('failed to open file %s for writing',fname_out);
      end;
      % write column headers
      fprintf(fid,'"FiberName","FiberNumber"');
      for g=1:length(groups)
        group_name = upper(groups{g});
        fprintf(fid,',"%s N","%s Mean","%s StdDev","%s StdErr","%s CV"',...
          group_name,group_name,group_name,group_name,group_name);
      end;
      for t=1:length(fiberstats.between_group(1).TTEST)
        groupA = upper(fiberstats.between_group(1).TTEST(t).groupA);
        groupB = upper(fiberstats.between_group(1).TTEST(t).groupB);
        fprintf(fid,',"%s VS %s tstat","%s VS %s pval","%s VS %s Cohen''s d"',...
          groupA,groupB,groupA,groupB,groupA,groupB);
      end;
      fprintf(fid,',"ANOVA Fstat","ANOVA pval","ANOVA eta squared"');
      if parms.age_flag
        fprintf(fid,',"ANCOVA Fstat","ANCOVA pval","ANCOVA eta squared","ANCOVA age pval"');
      end;
      fprintf(fid,'\n');

      for f=1:length(fiberstats.within_group)
        fprintf(fid,'"%s",%d',fiberstats.within_group(f).fibername,fiberstats.within_group(f).fibernum);
        for g=1:length(groups)
          fprintf(fid,',%d,%0.4f,%0.4f,%0.4f,%0.4f',...
            fiberstats.within_group(f,g).N,...
            fiberstats.within_group(f,g).mean,...
            fiberstats.within_group(f,g).stdv,...
            fiberstats.within_group(f,g).sterr,...
            fiberstats.within_group(f,g).cv);
        end;
        for t=1:length(fiberstats.between_group(f).TTEST)
          fprintf(fid,',%0.4f,%0.4f,%0.4f',...
            fiberstats.between_group(f).TTEST(t).t,...
            fiberstats.between_group(f).TTEST(t).plog,...
            fiberstats.between_group(f).TTEST(t).d);
        end;
        fprintf(fid,',%0.4f,%0.4f,%0.4f',...
          fiberstats.between_group(f).ANOVA.F,...
          fiberstats.between_group(f).ANOVA.plog,...
          fiberstats.between_group(f).ANOVA.eta2);
        if parms.age_flag
          fprintf(fid,',%0.4f,%0.4f,%0.4f,%0.4f',...
            fiberstats.between_group(f).ANCOVA.F,...
            fiberstats.between_group(f).ANCOVA.plog,...
            fiberstats.between_group(f).ANCOVA.eta2,...
            fiberstats.between_group(f).ANCOVA.plog_age);
        end;
        fprintf(fid,'\n');
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot fiberstats
if parms.plot_flag
  fprintf('%s: plotting fiberstats...\n',mfilename);
  figure('Visible','off');
  for m=1:length(parms.measlist)
    meas = parms.measlist{m};
    scalefact = parms.scalefacts(m);
    if parms.ICVnorm_flag & strcmp(meas,'vol')
      outstr = ['ICVnorm_' outstr];
      if ~isempty(figtitle)
        tmp_figtitle = [figtitle ' ICVnorm'];
      end;
    else
      tmp_figtitle = figtitle;
    end;
    fname_in = sprintf('%s/%s_%s_groupstats.mat',parms.outdir,...
      meas,outstr);
    clear fiberstats;
    load(fname_in);
    groups = unique({fiberstats.within_group.group});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % means
    fibernums = cell2mat({fiberstats.within_group(:,1).fibernum});
    ind_fibers = find(ismember(fibernums,parms.fibernums));
    data_names = {fiberstats.within_group(ind_fibers,1).fibername};
    data = zeros(size(fiberstats.within_group(ind_fibers,:)));
    err = zeros(size(fiberstats.within_group(ind_fibers,:)));
    for g=1:length(groups)
      data(:,g) = cell2mat({fiberstats.within_group(ind_fibers,g).mean});
      err(:,g) = cell2mat({fiberstats.within_group(ind_fibers,g).sterr});
    end;
    if parms.legend_flag
      if parms.annot_flag
        group_names = groups;
        for g=1:length(groups)
          str = upper(groups{g});
          N = cell2mat({fiberstats.within_group(ind_fibers,g).N});
          cv = cell2mat({fiberstats.within_group(ind_fibers,g).cv});
          cv_mean = mean(cv);
          cv_stderr = stderr(cv);
          str = sprintf('%s\nN=%d',str,min(N));
          if min(N)~=max(N), str = sprintf('%s to %d',str,max(N)); end;
          str = sprintf('%s\ncv = %0.2f +- %0.2f',str,cv_mean,cv_stderr);
          group_names{g} = str;
        end;
      end;
    else
      group_names = [];
    end;
    figxlabel = sprintf('mean %s',regexprep(meas,'_',' '));
    switch meas
      case 'FA'
        xlim = parms.xlim_FA;
      case 'MD'
        xlim = parms.xlim_MD;
      otherwise
        xlim = [];
    end;
    color_order = {'k','b','g'};
    clf;
    mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
      'barwidth',parms.barwidth,...
      'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
      'title',tmp_figtitle,'xlabel',figxlabel,'xlim',xlim,...
      'color_order',color_order);
    set(gcf,'Visible','off');
    fstem = sprintf('%s/%s_%s_mean',...
      parms.outdir,meas,outstr);
    if ~isempty(parms.outfix)
      fstem = [fstem '_' parms.outfix];
    end;
    fname_jpg = [fstem '.jpg'];
    set(gcf,'PaperUnits','inches',...
            'PaperPosition',[0 0 parms.paper_width parms.paper_height])
    print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
    if parms.plot_eps_flag
      fname_eps = [fstem '.eps'];
      print(gcf,'-deps',fname_eps);
    end;

    if isempty(fiberstats.between_group), return; end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % effect sizes for ONE GROUP (usually 'control') vs. other groups
    data = [];
    err = []; % what are error bars for effect size?
    group_names = [];
    for i=1:length(fiberstats.between_group(1).TTEST)
      groupA = fiberstats.between_group(1).TTEST(i).groupA;
      groupB = fiberstats.between_group(1).TTEST(i).groupB;
      if ~strcmp(lower(groupA),DEFAULT_GROUPA), continue; end;
      group_names{i} = groupB;
      tmp_vec = [];
      for j=1:length(ind_fibers)
        f = ind_fibers(j);
        tmp_vec = [tmp_vec,-fiberstats.between_group(f).TTEST(i).d];
      end;
      data = [data;tmp_vec];
    end;
    if ~isempty(data)
      data = data';
      if parms.legend_flag
        if parms.annot_flag
          for i=1:length(group_names)
            str = upper(group_names{i});
            g = find(strcmp(group_names{i},groups));
            N = cell2mat({fiberstats.within_group(:,g).N});
            str = sprintf('%s\nN=%d',str,min(N));
            if min(N)~=max(N), str = sprintf('%s to %d',str,max(N)); end;
            group_names{i} = str;
          end;
        end;
      else
        group_names = [];
      end;
      figxlabel = sprintf('%s effect size (Cohen''s d)',regexprep(meas,'_',' '));
      xlim = [-2,2];
      color_order = {'b','g'};
      clf;
      mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
        'barwidth',parms.barwidth,...
        'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
        'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
        'color_order',color_order);
      set(gcf,'Visible','off');
      fstem = sprintf('%s/%s_%s_effsz',...
        parms.outdir,meas,outstr);
      if ~isempty(parms.outfix)
        fstem = [fstem '_' parms.outfix];
      end;
      fname_jpg = [fstem '.jpg'];
      set(gcf,'PaperUnits','inches',...
              'PaperPosition',[0 0 parms.paper_width parms.paper_height])
      print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
      if parms.plot_eps_flag
        fname_eps = [fstem '.eps'];
        print(gcf,'-deps',fname_eps);
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % effect sizes for comparisons between all groups
    if parms.plot_flag>=2
      data = [];
      err = []; % what are error bars for effect size?
      group_names = [];
      for i=1:length(fiberstats.between_group(1).TTEST)
        groupA = fiberstats.between_group(1).TTEST(i).groupA;
        groupB = fiberstats.between_group(1).TTEST(i).groupB;
        group_names{i} = [groupA ' VS ' groupB];
        tmp_vec = [];
        for j=1:length(ind_fibers)
          f = ind_fibers(j);
          tmp_vec = [tmp_vec,-fiberstats.between_group(f).TTEST(i).d];
        end;
        data = [data;tmp_vec];
      end;
      if ~isempty(data)
        data = data';
        if ~parms.legend_flag, group_names = []; end;
        figxlabel = sprintf('%s effect size (Cohen''s d)',regexprep(meas,'_',' '));
        xlim = [-2,2];
        color_order = {'b','g','r'};
        clf;
        mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
          'barwidth',parms.barwidth,'text_v_offset',parms.text_v_offset,...
          'text_v_scale',parms.text_v_scale,...
          'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
          'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
          'color_order',color_order);
        set(gcf,'Visible','off');
        fstem = sprintf('%s/%s_%s_effszB',...
          parms.outdir,meas,outstr);
        if ~isempty(parms.outfix)
          fstem = [fstem '_' parms.outfix];
        end;
        fname_jpg = [fstem '.jpg'];
        set(gcf,'PaperUnits','inches',...
                'PaperPosition',[0 0 parms.paper_width parms.paper_height])
        print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
        if parms.plot_eps_flag
          fname_eps = [fstem '.eps'];
          print(gcf,'-deps',fname_eps);
        end;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % p-values for comparisons between all groups
    if parms.plot_flag>=2
      data = [];
      err = []; % what are error bars for effect size?
      group_names = [];
      for i=1:length(fiberstats.between_group(1).TTEST)
        groupA = fiberstats.between_group(1).TTEST(i).groupA;
        groupB = fiberstats.between_group(1).TTEST(i).groupB;
        group_names{i} = [groupA ' VS ' groupB];
        tmp_vec = [];
        for j=1:length(ind_fibers)
          f = ind_fibers(j);
          tmp_vec = [tmp_vec,fiberstats.between_group(f).TTEST(i).plog];
        end;
        data = [data;tmp_vec];
      end;
      if ~isempty(data)
        data = data';
        if ~parms.legend_flag, group_names = []; end;
        figxlabel = sprintf('%s p value (-log10)',regexprep(meas,'_',' '));
        xlim = [1,5];
        color_order = {'b','g','r'};
        clf;
        mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
          'barwidth',parms.barwidth,'text_v_offset',parms.text_v_offset,...
          'text_v_scale',parms.text_v_scale,...
          'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
          'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
          'color_order',color_order);
        set(gcf,'Visible','off');
        fstem = sprintf('%s/%s_%s_logpvalB',...
          parms.outdir,meas,outstr);
        if ~isempty(parms.outfix)
          fstem = [fstem '_' parms.outfix];
        end;
        fname_jpg = [fstem '.jpg'];
        set(gcf,'PaperUnits','inches',...
                'PaperPosition',[0 0 parms.paper_width parms.paper_height])
        print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
        if parms.plot_eps_flag
          fname_eps = [fstem '.eps'];
          print(gcf,'-deps',fname_eps);
        end;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANOVA p-values
    data = [];
    for j=1:length(ind_fibers)
      f = ind_fibers(j);
      data = [data;fiberstats.between_group(f).ANOVA.plog];
    end;
    err = [];
    group_names = [];
    figxlabel = sprintf('%s ANOVA p value (-log10)',regexprep(meas,'_',' '));
    xlim = [1 5];
    color_order = {'k'};
    clf;
    mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
      'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
      'text_v_scale',parms.text_v_scale,...
      'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
      'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
      'color_order',color_order);
    set(gcf,'Visible','off');
    fstem = sprintf('%s/%s_%s_ANOVA_logpval',...
      parms.outdir,meas,outstr);
    if ~isempty(parms.outfix)
      fstem = [fstem '_' parms.outfix];
    end;
    fname_jpg = [fstem '.jpg'];
    set(gcf,'PaperUnits','inches',...
            'PaperPosition',[0 0 parms.paper_width parms.paper_height])
    print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
    if parms.plot_eps_flag
      fname_eps = [fstem '.eps'];
      print(gcf,'-deps',fname_eps);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANOVA F-stats
    if parms.plot_flag>=2
      data = [];
      for j=1:length(ind_fibers)
        f = ind_fibers(j);
        data = [data;fiberstats.between_group(f).ANOVA.F];
      end;
      err = [];
      group_names = [];
      figxlabel = sprintf('%s ANOVA F',regexprep(meas,'_',' '));
      xlim = [0 15];
      color_order = {'k'};
      clf;
      mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
        'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
        'text_v_scale',parms.text_v_scale,...
        'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
        'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
        'color_order',color_order);
      set(gcf,'Visible','off');
      fstem = sprintf('%s/%s_%s_ANOVA_F',...
        parms.outdir,meas,outstr);
      if ~isempty(parms.outfix)
        fstem = [fstem '_' parms.outfix];
      end;
      fname_jpg = [fstem '.jpg'];
      set(gcf,'PaperUnits','inches',...
              'PaperPosition',[0 0 parms.paper_width parms.paper_height])
      print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
      if parms.plot_eps_flag
        fname_eps = [fstem '.eps'];
        print(gcf,'-deps',fname_eps);
      end;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANOVA eta2
    if parms.plot_flag>=2
      data = [];
      for j=1:length(ind_fibers)
        f = ind_fibers(j);
        data = [data;fiberstats.between_group(f).ANOVA.eta2];
      end;
      err = [];
      group_names = [];
      figxlabel = sprintf('%s ANOVA eta squared',regexprep(meas,'_',' '));
      xlim = [0 0.5];
      color_order = {'k'};
      clf;
      mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
        'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
        'text_v_scale',parms.text_v_scale,...
        'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
        'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
        'color_order',color_order);
      set(gcf,'Visible','off');
      fstem = sprintf('%s/%s_%s_ANOVA_eta2',...
        parms.outdir,meas,outstr);
      if ~isempty(parms.outfix)
        fstem = [fstem '_' parms.outfix];
      end;
      fname_jpg = [fstem '.jpg'];
      set(gcf,'PaperUnits','inches',...
              'PaperPosition',[0 0 parms.paper_width parms.paper_height])
      print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
      if parms.plot_eps_flag
        fname_eps = [fstem '.eps'];
        print(gcf,'-deps',fname_eps);
      end;
    end;

    if parms.age_flag
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ANCOVA p-values
      data = [];
      for j=1:length(ind_fibers)
        f = ind_fibers(j);
        data = [data;fiberstats.between_group(f).ANCOVA.plog];
      end;
      err = [];
      group_names = [];
      figxlabel = sprintf('%s ANCOVA p value (-log10)',regexprep(meas,'_',' '));
      xlim = [1 5];
      color_order = {'k'};
      clf;
      mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
        'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
        'text_v_scale',parms.text_v_scale,...
        'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
        'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
        'color_order',color_order);
      set(gcf,'Visible','off');
      fstem = sprintf('%s/%s_%s_ANCOVA_logpval',...
        parms.outdir,meas,outstr);
      if ~isempty(parms.outfix)
        fstem = [fstem '_' parms.outfix];
      end;
      fname_jpg = [fstem '.jpg'];
      set(gcf,'PaperUnits','inches',...
              'PaperPosition',[0 0 parms.paper_width parms.paper_height])
      print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
      if parms.plot_eps_flag
        fname_eps = [fstem '.eps'];
        print(gcf,'-deps',fname_eps);
      end;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ANCOVA F-stats
      if parms.plot_flag>=2
        data = [];
        for j=1:length(ind_fibers)
          f = ind_fibers(j);
          data = [data;fiberstats.between_group(f).ANCOVA.F];
        end;
        err = [];
        group_names = [];
        figxlabel = sprintf('%s ANCOVA F',regexprep(meas,'_',' '));
        xlim = [0 15];
        color_order = {'k'};
        clf;
        mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
          'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
          'text_v_scale',parms.text_v_scale,...
          'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
          'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
          'color_order',color_order);
        set(gcf,'Visible','off');
        fstem = sprintf('%s/%s_%s_ANCOVA_F',...
          parms.outdir,meas,outstr);
        if ~isempty(parms.outfix)
          fstem = [fstem '_' parms.outfix];
        end;
        fname_jpg = [fstem '.jpg'];
        set(gcf,'PaperUnits','inches',...
                'PaperPosition',[0 0 parms.paper_width parms.paper_height])
        print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
        if parms.plot_eps_flag
          fname_eps = [fstem '.eps'];
          print(gcf,'-deps',fname_eps);
        end;
      end;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ANCOVA eta2
      if parms.plot_flag>=2
        data = [];
        for j=1:length(ind_fibers)
          f = ind_fibers(j);
          data = [data;fiberstats.between_group(f).ANCOVA.eta2];
        end;
        err = [];
        group_names = [];
        figxlabel = sprintf('%s ANCOVA eta squared',regexprep(meas,'_',' '));
        xlim = [0 0.5];
        color_order = {'k'};
        clf;
        mmil_barh(data,'err',err,'data_names',data_names,'group_names',group_names,...
          'barwidth',parms.barwidth/2,'text_v_offset',parms.text_v_offset,...
          'text_v_scale',parms.text_v_scale,...
          'fontsize',parms.fontsize,'linewidth',parms.linewidth,...
          'title',figtitle,'xlabel',figxlabel,'xlim',xlim,...
          'color_order',color_order);
        set(gcf,'Visible','off');
        fstem = sprintf('%s/%s_%s_ANCOVA_eta2',...
          parms.outdir,meas,outstr);
        if ~isempty(parms.outfix)
          fstem = [fstem '_' parms.outfix];
        end;
        fname_jpg = [fstem '.jpg'];
        set(gcf,'PaperUnits','inches',...
                'PaperPosition',[0 0 parms.paper_width parms.paper_height])
        print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
        if parms.plot_eps_flag
          fname_eps = [fstem '.eps'];
          print(gcf,'-deps',fname_eps);
        end;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scatter plots
    if parms.plot_flag>=3 && parms.age_flag
      fname_in = matfiles{m};
      clear group_data;
      load(fname_in);
      group_data = select_subjects(group_data,parms.age_flag);
      nsubs = length(group_data);
      if nsubs==0
        error('no valid subjects found');
      end;
      for f=1:size(fiberstats.within_group,1)
        fibername = fiberstats.within_group(f,1).fibername;
        included_subjects = fiberstats.subject_inclusion(f).included;
        all_subjects = {group_data.SubjID};
        [c,ia,ib] = intersect(included_subjects,all_subjects);
        tmp_group_data = group_data(ib);
        subj_groups = {tmp_group_data.Group};
        subj_ages = cell2mat({tmp_group_data.Age});
        mkr_list = {'ks','bh','gp','rx','cd','mo','y*'};
        clf; hold on;
        h = 1;
        fiber_list=[];
        for i=1:size(fiberstats.within_group,1)
          fiber_list{i}=tmp_group_data(1).fiber_data(i).fibernum;
        end;
        fiber_list=cell2mat(fiber_list);
        for g=1:length(groups)
          if h>length(mkr_list)
            h=1;
          end;
          mkr = mkr_list{h}; h=h+1;
          ind = find(strcmp(groups{g},subj_groups));
          fib_ind = find(fiber_list==fiberstats.within_group(f,1).fibernum);
          vals = [];
          for i=ind
            val = tmp_group_data(i).fiber_data(fib_ind).avg*scalefact;
            vals = [vals val];
          end;
          ages = subj_ages(ind);
          plot(ages,vals,mkr);
        end;
        xlabel('Age');
        ylabel(meas);
        if parms.title_flag, title(fibername); end;
        legend('ACE','control','patient');
        tmp_name = regexprep(fibername,'\s+','_');
        fstem = sprintf('%s/%s_%s_%s_scat',...
          parms.outdir,meas,outstr,tmp_name);
        if ~isempty(parms.outfix)
          fstem = [fstem '_' parms.outfix];
        end;
        fname_jpg = [fstem '.jpg'];
        set(gcf,'PaperUnits','inches',...
                'PaperPosition',[0 0 parms.paper_width parms.paper_height])
        print(gcf,'-djpeg',fname_jpg,sprintf('-r%d',parms.dpi));
        if parms.plot_eps_flag
          fname_eps = [fstem '.eps'];
          print(gcf,'-deps',fname_eps);
        end;
      end;
    end;
  end;
  close(gcf);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function group_data = select_subjects(group_data,age_flag)
  if ~isfield(group_data,'StatsFlag')
    for s=1:length(group_data)
      group_data.StatsFlag = 1;
    end;
  end;
  for s=1:length(group_data)
    if age_flag && (isempty(group_data(s).Age) || group_data(s).Age==0)
      group_data(s).StatsFlag = 0;
    end;
    if strcmp(group_data(s).VisitNumber,'post')
      group_data(s).StatsFlag = 0;
    end;
  end;
  group_data = group_data(find(cell2mat({group_data.StatsFlag})));
return

