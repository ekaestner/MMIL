function BOLD_MMIL_GLM_ROI_GroupAvg(ProjID,varargin)
%function BOLD_MMIL_GLM_ROI_GroupAvg(ProjID,varargin)
%
% Purpose: Calculate cross-subject average of GLM ROI results
%   Plot individual subject and group results
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters:
%  'StudyInfo': struct array containing info for each subject
%    Must have 'VisitID' to indicate name of the orig data directory
%   'GLM_ROI_VisitID' and 'STRUCT_VisitID' may also be supplied
%    May control which subjects are included in average with
%      field 'GroupAvg_GLM_ROI'
%    May control which scans are used for each subject with BOLDScanNums
%    If empty, will use ProjID to get StudyInfo
%    {default = []}
%  'RootDirs': struct containing locations of root data dirs
%    must include this field: proc_bold
%      (e.g. define in MMIL_ProjInfo.csv)
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%    {default = []}
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    if ProjID is empty, will be relative to current working directory
%    {default = 'GLM_ROI_GroupAvg'}
%  'outstem': output file stem (relative to outdir)
%    {default = 'GLM_ROI'}
%  'qcflag': [0|1] exclude subjects with StudyInfo.QC = 0
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Optional Parameters that specify how analysis was done
%  'infix': string attached to processed BOLD file names
%    {default = 'corr_resBOLD'}
%
% Optional Input from ProjInfo, command line, or StudyInfo
%  'GLM_ROI_outstem': string added to output file names (appended to GLM stem)
%    {default = []}
%  'GLM_ROI_stem': file stem for batch of labels (e.g. lh.{stem}_{name}.label)
%    {default = []}
%
% Optional Parameters controlling how plots look
%  'contrast_levels': vector of contrast levels (must match results)
%    {default = [0.05,0.3,1]}
%  'contrast_label': x-axis label
%    {default 'Relative Stimulus Contrast'}
%  'normflag': normalize amplitudes to value for max contrast level
%     {default = 0}
%  'stderrflag': for average plots, use standard error for error bars
%     otherwise use standard deviation
%     {default = 1}
%  'ylim': y-axis limits
%    {default = []}
%  'xlim': x-axis limits
%    {default = [0,1.05]}
%  'logx_flag': [0|1] use log scale for x-axis
%    {default = 0}
%  'eps_flag': save plot in eps format (vector graphics)
%    {default = 0}
%  'visible_flag': [0|1] whether to display plots to screen
%    (otherwise save only)
%    {default = 1}
%
% Created:  03/05/11 by Don Hagler
% Last Mod: 12/09/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms_filter = {...
  'StudyInfo',[],[],...
  'RootDirs',[],[],...
...
  'outdir','GLM_ROI_GroupAvg',[],...
  'outstem','GLM_ROI',[],...
  'infix','corr_resBOLD',[],...
  'qcflag',true,[false true],...
  'forceflag',false,[false true],...  
...
  'infix','corr_resBOLD',[],...
...
  'GLM_ROI_outstem',[],[],...
  'GLM_ROI_stem',[],[],...
...
  'contrast_levels',[0.05,0.3,1],[],...
  'contrast_label','Relative Stimulus Contrast',[],...
  'normflag',false,[false true],...
  'stderrflag',true,[false true],...
  'ylim',[],[],...
  'xlim',[0,1.05],[],...
  'logx_flag',false,[false true],...
  'eps_flag',false,[false true],...
  'visible_flag',true,[false true],...
  'plot_subs_flag',true,[false true],...
...
  'roi_colors',{'b','g','r'},[],...
  'roi_linewidths',[],[],...
  'linewidth',1.5,[],...
  'fontsize',12,[],...
  'fontname','Arial',[],...
  'label_flag',true,[false true],...
  'units','% signal change',[],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'fnamestem','BOLD',[],...
  'required_containers',{'proc_bold','fsurf'},[],...
  'avg_tag','GroupAvg_GLM_ROI',[],...
...
  'legend_loc','SouthEast',[],...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RootDirs,StudyInfo,parms] = check_input(ProjID,parms_filter,varargin);

fname_results = [parms.outstem '_GroupAvg_results.mat'];
if ~exist(fname_results,'file') || parms.forceflag
  results = [];
  results.roi_avg_vals = single(0);
  results.nrois = [];
  results.ncontrasts = [];
  results.subjdata = [];
  for s=1:length(StudyInfo)
    fprintf('%s: getting GLM_ROI results for VisitID %s...\n',...
      mfilename,StudyInfo(s).VisitID);

    % get GLM_ROI results for this session
    subjresults = load_results(parms,RootDirs,StudyInfo(s));
%    if isempty(subjresults), continue; end;
    if isempty(subjresults)
      fprintf('%s: ERROR: missing data for %s\n',...
        mfilename,StudyInfo(s).VisitID);
      return;
    end;

    % extract thresholded values
    subjdata = get_subjdata(parms,subjresults,StudyInfo(s));

    % check that subjects have matching number of ROIs and contrast levels
    if isempty(results.nrois)
      results.nrois = subjdata.nrois;
      results.roinames = subjdata.roinames;
      results.ncontrasts = subjdata.ncontrasts;
      if results.ncontrasts ~= length(parms.contrast_levels)
        fprintf('%s: WARNING: mismatch between number of contrast_levels (%d) and results ncontrasts (%d)\n',...
          mfilename,length(parms.contrasts_levels),results.ncontrasts);
        parms.contrast_levels = [1:results.ncontrasts];
        parms.xlim = [1,results.ncontrasts];
      end;
      if parms.logx_flag
        parms.contrast_levels = log10(parms.contrast_levels);
        parms.xlim = [min(parms.contrast_levels),0];
      end;
    elseif results.nrois ~= subjdata.nrois
      error('number of ROIs for %s does not match\n',...
        mfilename,StudyInfo(s).VisitID);
    elseif results.ncontrasts ~= subjdata.ncontrasts
      error('number of contrast levels for %s does not match\n',...
        mfilename,StudyInfo(s).VisitID);
    end;

    % combine into struct array
    results.subjdata = [results.subjdata subjdata];

    % plot contrast response functions
    if parms.plot_subs_flag
      toplabl = sprintf('Amplitude vs. Contrast for %s',...
        regexprep(subjdata.VisitID,'_',' '));
      plot_subj_response(parms,subjdata,toplabl);
    end;
  end;

  % compile values for all subjects
  results.N = length(results.subjdata);
  if results.N==0, return; end;
  results.roi_avg_vals = zeros(results.nrois,results.ncontrasts,results.N);
  for i=1:results.N
    results.roi_avg_vals(:,:,i) = results.subjdata(i).roi_avg_vals;
  end;

  toplabl = 'Amplitude vs. Contrast for Multiple Subjects';
  plot_multisubj_response(parms,results,toplabl);

  toplabl = 'Amplitude vs. Contrast Group Average';
  plot_groupavg_response(parms,results,toplabl);
  
  save(fname_results,'results');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RootDirs,StudyInfo,parms] = check_input(ProjID,parms_filter,options)
  parms = mmil_args2parms(options,parms_filter);
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if ~isempty(ProjInfo)
    % For arg names present in both varargin and ProjInfo
    % the varargin values will appear in merged_args
    ProjInfo_args = MMIL_Args(ProjInfo,mfilename);
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  nsubs  = length(StudyInfo);
  if ~nsubs, error('no valid subjects in StudyInfo'); end;

  if isfield(StudyInfo,parms.avg_tag)
    for i=1:nsubs
      if isempty(StudyInfo(i).(parms.avg_tag))
        StudyInfo(i).(parms.avg_tag) = 0;
      end;
    end;
    ind_GroupAvg = find([StudyInfo.(parms.avg_tag)]);
    StudyInfo = StudyInfo(ind_GroupAvg);
    nsubs  = length(StudyInfo);
    if ~nsubs
      error('no valid subjects with StudyInfo.GroupAvg=1');
    end;
  end;

  for i=1:nsubs
    if ~isfield(StudyInfo,'GLM_ROI_stem') || isempty(StudyInfo(i).GLM_ROI_stem)
      StudyInfo(i).GLM_ROI_stem = parms.GLM_ROI_stem;
    end;
    if ~isfield(StudyInfo,'GLM_ROI_outstem') || isempty(StudyInfo(i).GLM_ROI_outstem)
      StudyInfo(i).GLM_ROI_outstem = parms.GLM_ROI_outstem;
    end;
  end;

  % change outstem to be full path
  if mmil_isrelative(parms.outdir)
    if isempty(ProjID)
      parms.outdir = [pwd '/' parms.outdir];
    else
      parms.outdir = [RootDirs.home '/MetaData/' ProjID '/' parms.outdir];
    end;
  end;
  mmil_mkdir(parms.outdir);
  parms.outstem = [parms.outdir '/' parms.outstem];

  [tmp,parms.ind_max_cont] = max(parms.contrast_levels);
  
  parms.fontargs = {'FontSize',parms.fontsize,'FontName',parms.fontname};
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subjdata = get_subjdata(parms,results,StudyInfo)
  subjdata = StudyInfo;
  subjdata.VisitID = StudyInfo.VisitID;
  subjdata.nrois = results.nroi;
  subjdata.ROI_data = results.ROI_data;
  subjdata.ncontrasts = results.nconds;
  subjdata.roinames = {subjdata.ROI_data.name};
  subjdata.roi_avg_vals = results.resp';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = load_results(parms,RootDirs,StudyInfo)
  results = [];
  % set GLM analysis file stem
  ContainerPath = [RootDirs.proc_bold '/' StudyInfo.proc_bold];
  snums_valid = mmil_getfield(StudyInfo,'BOLDScanNums',[]);
  [stem,errcode] = BOLD_MMIL_Set_GLM_Stem(ContainerPath,...
    'snums',snums_valid,'infix',parms.infix);
  if errcode, error('failed to set GLM stem for %s',ContainerPath); end;
  fname = stem;
  if ~isempty(StudyInfo.GLM_ROI_outstem)
    fname = [fname '_' StudyInfo.GLM_ROI_outstem];
  end;
  if ~isempty(StudyInfo.GLM_ROI_stem)
    fname = [fname '_' StudyInfo.GLM_ROI_stem];
  end;
  fname = [fname '_roi_data.mat'];
  if ~exist(fname,'file')
    fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
    return;
  end;
  load(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_labels(parms,roinames,toplabl)
  if ~isempty(parms.ylim), set(gca,'YLim',parms.ylim); end;
  set(gca,'XLim',parms.xlim);
  set(gca,parms.fontargs{:});
  xlabel(parms.contrast_label,parms.fontargs{:});
  ylabl = 'BOLD response';
  if parms.normflag
    ylabl = ['Relative ' ylabl];
  else
    ylabl = [ylabl ' (' parms.units ')'];
  end;
  ylabel(ylabl,parms.fontargs{:});
  title(toplabl,parms.fontargs{:});
  if parms.logx_flag
    tmp = get(gca);
    tmp_ticks = tmp.XTickLabel;
    for i=1:size(tmp_ticks,1)
      tmp_str = tmp_ticks(i,:);
      tmp_val = str2num(tmp_str);
      tmp_val = 10^tmp_val;
      tmp_str = num2str(tmp_val,'%0.2f');
      tmp_ticks(i,1:length(tmp_str)) = tmp_str;
    end;
    set(gca,'XTickLabel',tmp_ticks);
  end;
  legstrs = upper(regexprep(roinames,'_',' '));
  legend(legstrs,parms.fontargs{:},'Location',parms.legend_loc);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms,outstem)
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  print(gcf,'-dtiff',[outstem '.tif']);
  if parms.eps_flag, mmil_printeps(gcf,[outstem '.eps']); end;
  if ~parms.visible_flag, close(gcf); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_subj_response(parms,subjdata,toplabl)
  figure; clf; hold on;
  for r=1:subjdata.nrois
    roi_num = 1 + rem(r-1,length(parms.roi_colors));
    roi_color = parms.roi_colors{roi_num};
    if isempty(parms.roi_linewidths)
      roi_linewidth = parms.linewidth;
    else
      roi_num = 1 + rem(r-1,length(parms.roi_linewidths));
      roi_linewidth = parms.roi_linewidths(roi_num);
    end;    
    
    resps = subjdata.roi_avg_vals(r,:);
    if parms.normflag
      resps = resps/resps(parms.ind_max_cont);
    end;
    % plot amplitudes vs. contrasts
    if roi_linewidth>0
      plot(parms.contrast_levels,resps,'o-',...
        'Color',roi_color,'LineWidth',roi_linewidth);
    else
      plot(parms.contrast_levels,resps,'o','Color',roi_color);
    end;
  end;
  plot_labels(parms,subjdata.roinames,toplabl)
  outstem = [parms.outstem '_subj_' subjdata.VisitID];
  if parms.normflag, outstem = [outstem '_norm']; end;
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multisubj_response(parms,results,toplabl)
  figure; clf; hold on;
  for n=1:results.N
    for r=1:results.nrois
      roi_num = 1 + rem(r-1,length(parms.roi_colors));
      roi_color = parms.roi_colors{roi_num};
      if isempty(parms.roi_linewidths)
        roi_linewidth = parms.linewidth;
      else
        roi_num = 1 + rem(r-1,length(parms.roi_linewidths));
        roi_linewidth = parms.roi_linewidths(roi_num);
      end;    

      resps = results.subjdata(n).roi_avg_vals(r,:);
      if parms.normflag
        resps = resps/resps(parms.ind_max_cont);
      end;
      % plot amplitudes vs. contrasts
      if roi_linewidth>0
        plot(parms.contrast_levels,resps,'o-',...
          'Color',roi_color,'LineWidth',roi_linewidth);
      else
        plot(parms.contrast_levels,resps,'o','Color',roi_color);
      end;
    end;
  end;
  plot_labels(parms,results.roinames,toplabl)
  outstem = [parms.outstem '_multisubj'];
  if parms.normflag, outstem = [outstem '_norm']; end;
  save_plot(parms,outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_groupavg_response(parms,results,toplabl)
  figure; clf; hold on;
  for r=1:results.nrois
    roi_num = 1 + rem(r-1,length(parms.roi_colors));
    roi_color = parms.roi_colors{roi_num};
    if isempty(parms.roi_linewidths)
      roi_linewidth = parms.linewidth;
    else
      roi_num = 1 + rem(r-1,length(parms.roi_linewidths));
      roi_linewidth = parms.roi_linewidths(roi_num);
    end;    
  
    resps = squeeze(results.roi_avg_vals(r,:,:));
    if length(parms.contrast_levels)==1, resps = resps'; end;
    ind_good = find(max(isnan(resps),[],1)==0);
    nbad = results.N - length(ind_good);
    if nbad>0
      fprintf('%s: WARNING: missing %s value for %d/%d subjects\n',...
        mfilename,results.roinames{r},nbad,results.N);
    end;
    if isempty(ind_good), continue; end;
    resps = resps(:,ind_good);
    if parms.normflag
      normvals = ones(size(resps,1),1)*resps(parms.ind_max_cont,:);
      resps = resps./ normvals;
    end;
    mean_resps = mean(resps,2);
    std_resps = std(resps,[],2);
    if parms.stderrflag
      std_resps = std_resps/sqrt(results.N);
    end;
    % plot amplitudes vs. contrasts
    if roi_linewidth>0
      errorbar(parms.contrast_levels,mean_resps,std_resps,'o-',...
        'Color',roi_color,'LineWidth',roi_linewidth);
    else
      errorbar(parms.contrast_levels,mean_resps,std_resps,'o',...
        'Color',roi_color);
    end;
  end;
  plot_labels(parms,results.roinames,toplabl)
  outstem = [parms.outstem '_groupavg'];
  if parms.normflag, outstem = [outstem '_norm']; end;
  save_plot(parms,outstem);
return;

