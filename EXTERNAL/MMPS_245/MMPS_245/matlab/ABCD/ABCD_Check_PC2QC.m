function abcd_check_pc2qc(varargin)
%function abcd_check_pc2qc()
%
% Optional input:
%  'indir': input metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'qcdir': qc directory
%     {default = ''}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'fname_projinfo': file name of whole project info
%     {default = []}
%  'infix': file suffix of input file
%     containing info about classified incoming info
%     {default = 'merged_pcqcinfo'}
%  'infix_missing': file suffix of proc missing info
%     {default = 'proc_missingfiles_list'}
%  'outfix': file suffix of output file
%     {default = 'pc2qc'}
%  'combined_outfix': file suffix of combined eprime info
%     {default = 'combined_eprime'}
%  'targetseries': targeted servies list
%     {default = 'T1','T2','dMRI','dMRI_FM','dMRI_FM_AP','dMRI_FM_PA','fMRI_FM','fMRI_FM_AP','fMRI_FM_PA','fMRI_MID_task','fMRI_SST_task','fMRI_nBack_task','rsfMRI','combined_eprime'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/30/17 by Feng Xue
% Last Mod: 08/31/17 by Feng Xue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% procfail = 1 waitautoqc=0, needrawqc=0
% procfail = 0 waitautoqc=1, needrawqc=0
% exclude undefined seriestypes
% sep 1st?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check input parameters
parms = check_input(varargin);
%check lock files
fname_lck=sprintf('%s/MetaData/%s/.check_pc2qc.lck',getenv('HOME'),parms.outstem);
if exist(fname_lck,'file')
  fprintf('%s\n','lock files exist!.');
  return;
end

%Place lock file
fclose(fopen(fname_lck, 'w'));
while exist(sprintf('%s/.merged_pcqcinfo.lck',parms.outdir),'file')
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;


fname_series_info = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.infix);
fname_missingfiles_info = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.infix_missing);
fname_pc2qc_series = sprintf('%s/%s_%s_series.csv',parms.outdir,parms.outstem,parms.outfix);
fname_pc2qc_event = sprintf('%s/%s_%s_event.csv',parms.outdir,parms.outstem,parms.outfix);
if ~exist(fname_series_info,'file'), error('file %s not found',fname_series_info); end;
if ~exist(fname_missingfiles_info,'file'), error('file %s not found',fname_missingfiles_info); end;


series_info = abcd_load_csv(fname_series_info);
missingfiles_info = abcd_load_csv(fname_missingfiles_info);

%data reduction
ind_valid = find(~cellfun(@isempty,{series_info.ABCD_Compliant})); 
series_info = series_info(ind_valid);    

ind_valid = find(ismember({series_info.ABCD_Compliant},{'Yes','No'}));
series_info = series_info(ind_valid);

%filter all undefined series
ind_valid = cellfun('isempty',regexp({series_info.SeriesType},'Undefined'));
series_info = series_info(ind_valid);

%Completed
%series_info=series_info(find( cellfun(@(x)isequal(x,1),{series_info.Completed}) ));


[series_info_tmp(1:numel(series_info)).pGUID] = deal(series_info.pGUID);
[series_info_tmp(1:numel(series_info)).VisitID] = deal(series_info.VisitID);
[series_info_tmp(1:numel(series_info)).EventName] = deal(series_info.EventName);
[series_info_tmp(1:numel(series_info)).site] = deal(series_info.SiteName);
[series_info_tmp(1:numel(series_info)).SessionType] = deal(series_info.SessionType);
[series_info_tmp(1:numel(series_info)).SeriesType] = deal(series_info.SeriesType);
[series_info_tmp(1:numel(series_info)).ABCD_Compliant] = deal(series_info.ABCD_Compliant);
[series_info_tmp(1:numel(series_info)).SeriesDescription] = deal(series_info.SeriesDescription);
[series_info_tmp(1:numel(series_info)).Completed] = deal(series_info.Completed);
[series_info_tmp(1:numel(series_info)).fname_json] = deal(series_info.fname_json);
[series_info_tmp(1:numel(series_info)).fname_pc_json] = deal(series_info.fname_pc_json);
[series_info_tmp(1:numel(series_info)).fstem] = deal(series_info.fstem);
[series_info_tmp(1:numel(series_info)).nrev] = deal(series_info.nrev);
[series_info_tmp(1:numel(series_info)).revdisp] = deal(series_info.revdisp);
[series_info_tmp(1:numel(series_info)).QC] = deal(series_info.QC);

series_info=series_info_tmp';
clear series_info_tmp;
series_info = set_pguidevents(series_info);

series_info = summarize_byseries(series_info,'revdisp','gt(x,0)','iqc_pc2qc_revdisagree',1);

series_info = summarize_byseries(series_info,'nrev','gt(x,1)','iqc_pc2qc_multirev',1);
series_info = summarize_byseries(series_info,'iqc_pc2qc_revdisagree','gt(x,0)','iqc_pc2qc_multirev',[]);

series_info = summarize_byseries(series_info,'nrev','isempty','iqc_pc2qc_waitautoqc',1);

%series_info = find_autoqcfailure(series_info);
ind_procfail = find(ismember({series_info.VisitID},{missingfiles_info.VisitID}));
[series_info(ind_procfail).iqc_pc2qc_procfail] = deal(1);

[series_info(ind_procfail).iqc_pc2qc_waitautoqc] = deal(0);
[series_info(ind_procfail).iqc_pc2qc_revdisagree] = deal(0);
[series_info(ind_procfail).iqc_pc2qc_multirev] = deal(0);
series_info = find_needrawqc(series_info);
series_info = mmil_sortstruct(series_info,{'site','pGUID','EventName'});
% write output file
mmil_struct2csv(series_info,fname_pc2qc_series);

event_info = check_event(series_info,parms);
event_info = mmil_sortstruct(event_info,{'site','id_redcap','redcap_event_name'});
% write output file
mmil_struct2csv(event_info,fname_pc2qc_event);

%delete lock file
delete(fname_lck);

return;
%%%%%%
function series_info = summarize_byseries(series_info,field_source,formula,field_target,value_target);
  [~,idx_valid] = find_real_id(series_info,field_source,formula);
  if isempty(value_target)
    eval(sprintf('[series_info(idx_valid).%s] = deal([]);',field_target));
  else
    eval(sprintf('[series_info(idx_valid).%s] = deal(%d);',field_target,value_target));
  end
return

function [idx_empty,idx_valid] = find_real_id(series_info,field_source,formula,varargin)
  switch formula
    case 'isempty'
      eval(sprintf('idx_valid = cellfun(@isempty,{series_info.%s});',field_source));
      idx_valid = find(idx_valid);
      idx_empty = idx_valid;
    case 'notisempty'
      eval(sprintf('idx = ~cellfun(@isempty,{series_info.%s});',field_source));
      idx_valid = find(idx);
      idx_empty = find(~idx);
    otherwise
      eval(sprintf('idx_valid = ~cellfun(@isempty,{series_info.%s});',field_source));
      idx_empty = find(~idx_valid);
      tmp = series_info(idx_valid);
      idx_real = find(idx_valid);
      eval(sprintf('idx_valid = cellfun(@(x) %s,{tmp.%s});',formula,field_source));
      idx_valid = idx_real(idx_valid);
      if ~isempty(varargin)
        matrix_diff = varargin{1};
        if isempty(find(matrix_diff>1)), matrix_diff = find(matrix_diff); end;
        idx_valid = setdiff(idx_valid,matrix_diff);
      end
  end
return

function series_info = find_autoqcfailure(series_info)
%ind_completed = cellfun(@(x) isequal(x,1),{series_info.Completed});
%ind = find(ind_completed==1);
%for i=1:length(ind)
for i=1:length(series_info)
  %VisitID = series_info(ind(i)).VisitID;
  %fname_pc_json = series_info(ind(i)).fname_pc_json;
  VisitID = series_info(i).VisitID;
  fname_pc_json = series_info(i).fname_pc_json;
  parentroot = fname_pc_json(1,1:strfind(fname_pc_json,VisitID)-4);
  %SeriesType = series_info(ind(i)).SeriesType;
  SeriesType = series_info(i).SeriesType;
  switch SeriesType(1)
    case 'T'
      proc_root=sprintf('%s/proc',parentroot);
    case {'f','r'}
      proc_root=sprintf('%s/proc_bold',parentroot);
    case 'd'
      proc_root=sprintf('%s/proc_dti',parentroot);
    otherwise
      fprintf('Wrong SeriesType of %s',SeriesType);
      continue;
  end
  proc_base_pattern=sprintf('%s/*%s*',proc_root,VisitID);
  proc_base = dir(proc_base_pattern);
  if isempty(proc_base)
    fprintf('Can not find proc dir: %s',proc_base_pattern);
    continue;
  else
    proc_dir = sprintf('%s/%s',proc_root,proc_base.name);
    missing_files = sprintf('%s/*_missing_files.txt',proc_dir);
    %if length(dir(missing_files)) ==1, series_info(ind(i)).iqc_pc2qc_procfail = 1; end;
    if length(dir(missing_files)) ==1, series_info(i).iqc_pc2qc_procfail = 1; end;
  end
end
return

function series_info = find_needrawqc(series_info)

%ind_needrawqc = cell2mat(cellfun(@(x) lt(x,2),{series_info.nrev},'un',false));
%idx_needrawqc = find(ind_needrawqc == 1);

[~,idx_needrawqc] = find_real_id(series_info,'nrev','lt(x,2)');

%remove procfail
%ind_procfail = ~cellfun(@(x) isequal(x,1),{series_info.iqc_pc2qc_procfail},'un',false);
%idx_procfail = find(ind_procfail == 1);
[~,idx_procfail] = find_real_id(series_info,'iqc_pc2qc_procfail','isequal(x,1)');

%remove waitautoqc
%ind_waitautoqc = ~cellfun(@(x) isequal(x,1),{series_info.iqc_pc2qc_waitautoqc},'un',false);
%idx_waitautoqc = find(ind_waitautoqc == 1);
[~,idx_waitautoqc] = find_real_id(series_info,'iqc_pc2qc_waitautoqc','isequal(x,1)');

idx_needrawqc = setdiff(idx_needrawqc,idx_procfail);
idx_needrawqc = setdiff(idx_needrawqc,idx_waitautoqc);

[series_info(idx_needrawqc).iqc_pc2qc_needrawqc] = deal(1);
return

function event_info = check_event(series_info,parms)
  fprintf('%s: summarizing for each event...\n',mfilename);
  event_info = [];
  [pguidevent_uniq,ind_uniq] = unique({series_info.pguidevent},'first');
  for i=1:length(ind_uniq)
    j = ind_uniq(i);
    % copy event-specific info
    event_info(i).pguidevent = pguidevent_uniq{i};
    event_info(i).id_redcap = series_info(j).pGUID;
    event_info(i).redcap_event_name = series_info(j).EventName;
    event_info(i).site = series_info(j).site;
    % check QC
    ind_series = find(strcmp({series_info.pguidevent},pguidevent_uniq{i}));
    event_data = series_info(ind_series);

    %idx = cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_revdisagree},'un',false);
    [~,idx] = find_real_id(event_data,'iqc_pc2qc_revdisagree','isequal(x,1)');
    if isempty(idx)
      event_info(i).iqc_event_revdisagree = 0;
    else
      event_info(i).iqc_event_revdisagree = 1;
    end;

    %idx = cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_multirev},'un',false);
    [~,idx] = find_real_id(event_data,'iqc_pc2qc_multirev','isequal(x,1)');
    if length(idx) < length(event_data)
      event_info(i).iqc_event_multirev = 0;
    else
      event_info(i).iqc_event_multirev = 1;
    end;

    %idx = cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_waitautoqc},'un',false);
    [~,idx] = find_real_id(event_data,'iqc_pc2qc_waitautoqc','isequal(x,1)');
    if isempty(idx)
      event_info(i).iqc_event_waitautoqc = 0;
    else
      event_info(i).iqc_event_waitautoqc = 1;
      event_info = find_series_by_seriestype(event_info,i,'waitautoqc',event_data,parms);
    end;

    %idx = cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_needrawqc},'un',false);
    [~,idx] = find_real_id(event_data,'iqc_pc2qc_needrawqc','isequal(x,1)');
    if isempty(idx)
      event_info(i).iqc_event_needrawqc = 0;
    else
      event_info(i).iqc_event_needrawqc = 1;
      event_info = find_series_by_seriestype(event_info,i,'needrawqc',event_data,parms);
    end;

    %idx = cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_procfail},'un',false);
    [~,idx] = find_real_id(event_data,'iqc_pc2qc_procfail','isequal(x,1)');
    if isempty(idx)
      event_info(i).iqc_event_procfail = 0;
    else
      event_info(i).iqc_event_procfail = 1;
      event_info = find_series_by_seriestype(event_info,i,'procfail',event_data,parms);
    end;
    %idx = cellfun(@(x) isequal(x,1),{event_data.Completed},'un',false);
    [~,idx] = find_real_id(event_data,'Completed','isequal(x,1)',idx);
    event_data = event_data(idx);

    %idx = ~cellfun(@(x) isequal(x,1),{event_data.iqc_pc2qc_procfail},'un',false);
    %event_data = event_data(idx);

    %idx = cell2mat(cellfun(@(x) gt(x,0),{event_data.QC},'un',false));

%    %Special treatment for data collected in 2016 that are completed but with procfail.
%    if strcmp(pguidevent_uniq{i}, 'NDAR_INVCBV85UMM_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVVAMKAM75_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVHZL15962_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVT185D4UD_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INV64GAJKPV_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVC96XA4XL_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVDE8G5XB3_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVLKWJTNDZ_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INV5WZDJT4C_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVL1FL82DL_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INV9BPV6EYL_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVJ2CBCRYX_baseline_year_1_arm_1') || ...
%       strcmp(pguidevent_uniq{i}, 'NDAR_INVT2JJ42KE_baseline_year_1_arm_1')
%
%      [~,idx] = find_real_id(event_data,'iqc_pc2qc_procfail','isequal(x,0)');
%      event_data_completed = event_data(idx);
%    else
%      event_data_completed = event_data;
%    end

%%%%%%
%As on 08/25/17 we are excluding all procfail series
%%%%%%
      [~,idx] = find_real_id(event_data,'iqc_pc2qc_procfail','isequal(x,0)');
      event_data_completed = event_data(idx);


    %[~,idx] = find_real_id(event_data_completed,'QC','gt(x,0)');
    [~,idx] = find_real_id(event_data_completed,'QC','notisempty');
    if length(idx) < length(event_data_completed)
      event_info(i).iqc_event_qcd = 0;
    else
      event_info(i).iqc_event_qcd = 1;
    end;

    idx = ~cellfun(@isempty, regexp({event_data.SeriesType},'rsfMRI'));
    event_data = event_data(idx); %So we have all completed rsfmri records
    if isempty(event_data)
      event_info(i).iqc_event_rsfmri_autoqcd = 0;
    else
      event_info(i).iqc_event_rsfmri_autoqcd = 1;
      for rsid = 1:length(event_data)
        if isequal(event_data(rsid).iqc_pc2qc_waitautoqc,1) || isequal(event_data(rsid).iqc_pc2qc_procfail,1)
          event_info(i).iqc_event_rsfmri_autoqcd = 0;
        end
      end
    end

    %if regexp(event_info(i).pguidevent,'NDAR_INV08P1JKNE'), keyboard; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'qcdir',[],[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'outfix','pc2qc',[],...
    'infix','merged_pcqcinfo',[],...
    'infix_missing','proc_missingfiles_list',[],...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'fname_projinfo',[],[],...
    'forceflag',true,[false true],...
    'targetseries',{'T1','T2','dMRI','dMRI_FM','dMRI_FM_AP','dMRI_FM_PA','fMRI_FM','fMRI_FM_AP','fMRI_FM_PA','fMRI_MID_task','fMRI_SST_task','fMRI_nBack_task','rsfMRI'},[]...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;
  if parms.indir(1) ~= '/', parms.indir= parms.outdir; end;
  if isempty(parms.fname_projinfo), parms.fname_projinfo = sprintf('%s/ProjInfo/MMIL_ProjInfo_all.csv',getenv('HOME')); end;
  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;
  if isempty(parms.qcdir)
    % load project info
    projinfo = abcd_load_projinfo_all(parms);
    parms.qcdir = projinfo(end).qc; 
    clear projinfo;
  end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function series_info = set_pguidevents(series_info)
  if isfield(series_info,'pGUID')
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.pGUID},{series_info.EventName},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  else
    pguidevents = cellfun(@(x,y) sprintf('%s_%s',x,y),...
                    {series_info.id_redcap},{series_info.redcap_event_name},...
                    'UniformOutput',false);
    [series_info.pguidevent] = deal(pguidevents{:});
  end
return;

function event_info = find_series_by_seriestype(event_info,index,varname,series_info,parms)
  for i=1:length(parms.targetseries)
    idx = strcmp({series_info.SeriesType},parms.targetseries(i));
    if find(idx)
      eval(sprintf('[~,idx_valid] = find_real_id(series_info(idx),''iqc_pc2qc_%s'',''isequal(x,1)'');',varname));
      if ~isempty(idx_valid)
        eval(sprintf('event_info(index).iqc_%s_%s = 1;',regexprep(lower(parms.targetseries{i}),'fmri_(.*)_task','$1'),varname));
      else
        eval(sprintf('event_info(index).iqc_%s_%s = 0;',regexprep(lower(parms.targetseries{i}),'fmri_(.*)_task','$1'),varname));
      end;
    end;
  end
return
