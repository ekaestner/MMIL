function abcd_check_incoming(varargin)
%function abcd_check_incoming(varargin)
%1. remove duplicate records in incoming_info from all qc accounts
%2. combine incoming_info from all accounts and save a combined one.
%3. generate classified seriestype and save as combined_incoming_info_classified
%
% Optional input:
%   'outdir': output directory
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'fname_info_combined': output file name of combined info
%     if empty, will use outdir and outfix to construct
%     {default = []}
%   'fname_info_kspace_combined': output file name of combined kspace info
%     {default = []}
%   'fname_info_missingtgz_combined': output file name of combined missingtgz info
%     {default = []}
%   'fname_info_classified': output file name of classified info
%     {default = []}
%   'fname_projinfo': file name of whole project info
%     {default = []}
%   'site': run this script on site(s) only
%     runs on all sites if empty
%     {default = []}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'infix': file suffix of input file
%     containing info about unique events in incoming directory
%     {default = 'incoming_info'}
%   'outfix : file suffix for combined info file
%     {default = 'incoming_info'}
%    'valid_types' all valid SeriesTypes
%      {default = {'t1','t2','dmri','dmri_fm','dmri_fm_ap','dmri_fm_pa','rsfmri',...
%                  'fmri_fm','fmri_fm_ap','fmri_fm_pa','mid','sst','nback'}
%   'forceflag': ignore existing results and run
%     {default = 1}
%
% Created:  04/26/17 by Feng Xue
% Last Mod: 11/01/17 by Feng Xue
%

  parms=check_input(varargin);

%  %check lock files
%  fname_lck=sprintf('%s/MetaData/%s/.check_incoming.lck',getenv('HOME'),parms.outstem);
%  if exist(fname_lck,'file')
%    fprintf('%s\n','lock files exist!.');
%    return;
%  end
%  %Place lock file
%  fclose(fopen(fname_lck, 'w'));

  projinfo = abcd_load_projinfo_all(parms);
  warning off;
  delete(parms.fname_info_missingtgz_combined);
  warning on;

  incoming_info = [];
  for i=1:length(projinfo)
%    while exist(sprintf('/home/%s/MetaData/%s/.incoming_info.lck',projinfo(i).account,parms.instem),'file')
%      fprintf('%s\n','lock files exist, waiting for previous process to finish.');
%      pause(30);
%    end;

    fname_incoming_info = sprintf('/home/%s/MetaData/%s/%s_%s.csv',...
      projinfo(i).account,parms.instem,parms.instem,parms.infix);
    if exist(fname_incoming_info,'file')
      incoming_info_tmp = mmil_csv2struct(fname_incoming_info);
      if ~isempty(incoming_info_tmp)
        [qcroot,~,~] = fileparts(projinfo(i).pc);
        incoming_info_tmp = addindex(incoming_info_tmp,i,qcroot);
        incoming_info = [incoming_info;incoming_info_tmp];
      end
    else
      fprintf('info file %s not found',fname_incoming_info)
    end
  end

  incoming_info = mmil_sortstruct(incoming_info,{'jsonissue','iqc_tgz_missing','sourceID'});
  incoming_info = remove_dup(incoming_info);
  incoming_info = mmil_sortstruct(incoming_info,{'site','id_redcap','redcap_event_name'});

  if ~isempty(incoming_info)
    mmil_struct2csv(incoming_info,parms.fname_info_combined);
    fname_cache = sprintf('%s/cache/%s_combined_%s_%s.mat',...
      parms.outdir,parms.outstem,parms.outfix,datestr(now,'yyyymmdd'));
    abcd_info = incoming_info;
    save(fname_cache,'abcd_info');

    incoming_info_classified = classifySeries(incoming_info,parms,projinfo);
    mmil_struct2csv(incoming_info_classified,parms.fname_info_classified);

    fname_cache = sprintf('%s/%s/%s_combined_%s_classified_%s.mat',...
      parms.outdir,'cache',parms.outstem,parms.outfix,datestr(now,'yyyymmdd'));
    abcd_info = incoming_info_classified;
    save(fname_cache,'abcd_info');
  end

  incoming_info_missingtgz=incoming_info(cellfun(@(x)isequal(x,1),{incoming_info.iqc_tgz_missing}));
  incoming_info_kspace = incoming_info(cellfun(@(x)isequal(x,1),{incoming_info.jsonissue}));
  for i=1:length(projinfo)
    incoming_info_tmp = incoming_info(cellfun(@(x)isequal(x,i),{incoming_info.sourceID}));
    fname_incoming_info_local = sprintf('%s/%s_%s_%s.csv',...
      parms.outdir,parms.outstem,parms.infix,projinfo(i).account);
    if ~isempty(incoming_info_tmp), mmil_struct2csv(incoming_info_tmp,fname_incoming_info_local); end;
    
  end
  if ~isempty(incoming_info_kspace), mmil_struct2csv(incoming_info_kspace,parms.fname_info_kspace_combined); end;
  if ~isempty(incoming_info_missingtgz), mmil_struct2csv(incoming_info_missingtgz,parms.fname_info_missingtgz_combined); end;
  
%  %delete lock file
%  delete(fname_lck);
return

function incoming_info = addindex(incoming_info,sourceID,qcroot)
%%TODO: sourceID can be used as weight
  [incoming_info.jsonissue] = deal(0);
  [incoming_info.sourceID] = deal(sourceID);
  [incoming_info.qcroot] = deal(qcroot);
  idx=cellfun(@(x,y) all([isempty(x),isempty(y)]),{incoming_info.ClassifyType},{incoming_info.SeriesDescription});
  [incoming_info(find(idx==1)).jsonissue] = deal(1);
return

function incoming_info = remove_dup(incoming_info)
  [idx idx]=unique(regexprep({incoming_info.json},'[Ss]ession[^_]*_',''),'last');
  incoming_info = incoming_info(idx);
return;


function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'fname_info_combined',[],[],...
    'fname_info_kspace_combined',[],[],...
    'fname_info_missingtgz_combined',[],[],...
    'fname_info_classified',[],[],...
    'fname_projinfo',[],[],...
    'site',[],[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'infix','incoming_info',[],...
    'outfix','incoming_info',[],...
    ...
    'forceflag',true,[false true],...
    ...
    'valid_types',{'t1','t2','dmri','dmri_fm','dmri_fm_ap','dmri_fm_pa','rsfmri',...
                   'fmri_fm','fmri_fm_ap','fmri_fm_pa','mid','sst','nback'},[],...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;

  if isempty(parms.fname_projinfo), parms.fname_projinfo = sprintf('%s/ProjInfo/MMIL_ProjInfo_all.csv',getenv('HOME')); end;
  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;


  if isempty(parms.fname_info_combined)
    parms.fname_info_combined = sprintf('%s/%s_combined_%s.csv',...
      parms.outdir,parms.outstem,parms.outfix);
  end;
  if isempty(parms.fname_info_kspace_combined)
    parms.fname_info_kspace_combined = sprintf('%s/%s_combined_%s_kspace.csv',...
      parms.outdir,parms.outstem,parms.outfix);
  end;
  if isempty(parms.fname_info_missingtgz_combined)
    parms.fname_info_missingtgz_combined = sprintf('%s/%s_combined_%s_missingtgz.csv',...
      parms.outdir,parms.outstem,parms.outfix);
  end;
  if isempty(parms.fname_info_classified)
    parms.fname_info_classified = sprintf('%s/%s_combined_%s_classified.csv',...
      parms.outdir,parms.outstem,parms.outfix);
  end;
return;


function incoming_info_classified = classifySeries(incoming_info,parms,projinfo)
  incoming_info_classified = [];
  incoming_info_classified_pos = 0;

  %Data Reduction
  idx=cellfun(@(x,y) all([isempty(x),isempty(y)]),{incoming_info.ClassifyType},{incoming_info.SeriesDescription});
  incoming_info=incoming_info(find(idx==0));

  idx1=cellfun(@(x) regexp(lower(char(x)),'patient|setter|<mpr collection>|loc'), {incoming_info.SeriesDescription}, 'UniformOutput',false);
  idx=cellfun(@isempty,idx1);
  incoming_info=incoming_info(idx);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Switch to use normalized structural images
  % 9/28/2017
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %idx1=cellfun(@(x) regexp(lower(char(x)),'norm|physio'), {incoming_info.ClassifyType}, 'UniformOutput',false);

  idx1=cellfun(@(x) regexp(lower(char(x)),'physio'), {incoming_info.ClassifyType}, 'UniformOutput',false);
  idx=cellfun(@isempty,idx1);
  incoming_info=incoming_info(idx);
%  %%%Siemens
%  idx_Siemens = strcmp({incoming_info.Manufacturer, 'SIEMENS');
%  incoming_info_Siemens = incoming_info(idx_Siemens);
%  idx = ~strcmp(upper({incoming_info.ClassifyType}), 'ABCD_T1');
%  incoming_info_Siemens = incoming_info_Siemens(idx);
%  idx = ~strcmp(upper({incoming_info.ClassifyType}), 'ABCD_T2');
%  incoming_info_Siemens = incoming_info_Siemens(idx);
%
%  idx_NonSiemens = ~idx_Siemens;
%  incoming_info_NonSiemens = incoming_info(idx_NonSiemens);
%  idx1=cellfun(@(x) regexp(lower(char(x)),'norm'), {incoming_info_NonSiemens.ClassifyType}, 'UniformOutput',false);
%  idx=cellfun(@isempty,idx1);
%  incoming_info_NonSiemens = incoming_info_NonSiemens(idx);
%
%  incoming_info = [incoming_info_Siemens;incoming_info_NonSiemens];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~isempty(parms.site)
    % only for designated site.
    idx = ~cellfun('isempty',strfind({incoming_info.site},parms.site));
    incoming_info=incoming_info(idx);
  end


  for i=1:length(incoming_info)
    incoming_info_classified(i).pguidevent=incoming_info(i).pguidevent;
    incoming_info_classified(i).id_redcap=incoming_info(i).id_redcap;
    incoming_info_classified(i).redcap_event_name=incoming_info(i).redcap_event_name;
    incoming_info_classified(i).site=incoming_info(i).site;
    incoming_info_classified(i).StudyDate=incoming_info(i).StudyDate;
    incoming_info_classified(i).StudyInstanceUID=incoming_info(i).StudyInstanceUID;
    incoming_info_classified(i).SeriesInstanceUID=incoming_info(i).SeriesInstanceUID;
    incoming_info_classified(i).SeriesTime=incoming_info(i).SeriesTime;
    incoming_info_classified(i).SeriesDescription=incoming_info(i).SeriesDescription;
    incoming_info_classified(i).SeriesType=get_SeriesType(incoming_info(i),projinfo);
    incoming_info_classified(i).json=incoming_info(i).json;
    incoming_info_classified(i).sourceID=incoming_info(i).sourceID;
  end;
return;

function SeriesType = get_SeriesType(sinfo,projinfo)
  if ~isfield(sinfo,'ClassifyType') || isempty(sinfo.ClassifyType)
    SeriesType = 'undefined';
    return;
  end

  if isfield(sinfo,'ImageType') && ~isempty(regexp(sinfo.ImageType,'DERIVED'))
      SeriesType = 'undefined';
      return;
  end

  SeriesType = regexprep(sinfo.ClassifyType,'-','_');
  SeriesType = lower(regexprep(SeriesType,'ABCD_',''));
  %%%%%%%%%%%%%%%%%
  % Switch to norm
  %%%%%%%%%%%%%%%%
  if strcmp(upper(sinfo.Manufacturer),'SIEMENS') && (strcmp(SeriesType,'t1') || strcmp(SeriesType,'t2')), SeriesType = strcat(SeriesType,'_raw'); return; end;
  SeriesType = regexprep(SeriesType,'_norm','');
  %%%%%%%%%%%%%%%

  SeriesType = regexprep(SeriesType,'diffusion','dmri');
  SeriesType = regexprep(SeriesType,'dti','dmri');
  SeriesType = regexprep(SeriesType,'sst_fmri','sst');
  SeriesType = regexprep(SeriesType,'nback_fmri','nback');
  SeriesType = regexprep(SeriesType,'(mid_fmri|fmri_mid)','mid');
%  if ~isempty(regexp(SeriesType,'norm')) ||...
%     ~isempty(regexp(SeriesType,'physio'))
%      SeriesType = 'undefined';
%      return;
%     %~isempty(regexp(SeriesType,'error'))
%     %~isempty(regexp(SeriesType,'qa')) ||...
%  end
  %if strcmp(sinfo.id_redcap,'NDAR_INVZKP2G8H4') && regexp(sinfo.SeriesDescription,'ABCD_T1'), keyboard; end;
  %if strcmp(sinfo.id_redcap,'NDAR_INVH2D9U5VD'), keyboard; end;
  if ~strcmp('dmri_qa',SeriesType) &&...
     ~strcmp('original',SeriesType) &&...
     ~strcmp('coil error',SeriesType) &&...
     isempty(regexp(SeriesType,'mosaic')) &&...
     isempty(regexp(SeriesType,'ge')) &&...
     isempty(regexp(SeriesType,'philips')) &&...
     isempty(regexp(SeriesType,'siemens'))

      return;
  end;
  SeriesDescription = sinfo.SeriesDescription;
  if isempty(SeriesDescription)
    SeriesType = 'undefined';
    return;
  end
  if ~isempty(regexp(SeriesDescription,'Sag_MPRAGE_T1')) ||...
     ~isempty(regexp(SeriesDescription,'ABCD_T1')) ||...
     ~isempty(regexp(SeriesDescription,'3D T1'))
    if strcmp(upper(sinfo.Manufacturer),'SIEMENS')
      fname_json = sprintf('%s/%s/%s',projinfo(sinfo.sourceID).incoming,sinfo.site,sinfo.json);
      if exist(fname_json)
        jinfo = abcd_check_json(sprintf('%s/%s/%s',projinfo(sinfo.sourceID).incoming,sinfo.site,sinfo.json));
      else
        fname_tgz = regexprep(fname_json,'json$','tgz');
        cmd = sprintf('tar xf %s ''*.json'' -O > /tmp/%s', fname_tgz, sinfo.json);
        [status, result] = unix(cmd);
        if ~status
          jinfo = abcd_check_json(sprintf('/tmp/%s',sinfo.json)); 
          cmd = sprintf('rm -f /tmp/%s', sinfo.json);
        end;
      end
      try
        for jid=1:length(jinfo.ClassifyType)
          if strcmp(jinfo.ClassifyType(jid),'ABCD-T1-NORM')
            SeriesType = 't1';
          elseif strcmp(jinfo.ClassifyType(jid),'ABCD-T1')
            SeriesType = 't1_raw';
          end 
        end
      catch ME
        fprintf('ERROR when processing %s\n', sinfo.id_redcap);
        warning(ME.message);
      end
    else
      SeriesType = 't1';
    end
  elseif ~isempty(regexp(SeriesDescription,'ABCD_T2')) ||...
     ~isempty(regexp(SeriesDescription,'3D T2'))
    if strcmp(upper(sinfo.Manufacturer),'SIEMENS')
      fname_json = sprintf('%s/%s/%s',projinfo(sinfo.sourceID).incoming,sinfo.site,sinfo.json);
      if exist(fname_json)
        jinfo = abcd_check_json(sprintf('%s/%s/%s',projinfo(sinfo.sourceID).incoming,sinfo.site,sinfo.json));
      else
        fname_tgz = regexprep(fname_json,'json$','tgz');
        cmd = sprintf('tar xf %s ''*.json'' -O > /tmp/%s', fname_tgz, sinfo.json);
        [status, result] = unix(cmd);
        if ~status
          jinfo = abcd_check_json(sprintf('/tmp/%s',sinfo.json));
          cmd = sprintf('rm -f /tmp/%s', sinfo.json);
        end;
      end 
      try
        for jid=1:length(jinfo.ClassifyType)
          if strcmp(jinfo.ClassifyType(jid),'ABCD-T2-NORM')
            SeriesType = 't2';
          elseif strcmp(jinfo.ClassifyType(jid),'ABCD-T2')
            SeriesType = 't2_raw';
          end 
        end
      catch ME
        fprintf('ERROR when processing %s\n', sinfo.id_redcap);
        warning(ME.message);
      end
    else
      SeriesType = 't2';
    end
  elseif (~isempty(regexp(SeriesDescription,'rsfMRI')) ||...
         ~isempty(regexp(SeriesDescription,'rest'))) &&...
         isempty(regexp(SeriesDescription,'Field'))
    SeriesType = 'rsfmri';
  elseif ~isempty(regexp(SeriesDescription,'rsfMRI')) &&...
         ~isempty(regexp(SeriesDescription,'Field'))
    SeriesType = 'fmri_fm';
  elseif ~isempty(regexp(SeriesDescription,'DTI\s*[12]')) ||...
         (~isempty(regexp(SeriesDescription,'ABCD_dMRI')) &&...
         isempty(regexp(SeriesDescription,'DistortionMap')))
    SeriesType = 'dmri';
  elseif ~isempty(regexp(SeriesDescription,'DTI Fieldmap P')) ||...
         ~isempty(regexp(SeriesDescription,'dMRI_DistortionMap_P'))
    SeriesType = 'dmri_fm_pa';
  elseif ~isempty(regexp(SeriesDescription,'DTI Fieldmap A')) ||...
         ~isempty(regexp(SeriesDescription,'dMRI_DistortionMap_A'))
    SeriesType = 'dmri_fm_ap';
  elseif ~isempty(regexp(SeriesDescription,'n-back'))
    SeriesType = 'nback';
  elseif ~isempty(regexp(SeriesDescription,'MID')) ||...
         ~isempty(regexp(SeriesDescription,'fMRI.*Monetary_Incentive.*'))
    SeriesType = 'mid';
  elseif ~isempty(regexp(SeriesDescription,'SST')) ||...
         ~isempty(regexp(SeriesDescription,'fMRI_task_Stop')) ||...
         ~isempty(regexp(SeriesDescription,'fMRI_Stop_Task_Run'))
    SeriesType = 'sst';
  elseif ~isempty(regexp(SeriesDescription,'fMRI')) &&...
         ~isempty(regexp(SeriesDescription,'DistortionMap_P'))
    SeriesType = 'fmri_fm_pa';
  elseif ~isempty(regexp(SeriesDescription,'fMRI')) &&...
         ~isempty(regexp(SeriesDescription,'DistortionMap_A'))
    SeriesType = 'fmri_fm_ap';
  else
    SeriesType = 'undefined';
  end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
