function abcd_check_import(varargin)
%function abcd_check_import(varargin)
%
% Optional input:
%   'incoming_dir': incoming directory
%     {default = ''}
%   'unpack_dir': unpack directory
%     {default = ''}
%   'indir': metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'outdir': output directory
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'fname_info': input file containing info about
%     unique events in incoming directory
%     if empty, will use indir and infix to 
%     {default = []}
%   'infix': file suffix of input file
%     containing info about unique events in incoming directory
%     {default = 'combined_incoming_info_classified'}
%   'fname_projinfo': input file containing info about whole project
%     if empty, will use indir and infix to 
%     {default = []}
%   'outfix : file suffix for output file
%     containing info about status of unpack / recon
%     {default = 'import_info'}
%   'site': run this script on site(s) only
%     runs on all sites if empty
%     {default = []}
%   'subdir' dir to save status of import_report, such as skipped series
%            and checked pguidevents
%     {default = 'import_status'}
%    'valid_types' all valid SeriesTypes
%      {default = {'t1','t2','dmri','dmri_fm','dmri_fm_ap','dmri_fm_pa','rsfmri',...
%                  'fmri_fm','fmri_fm_ap','fmri_fm_pa','mid','sst','nback'}
%   'forceflag': ignore checkpoints and run on all series
%     {default = 1}
%
% Created:  02/17/17 by Don Hagler
% Last Mod: 07/27/17 by Feng Xue
%

%% TODO: include results for series that have already been checked in output
%%       Might be danger? what if there is correction of raw data?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
parms = check_input(varargin)

% get and write import info
if isempty(parms.site)
  fname_lck=sprintf('%s/.import_info_all.lck',parms.indir);

  parms.fname_import = sprintf('%s/%s_%s_all.csv',...
    parms.outdir,parms.outstem,parms.outfix);
  parms.fname_breakpoint = sprintf('%s/%s/%s.mat',...
    parms.outdir,parms.subdir,'pguid_checked_all');
  parms.fname_skippedseries = sprintf('%s/%s/%s.txt',...
    parms.outdir,parms.subdir,'skippedseries_all');
else
  fname_lck=sprintf('%s/.import_info_%s.lck',parms.indir,parms.site);

  parms.fname_import = sprintf('%s/%s_%s_%s.csv',...
    parms.outdir,parms.outstem,parms.outfix,parms.site);
  parms.fname_breakpoint = sprintf('%s/%s/%s_%s.mat',...
    parms.outdir,parms.subdir,'pguid_checked',parms.site);
  parms.fname_skippedseries = sprintf('%s/%s/%s_%s.txt',...
    parms.outdir,parms.subdir,'skippedseries',parms.site);
end;
%check lock files
if exist(fname_lck,'file')
  fprintf('lock files exist!.');
  return;
end
fclose(fopen(fname_lck, 'w'));

if exist(sprintf('%s/.check_incoming.lck',parms.indir),'file'), fprintf('Waiting for previous process to finish.');end;
while exist(sprintf('%s/.check_incoming.lck',parms.indir),'file')
  fprintf('.');
  pause(30);
end;
fprintf('\n');

% load incoming csv with unique pguids
incoming_info = abcd_load_csv(parms.fname_info);

idx = cellfun(@(x)ismember(x,parms.valid_types), {incoming_info.SeriesType});
incoming_info=incoming_info(idx);

projinfo = abcd_load_projinfo_all(parms);


import_info = get_import_info(incoming_info,projinfo,parms);
% write output file
if length(import_info)>0, mmil_struct2csv(import_info,parms.fname_import); end

%delete lock file
delete(fname_lck);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'incoming_dir',[],[],...
    'unpack_dir',[],[],...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'outdir','MetaData/DAL_ABCD_QC/import_reports',[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'fname_info',[],[],...
    'infix','combined_incoming_info_classified',[],...
    'fname_projinfo',[],[],...
    'outfix','import_info',[],...
    'site',[],[],...
    ...
    'subdir','import_status',[],...
    'forceflag',true,[false true],...
    ...
    'valid_types',{'t1','t2','dmri','dmri_fm','dmri_fm_ap','dmri_fm_pa','rsfmri',...
                   'fmri_fm','fmri_fm_ap','fmri_fm_pa','mid','sst','nback'},[],...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;
  if parms.indir(1) ~= '/', parms.indir= sprintf('%s/%s',getenv('HOME'),parms.indir); end;
  if ~exist(parms.outdir), mmil_mkdir(parms.outdir); end;

  if isempty(parms.fname_projinfo), parms.fname_projinfo = sprintf('%s/ProjInfo/MMIL_ProjInfo_all.csv',getenv('HOME')); end;
  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;

  if isempty(parms.fname_info)
    parms.fname_info = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.infix);
  end;
  if ~exist(parms.fname_info,'file'), error('info file %s not found',parms.fname_info); end;
  mmil_mkdir([parms.outdir '/' parms.subdir]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function import_info = get_import_info(incoming_info,projinfo,parms)
  import_info = [];
  import_info_pos = 0;

%  pguid_checked = [];
%  if parms.forceflag==0 && exist(parms.fname_breakpoint,'file')
%    isChecked=1;
%    tmp = load(parms.fname_breakpoint);
%    pguid_checked = tmp.pguidevent;
%  else
%    isChecked=0;
%  end
  %skippedseries=fopen(parms.fname_skippedseries,'a');

  if ~isempty(parms.site)
    % only for designated site.
    idx = ~cellfun('isempty',strfind({incoming_info.site},parms.site));
    incoming_info=incoming_info(idx);
  end


  pguidevent_uniq=unique({incoming_info.pguidevent});

  current_pos=0;
  totallen=length(incoming_info);

  % loop over unique events
  for i=1:length(pguidevent_uniq)
    % initialize import info
    pguidevent = pguidevent_uniq{i};

%    % Start checking json files if pguidevent equals pguid_checked
%    if isChecked==1
%      if strcmp(pguid_checked,pguidevent)
%	      isChecked = 0;
%      end
%      fprintf('%s is already checked, skipping... \n',pguidevent)
%      continue;
%    end

    idx = ~cellfun('isempty',strfind({incoming_info.pguidevent},pguidevent));
    if isempty(find(idx==1))
      continue;
    end

    % initialize import info
    import_info_pos=import_info_pos+1;

    import_info = init_info(import_info,import_info_pos);

    series_info=incoming_info(idx);
    import_info(length(import_info)).pguidevent=pguidevent;
    import_info(length(import_info)).site=series_info(1).site;
    import_info(length(import_info)).id_redcap=series_info(1).id_redcap;
    import_info(length(import_info)).redcap_event_name=series_info(1).redcap_event_name;
    site = series_info(1).site;

    % loop over series
    for j=1:length(series_info)

      fprintf('%s: checking series for %s...\n',mfilename,site);
      incoming_dir = sprintf('%s/%s',projinfo(series_info(j).sourceID).incoming,site);
      unpack_dir = sprintf('%s/%s',projinfo(series_info(j).sourceID).unpack,site);

      [fdir,fstem_json,fext] = fileparts(series_info(j).json); 

      current_pos=current_pos+1;
      fprintf('\n%d of %d finished \n', current_pos, totallen)
      fprintf('\nProcessing %s\n',fstem_json)

      % get SeriesType from SeriesDescription
      %SeriesType = get_SeriesType(series_info(j));
      SeriesType = series_info(j).SeriesType;

      tag = sprintf('iqc_%s_received',SeriesType);
      import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;

      % check for .unpacked file
      fname_check = sprintf('%s/%s/%s.unpacked',unpack_dir,fstem_json,fstem_json);
      if exist(fname_check,'file')
        tag = sprintf('iqc_%s_unpacked',SeriesType);
        import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
      end;
      % check for .unpack_error
      fname_check = sprintf('%s/%s/%s.unpack_error',unpack_dir,fstem_json,fstem_json);
      if exist(fname_check,'file')
        tag = sprintf('iqc_%s_unpack_err',SeriesType);
        import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
        import_info(length(import_info)).iqc_unpack_err = import_info(length(import_info)).iqc_unpack_err + 1;
      end;

      % check recon_required for multiband data
      %if series_info.multiband_flag
      tag = sprintf('iqc_%s_recon_req',SeriesType);
      if isfield(import_info(length(import_info)),tag)
        %First, look for json in unpack_dir, then look for it in /tmp, otherwise incoming_dir.
        fname_json = sprintf('%s/%s.json',unpack_dir,fstem_json);
        if ~exist(fname_json,'file'), fname_json = sprintf('/dev/shm/jq_cache/jq_tmp_hastgz_%s_%s/%s.json',projinfo(series_info(j).sourceID).account,site,fstem_json); end
        if ~exist(fname_json,'file'), fname_json = sprintf('/dev/shm/jq_cache/jq_tmp_missingtgz_%s_%s/%s.json',projinfo(series_info(j).sourceID).account,site,fstem_json); end
        if ~exist(fname_json,'file'), fname_json = sprintf('%s/%s.json',incoming_dir,fstem_json); end
        if exist(fname_json,'file')
          % load json file
          [jinfo,fstem_json,fname_tar,fname_pfile_json] = abcd_check_json(fname_json,incoming_dir);
  
          if ~isempty(jinfo)
            if jinfo.recon_required_flag
              tag = sprintf('iqc_%s_recon_req',SeriesType);
              import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
            end;
            if jinfo.missing_kspace_flag
              tag = sprintf('iqc_%s_missing_ks',SeriesType);
              import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
              import_info(length(import_info)).iqc_missing_ks = import_info(length(import_info)).iqc_missing_ks + 1;
            end;
    
            % check for .reconned
            fname_check = sprintf('%s/%s/%s.reconned',unpack_dir,fstem_json,fstem_json);
            if exist(fname_check,'file')
              tag = sprintf('iqc_%s_reconned',SeriesType);
              import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
            end;
    
            % check for .recon_error
            fname_check = sprintf('%s/%s/%s.recon_error',unpack_dir,fstem_json,fstem_json);
            if exist(fname_check,'file')
              tag = sprintf('iqc_%s_recon_err',SeriesType);
              import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
              import_info(length(import_info)).iqc_recon_err = import_info(length(import_info)).iqc_recon_err + 1;
            end;
          end;
        end;
      end;

      % check for .renamed
      fname_check = sprintf('%s/%s/%s.renamed',unpack_dir,fstem_json,fstem_json);
      if exist(fname_check,'file')
        tag = sprintf('iqc_%s_imported',SeriesType);
        import_info(length(import_info)).(tag) = import_info(length(import_info)).(tag) + 1;
      end;

    end;

%    % save last checked pguid
%    save(parms.fname_breakpoint,'pguidevent');

  end;
%  fclose(skippedseries);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function import_info = init_info(import_info,import_info_pos)
  import_info(import_info_pos).pguidevent = '';
  import_info(import_info_pos).site = '';
  import_info(import_info_pos).id_redcap = '';
  import_info(import_info_pos).redcap_event_name = '';
  % T1
  import_info(import_info_pos).iqc_t1_received = 0;
  import_info(import_info_pos).iqc_t1_unpacked = 0;
  import_info(import_info_pos).iqc_t1_unpack_err = 0;
  import_info(import_info_pos).iqc_t1_imported = 0;
  % T2
  import_info(import_info_pos).iqc_t2_received = 0;
  import_info(import_info_pos).iqc_t2_unpacked = 0;
  import_info(import_info_pos).iqc_t2_unpack_err = 0;
  import_info(import_info_pos).iqc_t2_imported = 0;
  % dmri
  import_info(import_info_pos).iqc_dmri_received = 0;
  import_info(import_info_pos).iqc_dmri_unpacked = 0;
  import_info(import_info_pos).iqc_dmri_unpack_err = 0;
  import_info(import_info_pos).iqc_dmri_recon_req = 0;
  import_info(import_info_pos).iqc_dmri_missing_ks = 0;
  import_info(import_info_pos).iqc_dmri_reconned = 0;
  import_info(import_info_pos).iqc_dmri_recon_err = 0;
  import_info(import_info_pos).iqc_dmri_imported = 0;
  % dmri_FM
  import_info(import_info_pos).iqc_dmri_fm_received = 0;
  import_info(import_info_pos).iqc_dmri_fm_unpacked = 0;
  import_info(import_info_pos).iqc_dmri_fm_unpack_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_recon_req = 0;
  import_info(import_info_pos).iqc_dmri_fm_missing_ks = 0;
  import_info(import_info_pos).iqc_dmri_fm_reconned = 0;
  import_info(import_info_pos).iqc_dmri_fm_recon_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_imported = 0;
  % dmri_fm_ap
  import_info(import_info_pos).iqc_dmri_fm_ap_received = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_unpacked = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_unpack_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_recon_req = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_missing_ks = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_reconned = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_recon_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_ap_imported = 0;
  % dmri_fm_pa
  import_info(import_info_pos).iqc_dmri_fm_pa_received = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_unpacked = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_unpack_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_recon_req = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_missing_ks = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_reconned = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_recon_err = 0;
  import_info(import_info_pos).iqc_dmri_fm_pa_imported = 0;
  % rsfmri
  import_info(import_info_pos).iqc_rsfmri_received = 0;
  import_info(import_info_pos).iqc_rsfmri_unpacked = 0;
  import_info(import_info_pos).iqc_rsfmri_unpack_err = 0;
  import_info(import_info_pos).iqc_rsfmri_recon_req = 0;
  import_info(import_info_pos).iqc_rsfmri_missing_ks = 0;
  import_info(import_info_pos).iqc_rsfmri_reconned = 0;
  import_info(import_info_pos).iqc_rsfmri_recon_err = 0;
  import_info(import_info_pos).iqc_rsfmri_imported = 0;
  % fmri_FM
  import_info(import_info_pos).iqc_fmri_fm_received = 0;
  import_info(import_info_pos).iqc_fmri_fm_unpacked = 0;
  import_info(import_info_pos).iqc_fmri_fm_unpack_err = 0;
  import_info(import_info_pos).iqc_fmri_fm_imported = 0;
  % fmri_fm_ap
  import_info(import_info_pos).iqc_fmri_fm_ap_received = 0;
  import_info(import_info_pos).iqc_fmri_fm_ap_unpacked = 0;
  import_info(import_info_pos).iqc_fmri_fm_ap_unpack_err = 0;
  import_info(import_info_pos).iqc_fmri_fm_ap_imported = 0;
  % fmri_fm_pa
  import_info(import_info_pos).iqc_fmri_fm_pa_received = 0;
  import_info(import_info_pos).iqc_fmri_fm_pa_unpacked = 0;
  import_info(import_info_pos).iqc_fmri_fm_pa_unpack_err = 0;
  import_info(import_info_pos).iqc_fmri_fm_pa_imported = 0;
  % mid
  import_info(import_info_pos).iqc_mid_received = 0;
  import_info(import_info_pos).iqc_mid_unpacked = 0;
  import_info(import_info_pos).iqc_mid_unpack_err = 0;
  import_info(import_info_pos).iqc_mid_recon_req = 0;
  import_info(import_info_pos).iqc_mid_missing_ks = 0;
  import_info(import_info_pos).iqc_mid_reconned = 0;
  import_info(import_info_pos).iqc_mid_recon_err = 0;
  import_info(import_info_pos).iqc_mid_imported = 0;
  % sst
  import_info(import_info_pos).iqc_sst_received = 0;
  import_info(import_info_pos).iqc_sst_unpacked = 0;
  import_info(import_info_pos).iqc_sst_unpack_err = 0;
  import_info(import_info_pos).iqc_sst_recon_req = 0;
  import_info(import_info_pos).iqc_sst_missing_ks = 0;
  import_info(import_info_pos).iqc_sst_reconned = 0;
  import_info(import_info_pos).iqc_sst_recon_err = 0;
  import_info(import_info_pos).iqc_sst_imported = 0;
  % nback
  import_info(import_info_pos).iqc_nback_received = 0;
  import_info(import_info_pos).iqc_nback_unpacked = 0;
  import_info(import_info_pos).iqc_nback_unpack_err = 0;
  import_info(import_info_pos).iqc_nback_recon_req = 0;
  import_info(import_info_pos).iqc_nback_missing_ks = 0;
  import_info(import_info_pos).iqc_nback_reconned = 0;
  import_info(import_info_pos).iqc_nback_recon_err = 0;
  import_info(import_info_pos).iqc_nback_imported = 0;

  %stats
  import_info(import_info_pos).iqc_unpack_err = 0;
  import_info(import_info_pos).iqc_recon_err = 0;
  import_info(import_info_pos).iqc_missing_ks = 0;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
