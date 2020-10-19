function abcd_check_incomingtopc(varargin)
%function abcd_check_incomingtopc(varargin)
%
% Optional input:
%  'indir': input metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'pc_infix': file suffix of pcinfo file
%       {default = 'combined_pcinfo'}
%  'import_infix': file suffix of import file
%       {default = 'import_info'}
%  'outfix': file suffix of output file
%     {default = 'incomingtopc_info'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/01/17 by Feng Xue
% Prev Mod: 03/08/17 by Don Hagler
% Last Mod: 05/19/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check input parameters
parms = check_input(varargin);
%check lock files
fname_lck=sprintf('%s/.incomingtopc.lck',parms.outdir);
if exist(fname_lck,'file')
  fprintf('%s\n','lock files exist!.');
  return;
end
while any([exist(sprintf('%s/.merged_pcqcinfo.lck',parms.outdir),'file'),exist(sprintf('%s/.combine_import_info.lck',parms.outdir),'file')])
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;

%Place lock file
fclose(fopen(fname_lck, 'w'));


% combine import reports from all sites
fname_import = sprintf('%s/%s_%s.csv',...
  parms.indir,parms.instem,parms.import_infix);
fname_pcinfo = sprintf('%s/%s_%s.csv',...
  parms.outdir,parms.outstem,parms.pc_infix);
if ~exist(fname_import,'file')
  error('file %s not found',fname_import);
end;
if ~exist(fname_pcinfo,'file')
  error('file %s not found',fname_pcinfo);
end;

% output file
fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.outfix);
if ~exist(fname_out,'file') || parms.forceflag
  import_data = mmil_csv2struct(fname_import);
  pc_info = mmil_csv2struct(fname_pcinfo);

  import_data_missing=[];
  for i=1:length(import_data)
    import_data_missing(i).id_redcap=import_data(i).id_redcap;
    import_data_missing(i).redcap_event_name=import_data(i).redcap_event_name;
    import_data_missing(i).site=import_data(i).site;
    if ((import_data(i).iqc_t1_received==1 && import_data(i).iqc_t1_imported==0 && import_data(i).iqc_t1_unpack_err==0) || ...
        (import_data(i).iqc_t2_received==1 && import_data(i).iqc_t2_imported==0 && import_data(i).iqc_t2_unpack_err==0) || ...
        (import_data(i).iqc_dmri_received==1 && import_data(i).iqc_dmri_imported==0 && import_data(i).iqc_dmri_unpack_err==0 && import_data(i).iqc_dmri_recon_err==0) || ...
        (import_data(i).iqc_dmri_fm_received==1 && import_data(i).iqc_dmri_fm_imported==0 && import_data(i).iqc_dmri_fm_unpack_err==0 && import_data(i).iqc_dmri_fm_recon_err==0) || ...
        (import_data(i).iqc_dmri_fm_ap_received==1 && import_data(i).iqc_dmri_fm_ap_imported==0 && import_data(i).iqc_dmri_fm_ap_unpack_err==0 && import_data(i).iqc_dmri_fm_ap_recon_err==0) || ...
        (import_data(i).iqc_dmri_fm_pa_received==1 && import_data(i).iqc_dmri_fm_pa_imported==0 && import_data(i).iqc_dmri_fm_pa_unpack_err==0 && import_data(i).iqc_dmri_fm_pa_recon_err==0) || ...
        (import_data(i).iqc_rsfmri_received==1 && import_data(i).iqc_rsfmri_imported==0 && import_data(i).iqc_rsfmri_unpack_err==0 && import_data(i).iqc_rsfmri_recon_err==0) || ...
        (import_data(i).iqc_fmri_fm_received==1 && import_data(i).iqc_fmri_fm_imported==0 && import_data(i).iqc_fmri_fm_unpack_err==0) || ...
        (import_data(i).iqc_fmri_fm_ap_received==1 && import_data(i).iqc_fmri_fm_ap_imported==0 && import_data(i).iqc_fmri_fm_ap_unpack_err==0) || ...
        (import_data(i).iqc_fmri_fm_pa_received==1 && import_data(i).iqc_fmri_fm_pa_imported==0 && import_data(i).iqc_fmri_fm_pa_unpack_err==0) || ...
        (import_data(i).iqc_mid_received==1 && import_data(i).iqc_mid_imported==0 && import_data(i).iqc_mid_unpack_err==0 && import_data(i).iqc_mid_recon_err==0) || ...
        (import_data(i).iqc_sst_received==1 && import_data(i).iqc_sst_imported==0 && import_data(i).iqc_sst_unpack_err==0 && import_data(i).iqc_sst_recon_err==0) || ...
        (import_data(i).iqc_nback_received==1 && import_data(i).iqc_nback_imported==0 && import_data(i).iqc_nback_unpack_err==0 && import_data(i).iqc_nback_recon_err==0));
      import_data_missing(i).iqc_incomingtopc_event = 0;
    else
      idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
      if isempty(find(idx==1))
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      pc_info_tmp=pc_info(idx);
      if (import_data(i).iqc_t1_imported==1) && ~hasSeries(pc_info_tmp,'T1')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_t2_imported==1) && ~hasSeries(pc_info_tmp,'T2')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_dmri_imported==1) && ~hasSeries(pc_info_tmp,'dMRI')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_dmri_fm_imported==1) && ~hasSeries(pc_info_tmp,'dMRI_FM')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_dmri_fm_ap_imported==1) && ~hasSeries(pc_info_tmp,'dMRI_FM_AP')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_dmri_fm_pa_imported==1) && ~hasSeries(pc_info_tmp,'dMRI_FM_PA')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_rsfmri_imported==1) && ~hasSeries(pc_info_tmp,'rsfMRI')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_fmri_fm_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_FM')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_fmri_fm_ap_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_FM_AP')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_fmri_fm_pa_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_FM_PA')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_mid_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_MID_task')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_sst_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_SST_task')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      if (import_data(i).iqc_nback_imported==1) && ~hasSeries(pc_info_tmp,'fMRI_nBack_task')
        import_data_missing(i).iqc_incomingtopc_event = 0;
        continue
      end
      import_data_missing(i).iqc_incomingtopc_event = 1;
    end
  end
  % write file
  import_data_missing=set_pguidevents(import_data_missing);
  import_data_missing = mmil_sortstruct(import_data_missing,{'site','id_redcap','redcap_event_name'});
  mmil_struct2csv(import_data_missing,fname_out);
end

%delete lock file
delete(fname_lck);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'pc_infix','combined_pcinfo',[],...
    'outfix','incomingtopc_info',[],...
    'import_infix','import_info',[],...
    'forceflag',true,[false true],...
  });
  if parms.outdir(1) ~= '/'
     parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir);
  end
  if parms.indir(1) ~= '/', parms.indir=parms.outdir; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hasSeries = hasSeries(pc_info,sType_source)
  idx = ~cellfun('isempty',strfind({pc_info.SeriesType},sType_source)); 
  if isempty(find(idx==1))
    hasSeries=0;
  else
    hasSeries=1;
  end
return;

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
