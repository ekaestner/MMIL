function abcd_check_missingseries(varargin)
%function abcd_check_missingseries(varargin)
%
% Optional input:
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'indir': directory to the import_repots
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'ra_infix': file suffix of ra_checklist file
%       {default = 'ra_checklist'}
%  'import_infix': file suffix of import file
%       {default = 'import_info'}
%  'outfix': file suffix of output file
%     {default = 'missingseries'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/10/17 by Feng Xue
% Last Mod: 05/19/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check lock files
fname_lck=sprintf('%s/MetaData/DAL_ABCD_QC/.missingseries.lck',getenv('HOME'));
if exist(fname_lck,'file')
  fprintf('%s\n','lock files exist!.');
  return;
end
while any([exist(sprintf('%s/MetaData/DAL_ABCD_QC/.combine_import_info.lck',getenv('HOME')),'file'),exist(sprintf('%s/batchdirs/DAL_ABCD_QC_download_ra_checklist/pbsout/.isrunning',getenv('HOME')),'file')])
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;

%Place lock file
fclose(fopen(fname_lck, 'w'));

%check input parameters
parms = check_input(varargin);

fname_import = sprintf('%s/%s_%s.csv',...
  parms.indir,parms.instem,parms.import_infix);
if ~exist(fname_import,'file')
    error('no import reports found in %s\n',indir);
end
fname_ra_checklist = sprintf('%s/%s_%s.csv',...
  parms.indir,parms.instem,parms.ra_infix);
if ~exist(fname_ra_checklist,'file')
  error('file %s not found',fname_ra_checklist);
end;
% output file
fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.outfix);
if ~exist(fname_out,'file') || parms.forceflag
  import_data = mmil_csv2struct(fname_import);
  idx = ~cellfun(@isempty, regexp({import_data.id_redcap},'^NDAR_INV.*'));
  import_data = import_data(idx);
  ra_checklist_data = mmil_csv2struct(fname_ra_checklist);

  % data reduction
  ra_checklist_data = ra_checklist_data(~cellfun('isempty',{ra_checklist_data.ra_scan_check_list_t1}));

  import_data_missing=[];
  for i=1:length(import_data)
    import_data_missing(i).pguidevent=import_data(i).pguidevent;
    import_data_missing(i).id_redcap=import_data(i).id_redcap;
    import_data_missing(i).redcap_event_name=import_data(i).redcap_event_name;
    import_data_missing(i).site=import_data(i).site;
    idx = ~cellfun('isempty',strfind({ra_checklist_data.id_redcap},import_data(i).id_redcap)); 
    if find(idx==1)
      tmp=ra_checklist_data(idx);
      idx = ~cellfun('isempty',strfind({tmp.redcap_event_name},import_data(i).redcap_event_name)); 
      validdata=find(idx==1);
      if ~isempty(validdata)
        ra_data=tmp(idx);
      end
    else
      validdata=find(idx==1);
    end
    if length(idx==2)
      if strcmp(ra_data(1).site,'OAHU')
        ra_data = ra_data(2);
      else
        ra_data = ra_data(1);
      end
      idx=[1];
    end

    %Just in case
    if (isempty(import_data(i).iqc_t1_received));import_data(i).iqc_t1_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_t1_missing = -import_data(i).iqc_t1_received;
      import_data_missing(i).iqc_t1_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_t1));ra_data(idx).ra_scan_check_list_t1=0;end;
      import_data_missing(i).iqc_t1_missing=ra_data(idx).ra_scan_check_list_t1 - import_data(i).iqc_t1_received;
    end

    %Just in case
    if (isempty(import_data(i).iqc_t2_received));import_data(i).iqc_t2_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_t2_missing = -import_data(i).iqc_t2_received;
      import_data_missing(i).iqc_t2_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_t2));ra_data(idx).ra_scan_check_list_t2=0;end;
      import_data_missing(i).iqc_t2_missing=ra_data(idx).ra_scan_check_list_t2 - import_data(i).iqc_t2_received;
    end

    %Just in case
    if (isempty(import_data(i).iqc_dmri_received));import_data(i).iqc_dmri_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_dmri_missing = -import_data(i).iqc_dmri_received;
      import_data_missing(i).iqc_dmri_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_dsc));ra_data(idx).ra_scan_check_list_dsc=0;end;
      import_data_missing(i).iqc_dmri_missing=ra_data(idx).ra_scan_check_list_dsc - import_data(i).iqc_dmri_received;
    end

    %Just in case
    if (isempty(import_data(i).iqc_rsfmri_received));import_data(i).iqc_rsfmri_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_rsfmri_missing = -import_data(i).iqc_rsfmri_received;
      import_data_missing(i).iqc_rsfmri_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_rsc));ra_data(idx).ra_scan_check_list_rsc=0;end;
      import_data_missing(i).iqc_rsfmri_missing=ra_data(idx).ra_scan_check_list_rsc - import_data(i).iqc_rsfmri_received;
    end

    %Just in case
    if (isempty(import_data(i).iqc_mid_received));import_data(i).iqc_mid_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_mid_missing = -import_data(i).iqc_mid_received;
      import_data_missing(i).iqc_mid_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_rcom));ra_data(idx).ra_scan_check_list_rcom=0;end;
      if (isempty(ra_data(idx).ra_scan_cl_mid_scan_lap));ra_data(idx).ra_scan_cl_mid_scan_lap=0;end;
      if ra_data(idx).ra_scan_cl_mid_scan_lap==2
        import_data_missing(i).iqc_mid_missing=0;
      else
        import_data_missing(i).iqc_mid_missing=ra_data(idx).ra_scan_check_list_rcom - import_data(i).iqc_mid_received;
      end
    end

    %Just in case
    if (isempty(import_data(i).iqc_sst_received));import_data(i).iqc_sst_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_sst_missing = -import_data(i).iqc_sst_received;
      import_data_missing(i).iqc_sst_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_sstrc));ra_data(idx).ra_scan_check_list_sstrc=0;end;
      if (isempty(ra_data(idx).ra_scan_cl_sst_scan_lap));ra_data(idx).ra_scan_cl_sst_scan_lap=0;end;
      if ra_data(idx).ra_scan_cl_sst_scan_lap == 2
        import_data_missing(i).iqc_sst_missing=0;
      else
        import_data_missing(i).iqc_sst_missing=ra_data(idx).ra_scan_check_list_sstrc - import_data(i).iqc_sst_received;
      end
    end

    %Just in case
    if (isempty(import_data(i).iqc_nback_received));import_data(i).iqc_nback_received=0;end;
    if isempty(validdata)
      %import_data_missing(i).iqc_nback_missing = -import_data(i).iqc_nback_received;
      import_data_missing(i).iqc_nback_missing = 0;
    else
      %Just in case
      if (isempty(ra_data(idx).ra_scan_check_list_vemorc));ra_data(idx).ra_scan_check_list_vemorc=0;end;
      if (isempty(ra_data(idx).ra_scan_cl_nbac_scan_lap));ra_data(idx).ra_scan_cl_nbac_scan_lap=0;end;
      if ra_data(idx).ra_scan_cl_nbac_scan_lap ==2
        import_data_missing(i).iqc_nback_missing=0;
      else
        import_data_missing(i).iqc_nback_missing=ra_data(idx).ra_scan_check_list_vemorc - import_data(i).iqc_nback_received;
      end
    end

    if (import_data_missing(i).iqc_t1_missing >0 ||...
       import_data_missing(i).iqc_t2_missing >0 ||...
       import_data_missing(i).iqc_dmri_missing >0 ||...
       import_data_missing(i).iqc_rsfmri_missing >0 ||...
       import_data_missing(i).iqc_mid_missing >0 ||...
       import_data_missing(i).iqc_sst_missing >0 ||...
       import_data_missing(i).iqc_nback_missing >0)

       import_data_missing(i).iqc_series_missing = 1;
    else
       import_data_missing(i).iqc_series_missing = 0;
    end

    if (import_data_missing(i).iqc_t1_missing <0 ||...
       import_data_missing(i).iqc_t2_missing <0 ||...
       import_data_missing(i).iqc_dmri_missing <0 ||...
       import_data_missing(i).iqc_rsfmri_missing <0 ||...
       import_data_missing(i).iqc_mid_missing <0 ||...
       import_data_missing(i).iqc_sst_missing <0 ||...
       import_data_missing(i).iqc_nback_missing <0)

       import_data_missing(i).iqc_ra_series_missing_memo = 1;
    else
       import_data_missing(i).iqc_ra_series_missing_memo = 0;
    end
  end
  import_data_missing = mmil_sortstruct(import_data_missing,{'site','id_redcap','redcap_event_name'});
  % write file
  mmil_struct2csv(import_data_missing,fname_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delete lock file
delete(fname_lck);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'ra_infix','ra_checklist',[],...
    'outfix','missingseries',[],...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'import_infix','import_info',[],...
    'forceflag',true,[false true],...
  });
  if parms.outdir(1) ~= '/'
     parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir);
  end
  if parms.indir(1) ~= '/', parms.indir=parms.outdir; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
