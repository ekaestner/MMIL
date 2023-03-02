function abcd_incomingtopc_report(varargin)
%function abcd_incomingtopc_report(varargin)
%
% Optional input:
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'import_dir': directory to the import_repots
%     {default = 'import_reports'}
%  'pc_infix': file suffix of pcinfo file
%       {default = 'pcinfo'}
%  'import_infix': file suffix of import file
%       {default = 'import_info'}
%  'output_infix': file suffix of output file
%     {default = 'incomingtopc_info'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  12/05/16 by Jose Teruel
% Prev Mod: 03/08/17 by Don Hagler
% Last Mod: 03/09/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input parameters
parms = check_input(varargin);

% combine import reports from all sites
fname_import = sprintf('%s/%s_%s.csv',...
  parms.outdir,parms.outstem,parms.import_infix);
if ~exist(fname_import,'file') || parms.forceflag
  indir = sprintf('%s/%s',parms.outdir,parms.import_dir);
  if ~exist(indir,'dir')
    error('import reports dir %s not found',indir);
  end;
  flist = dir(sprintf('%s/%s_%s_*.csv',...
    indir,parms.outstem,parms.import_infix));
  if isempty(flist)
    error('no import reports found in %s\n',indir);
  end;
  import_info = [];
  for i=1:length(flist)
    fname_tmp = sprintf('%s/%s',...
      indir,flist(i).name);
    tmp_info = mmil_readtext(fname_tmp);
    if isempty(import_info)
      import_info = tmp_info;
    else
      import_info = cat(1,import_info,tmp_info(2:end,:));
    end;
  end;
  mmil_write_csv(fname_import,import_info);
end;
fname_pcinfo = sprintf('%s/%s_%s.csv',...
  parms.outdir,parms.outstem,parms.pc_infix);
if ~exist(fname_pcinfo,'file')
  error('file %s not found',fname_pcinfo);
end;

% output file
fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.output_infix);
if ~exist(fname_out,'file') || parms.forceflag
  waiting_process = {'pguidevent' 'id_redcap', 'redcap_event_name'};
  import_data = mmil_csv2struct(fname_import);
  pc_info = mmil_csv2struct(fname_pcinfo);
  % data reduction
  import_data = import_data(~cellfun('isempty',{import_data.iqc_t1_imported}));
  import_data_missing=[];
  for i=1:length(import_data)
    import_data_missing(i).id_redcap=import_data(i).id_redcap;
    import_data_missing(i).redcap_event_name=import_data(i).redcap_event_name;
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
      if (import_data(i).iqc_t1_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'T1')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_t2_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'T2')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_dmri_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'dMRI')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_dmri_fm_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'dMRI_FM')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_dmri_fm_ap_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'dMRI_FM_AP')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_dmri_fm_pa_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'dMRI_FM_PA')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_rsfmri_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'rsfMRI')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_fmri_fm_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_FM')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_fmri_fm_ap_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_FM_AP')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_fmri_fm_pa_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_FM_PA')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_mid_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_MID_task')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_sst_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_SST_task')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      if (import_data(i).iqc_nback_imported==1)
        idx = ~cellfun('isempty',strfind({pc_info.pGUID},import_data(i).id_redcap)); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
        pc_info_tmp=pc_info(idx);
        idx = ~cellfun('isempty',strfind({pc_info_tmp.SeriesType},'fMRI_nBack_task')); 
        if isempty(find(idx==1))
          import_data_missing(i).iqc_incomingtopc_event = 0;
          continue
        end
      end
      import_data_missing(i).iqc_incomingtopc_event = 1;
    end
  end
  % write file
  mmil_struct2csv(import_data_missing,fname_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'pc_infix','pcinfo',[],...
    'output_infix','incomingtopc_info',[],...
    'import_dir','import_reports',[],...
    'import_infix','import_info',[],...
    'forceflag',true,[false true],...
  });
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
