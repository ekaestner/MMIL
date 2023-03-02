function abcd_combine_import_report(varargin)
%function abcd_combine_import_report()
%
% Optional input:
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'import1_dir': directory to the first import report
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%  'import2_dir': directory to the second import report
%     {default = '/home/mmilrec18/MetaData/DAL_ABCD_QC'}
%  'import_infix': file suffix of import file
%       {default = 'import_info_all'}
%  'output_infix': file suffix of output file
%     {default = 'combined_import_info'}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/10/17 by Feng Xue
% Last Mod: 03/15/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check lock files
fname_lck='/home/mmilrec14/MetaData/DAL_ABCD_QC/.combine_import_info.lck';
if exist(fname_lck,'file')
  fprintf('%s\n','lock files exist!.');
  return;
end
while ~isempty(dir('/home/mmilrec14/MetaData/DAL_ABCD_QC/.import_info*.lck')) ||...
      ~isempty(dir('/home/mmilrec18/MetaData/DAL_ABCD_QC/.import_info*.lck'))
  fprintf('%s\n','lock files exist, waiting for previous process to finish.');
  pause(30);
end;

%Place lock file
fclose(fopen(fname_lck, 'w'));

%check input parameters
parms = check_input(varargin);
import1_info_fname=sprintf('%s/%s_%s.csv',parms.import1_dir,parms.instem,parms.import_infix);
import2_info_fname=sprintf('%s/%s_%s.csv',parms.import2_dir,parms.instem,parms.import_infix);
combined_info_fname=sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,parms.output_infix);

import1_info = mmil_csv2struct(import1_info_fname);
import2_info = mmil_csv2struct(import2_info_fname);

combined_info=combine_info(import1_info,import2_info);
% write output file
mmil_struct2csv(combined_info,combined_info_fname);

%delete lock file
delete(fname_lck);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combined_info=combine_info(import1_info,import2_info)
  idx3=[];
  for i=1:length(import2_info)
      pguidevent = import2_info(i).pguidevent;
      idx = ~cellfun('isempty',strfind({import1_info.pguidevent},pguidevent));
      if isempty(find(idx==1))
        idx3=[idx3 i];
        continue;
      end
      
      import1_info(find(idx==1)).iqc_t1_received=import1_info(find(idx==1)).iqc_t1_received+import2_info(i).iqc_t1_received;
      import1_info(find(idx==1)).iqc_t2_received=import1_info(find(idx==1)).iqc_t2_received+import2_info(i).iqc_t2_received;
      import1_info(find(idx==1)).iqc_dmri_received=import1_info(find(idx==1)).iqc_dmri_received+import2_info(i).iqc_dmri_received;
      import1_info(find(idx==1)).iqc_dmri_fm_received=import1_info(find(idx==1)).iqc_dmri_fm_received+import2_info(i).iqc_dmri_fm_received;
      import1_info(find(idx==1)).iqc_dmri_fm_ap_received=import1_info(find(idx==1)).iqc_dmri_fm_ap_received+import2_info(i).iqc_dmri_fm_ap_received;
      import1_info(find(idx==1)).iqc_dmri_fm_pa_received=import1_info(find(idx==1)).iqc_dmri_fm_pa_received+import2_info(i).iqc_dmri_fm_pa_received;
      import1_info(find(idx==1)).iqc_rsfmri_received=import1_info(find(idx==1)).iqc_rsfmri_received+import2_info(i).iqc_rsfmri_received;
      import1_info(find(idx==1)).iqc_fmri_fm_received=import1_info(find(idx==1)).iqc_fmri_fm_received+import2_info(i).iqc_fmri_fm_received;
      import1_info(find(idx==1)).iqc_fmri_fm_ap_received=import1_info(find(idx==1)).iqc_fmri_fm_ap_received+import2_info(i).iqc_fmri_fm_ap_received;
      import1_info(find(idx==1)).iqc_fmri_fm_pa_received=import1_info(find(idx==1)).iqc_fmri_fm_pa_received+import2_info(i).iqc_fmri_fm_pa_received;
      import1_info(find(idx==1)).iqc_sst_received=import1_info(find(idx==1)).iqc_sst_received+import2_info(i).iqc_sst_received;
      import1_info(find(idx==1)).iqc_nback_received=import1_info(find(idx==1)).iqc_nback_received+import2_info(i).iqc_nback_received;
      import1_info(find(idx==1)).iqc_mid_received=import1_info(find(idx==1)).iqc_mid_received+import2_info(i).iqc_mid_received;
  end
  combined_info=[import1_info;import2_info(idx3)];

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'output_infix','combined_import_info',[],...
    'import1_dir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'import2_dir','/home/mmilrec18/MetaData/DAL_ABCD_QC',[],...
    'import_infix','import_info_all',[],...
    'forceflag',true,[false true],...
  });
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
