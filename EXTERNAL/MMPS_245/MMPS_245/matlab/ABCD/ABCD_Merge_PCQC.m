function abcd_merge_pcqc(varargin)
%function abcd_merge_pcqc(varargin)
%
% Optional input:
%  'outdir': output metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%  'instem': input file stem
%     {default = 'DAL_ABCD_QC'}
%  'outstem': output file stem
%     {default = 'DAL_ABCD_QC'}
%  'pc_infix': file suffix of pcinfo file
%       {default = 'pcinfo'}
%  'qc_infix': file suffix of qcinfo file
%       {default = 'qcinfo'}
%  'outfix': file suffix of output file
%     {default = 'merged_pcqcinfo'}
%  'fname_projinfo': file name of whole project info
%     {default = []}
%  'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  03/21/17 by Feng Xue
% Last Mod: 07/30/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%check input parameters
parms = check_input(varargin);

%check lock files
fname_lck=sprintf('%s/.%s.lck',parms.outdir,parms.outfix);
if exist(fname_lck,'file')
  fprintf('%s','lock files exist!.');
  return;
else
  %Place lock file
  fclose(fopen(fname_lck, 'w'));
end

fname_pcqcinfo_combined = sprintf('%s/%s_%s.csv',...
  parms.outdir,parms.outstem,parms.outfix);

if ~exist(fname_pcqcinfo_combined,'file') || parms.forceflag

  projinfo = abcd_load_projinfo_all(parms);
 
  pc_info = [];
  qc_info = []; 
  % combine import reports from all sites
  for i=1:length(projinfo)
    fname_pcinfo = sprintf('/home/%s/MetaData/%s/%s_%s.csv',...
      projinfo(i).account,parms.instem,parms.instem,parms.pc_infix);
    fname_qcinfo = sprintf('/home/%s/MetaData/%s/%s_%s.csv',...
      projinfo(i).account,parms.instem,parms.instem,parms.qc_infix);
    if exist(fname_pcinfo,'file')
      try
        pc_info_tmp = abcd_load_csv(fname_pcinfo);
        %kick out duplicates
        [idx,idx]=unique({pc_info_tmp.VisitID});
        pc_info_tmp=pc_info_tmp(sort(idx));
        pc_info = [pc_info;pc_info_tmp];
      catch ME
        warning(ME.message);
      end
    end
    if exist(fname_qcinfo,'file')
      try
        qc_info_tmp = abcd_load_csv(fname_qcinfo);
        %kick out duplicates
        [idx,idx]=unique({qc_info_tmp.VisitID});
        qc_info_tmp=qc_info_tmp(sort(idx));
        qc_info = [qc_info;qc_info_tmp];
      catch ME
        warning(ME.message);
      end
    end
  end
  
  fname_pcinfo_combined = sprintf('%s/%s_combined_%s.csv',...
    parms.outdir,parms.outstem,parms.pc_infix);
  fname_qcinfo_combined = sprintf('%s/%s_combined_%s.csv',...
    parms.outdir,parms.outstem,parms.qc_infix);
  
  [idx,idx]=unique({pc_info.VisitID},'last');
  pc_info=pc_info(sort(idx));
  
  [idx,idx]=unique({qc_info.VisitID},'last');
  qc_info=qc_info(sort(idx));

  
  if ~isempty(pc_info)
    pc_info = mmil_sortstruct(pc_info,{'SiteName','pGUID','EventName'});
    mmil_struct2csv(pc_info,fname_pcinfo_combined); 
    fname_cache = sprintf('%s/%s/%s_combined_%s_%s.mat',...
      parms.outdir,'cache',parms.outstem,parms.pc_infix,datestr(now,'yyyymmdd'));
    abcd_info = pc_info;
    save(fname_cache,'abcd_info');
  end;
  if ~isempty(qc_info)
    qc_info = mmil_sortstruct(qc_info,{'VisitID'});
    mmil_struct2csv(qc_info,fname_qcinfo_combined);
  end;

  clear pc_info_tmp pc_info_tmp; 

  pc_info_fieldnames=fieldnames(pc_info);
  qc_info_fieldnames=fieldnames(qc_info);
  [idx1 idx2]=ismember({qc_info.VisitID},{pc_info.VisitID});
  [idx3 idx4]=ismember(qc_info_fieldnames,pc_info_fieldnames);
  qc_info_common=qc_info(idx1);
  pc_info_common=pc_info(idx2(idx2>0));
  pc_info_cell=[pc_info_fieldnames';struct2cell(pc_info_common)'];
  qc_info_cell=[qc_info_fieldnames';struct2cell(qc_info_common)'];
  pcqc_info_cell=[pc_info_cell qc_info_cell(:,find(idx3==0))];


  %Put back non-common data
  pc_info_uniq=pc_info(setdiff(1:length(pc_info),idx2)); 
  pc_info_uniq_cell=struct2cell(pc_info_uniq)';
  pcqc_info_cell(length(pcqc_info_cell)+1:length(pcqc_info_cell)+length(pc_info_uniq),1:size(pc_info_uniq_cell,2))=pc_info_uniq_cell;

  %Sort by site then by VisitID
  pcqc_info_cell=[pcqc_info_cell(1,:);sortrows(pcqc_info_cell(2:end,:),[5 2])];

  qc_info_uniq=qc_info(~idx1);
  qc_info_uniq_cell=struct2cell(qc_info_uniq)';

  %Sort by VisitID
  qc_info_uniq_cell=sortrows(qc_info_uniq_cell,1);

  qc_info_space=length(pcqc_info_cell)+1:length(pcqc_info_cell)+length(qc_info_uniq);

  pcqc_info_cell(qc_info_space,idx4(idx4>0))=qc_info_uniq_cell(:,find(idx3==1));
  pcqc_info_cell(qc_info_space,size(pc_info_uniq_cell,2)+1:end)=qc_info_uniq_cell(:,find(idx3==0));

  if ~isempty(pcqc_info_cell)
    pcqc_info = cell2struct(pcqc_info_cell(2:end,:),pcqc_info_cell(1,:),2);
    mmil_struct2csv(pcqc_info,fname_pcqcinfo_combined);
    fname_cache = sprintf('%s/%s/%s_%s_%s.mat',...
      parms.outdir,'cache',parms.outstem,parms.outfix,datestr(now,'yyyymmdd'));
    abcd_info = pcqc_info;
    save(fname_cache,'abcd_info');
  end;

end

%delete lock file
delete(fname_lck);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'instem','DAL_ABCD_QC',[],...
    'pc_infix','pcinfo',[],...
    'qc_infix','qcinfo',[],...
    'outfix','merged_pcqcinfo',[],...
    'fname_projinfo',[],[],...
    'forceflag',true,[false true],...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end

  if isempty(parms.fname_projinfo), parms.fname_projinfo = sprintf('%s/ProjInfo/MMIL_ProjInfo_all.csv',getenv('HOME')); end;
  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
