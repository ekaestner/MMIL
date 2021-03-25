function abcd_summarize_pc_report(varargin)
%function abcd_summarize_pc_report([options])
%
% optional input:
%   'indir': input directory containing pc_info spreadsheet
%     {default = '/home/mmilrec14/MetaData/DAL_ABCD_QC'}
%   'outdir': output directory
%     {default = './'}
%   'fname': pcinfo summary csv file
%     {default = 'DAL_ABCD_QC_pcinfo.csv'}
%   'sel_field': field for compliance
%     {default = 'ABCD_Compliant'}
%   'sel_field2': field for complete scan
%     {default = 'Completed'}
%   'date_field': field with study date
%     {default = 'StudyDate'}
%   'outstem': output file stem
%     {default = 'pcinfo_nc'}
%   'outstem2': output file stem
%     {default = 'pcinfo_incomplete'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%
% Created:  12/08/16 by Jose Teruel
% Prev Mod: 01/16/17 by Don Hagler
% Last Mod: 02/17/17 by Jose Teruel
%

% NOTE: based on abcd_summarize_pc, last mod 12/08/16, created by Jose
% Teruel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input
parms = check_input(varargin);

% check output, maybe quit depending on forceflag
output_exists_flag = check_output(parms);
if output_exists_flag && ~parms.forceflag
  fprintf('%s: WARNING: not overwriting existing output\n',mfilename);
  return;
end;

mmil_mkdir(parms.outdir);

% loac csv summary
pc_summary = load_summary(parms);

% summary non compliant and incomplete (current and previous month)
[pc_summary_nc_current, pc_summary_incomplete_current,...
    pc_summary_nc_past, pc_summary_incomplete_past, p_date_f, c_date_f] = sort_summary(pc_summary, parms);

% write csv files with one row for each series
write_outputs(pc_summary_nc_current, pc_summary_incomplete_current,...
    pc_summary_nc_past, pc_summary_incomplete_past, p_date_f, c_date_f, parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'rootdir','/space/syn05/1/data/MMILDB',[],...
    'indir','/home/mmilrec14/MetaData/DAL_ABCD_QC',[],...
    'outdir','/home/mmilrec14/MetaData/DAL_ABCD_QC/compliance_reports',[],...    
    'fname', 'DAL_ABCD_QC_pcinfo.csv',[],... 
    'sel_field', 'ABCD_Compliant',[],...
    'sel_field2', 'Completed',[],...
    'date_field', 'StudyDate',[],...
    'outstem','DAL_ABCD_QC_pcinfo_notcompliant',[],...
    'outstem2','DAL_ABCD_QC_pcinfo_incomplete',[],...    
    'forceflag',false,[false true],...
  });
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output_exists_flag = check_output(parms)
  output_exists_flag = true;
  suffix_list = {'series','sessions'};
  for i=1:length(suffix_list)
    suffix = suffix_list{i};
    fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,suffix);
    if ~exist(fname_out,'file')
      output_exists_flag = false;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pc_summary = load_summary(parms)
  fname = fullfile(parms.indir,parms.fname);
  pc_summary=mmil_readtext(fname);  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pc_summary_nc_current, pc_summary_incomplete_current,...
    pc_summary_nc_past, pc_summary_incomplete_past, p_date_f,...
    c_date_f] = sort_summary(pc_summary, parms)
  %select date field and exclude anything not acquired current or
  %previous month
  [p_date_i, p_date_f, c_date_i, c_date_f] = setdatelimits();
  date_index = find(strcmp(pc_summary(1,:),parms.date_field));
  pc_index = find(strcmp(pc_summary(1,:),parms.sel_field));
  complete_index = find(strcmp(pc_summary(1,:),parms.sel_field2));
  pc_summary_nc_current = pc_summary(1,:);
  pc_summary_incomplete_current = pc_summary(1,:);
  pc_summary_nc_past = pc_summary(1,:);
  pc_summary_incomplete_past = pc_summary(1,:);    
  nc_index_c= 2; 
  ncomplete_index_c = 2;
  nc_index_p= 2; 
  ncomplete_index_p = 2;    
  for i=2:size(pc_summary,1)
    if str2double(pc_summary{i,date_index})<= str2double(c_date_f) &&...
       str2double(pc_summary{i,date_index})>= str2double(c_date_i)         
      if (strcmp(pc_summary{i,pc_index},'No'))
        pc_summary_nc_current(nc_index_c,:) = pc_summary(i,:);
        nc_index_c = nc_index_c+1;
      elseif (strcmp(pc_summary{i,pc_index},'Yes')) && (pc_summary{i,complete_index}==0)
        pc_summary_incomplete_current(ncomplete_index_c,:) = pc_summary(i,:);
        ncomplete_index_c = ncomplete_index_c+1;
      end
    elseif str2double(pc_summary{i,date_index})<= str2double(p_date_f) &&...
           str2double(pc_summary{i,date_index})>= str2double(p_date_i)
      if (strcmp(pc_summary{i,pc_index},'No'))
        pc_summary_nc_past(nc_index_p,:) = pc_summary(i,:);
        nc_index_p = nc_index_p+1;
      elseif (strcmp(pc_summary{i,pc_index},'Yes')) && (pc_summary{i,complete_index}==0)
        pc_summary_incomplete_past(ncomplete_index_p,:) = pc_summary(i,:);
        ncomplete_index_p = ncomplete_index_p+1;
      end                     
    end
  end
  p_date_f = p_date_f(1:6);
  c_date_f = c_date_f(1:6);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_date_i, p_date_f, c_date_i, c_date_f] = setdatelimits()
  c_date = datevec(date);
  month = num2str(c_date(2));
  day = num2str(c_date(3));
  if length(month)==1, month = sprintf('0%s',month); end
  if length(day)==1, day = sprintf('0%s',day); end    
  c_date_f = sprintf('%d%s%s',c_date(1),month,day);
  c_date_i = sprintf('%d%s01',c_date(1),month);
  p_year = c_date(1);
  p_month = c_date(2)-1;
  if p_month==0
    p_month=12;
    p_year=p_year-1;
  end
  month = num2str(p_month);
  if length(month)==1, month = sprintf('0%s',month); end
  p_date_f = sprintf('%d%s%d',p_year,month,31);
  p_date_i = sprintf('%d%s01',p_year,month);  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_outputs(pc_summary_nc_current, pc_summary_incomplete_current,...
    pc_summary_nc_past, pc_summary_incomplete_past, p_date_f, c_date_f, parms)

  fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,c_date_f);
  fprintf('%s: writing output to %s...\n',mfilename,fname_out);
  mmil_write_csv(fname_out,pc_summary_nc_current);
  
  fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem2,c_date_f);
  fprintf('%s: writing output to %s...\n',mfilename,fname_out);
  mmil_write_csv(fname_out,pc_summary_incomplete_current);
  
  fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,p_date_f);
  fprintf('%s: writing output to %s...\n',mfilename,fname_out);
  mmil_write_csv(fname_out,pc_summary_nc_past);
  
  fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem2,p_date_f);
  fprintf('%s: writing output to %s...\n',mfilename,fname_out);
  mmil_write_csv(fname_out,pc_summary_incomplete_past);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

