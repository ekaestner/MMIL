% Last Mod: 08/03/11 by Don Hagler

qcflag = 0;
sssflag = 0;

if sssflag
  prefix = 'proc_sss';
else
  prefix = 'proc';
end;

batchname = [prefix '_MEG'];

MMIL_Process_MEG_Exams('TEST',...
  'RAW_sssflag',sssflag,'PROC_prefix',prefix,...
  'qcflag',qcflag,...
  'batchname',batchname);

