fname_studyinfo = 'data/groupavg_studyinfo.csv';

outdir = 'groupavg';
hemilist = 'lh';
%hemilist = {'rh','lh'};
mbmask_flag = 0;
%cond_contrasts = [];
%cond_contrasts = {'6','8'};
cond_contrasts = {'8','6';'8','4';'8','5'};
%cond_contrasts = 'all';

studyinfo = ts_load_studyinfo(fname_studyinfo);
tic
ts_groupavg(studyinfo,...
  'outdir',outdir,...
  'cond_contrasts',cond_contrasts,...
  'hemilist',hemilist,...
  'mbmask_flag',mbmask_flag...
);
toc


