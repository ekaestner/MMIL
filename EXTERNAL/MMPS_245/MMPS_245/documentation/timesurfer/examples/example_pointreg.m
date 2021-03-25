hptsfile = 'eeglocs.hpts';

surffile = '/space/md5/4/data/MMILDB/RECHARGE/ARO_EEG/FSContainers/FREESURFERRECON_ARO_EEG_F104_20080321.131139_1/bem/outer_scalp4.tri';

intransfile = 'mri2head.trans';
outtransfile = 'mri2head-new.trans';

colorval = 10;
face_alpha = 0.6;
edge_alpha = 0.1;

ts_pointreg(hptsfile,surffile,'outtransfile',outtransfile,...
         'intransfile',intransfile,'colorval',colorval,...
         'face_alpha',face_alpha,'edge_alpha',edge_alpha,...
         'trans_type','mri2head');

return;

