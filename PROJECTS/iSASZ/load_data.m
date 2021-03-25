%% Read in & save Telemetry Data for the ieeg SA/SZ project
clear;clc;
addpath /home/ekaestne/fieldtrip-20131031/;
ft_defaults

outdir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/raw_data';

%% %%%%%%%%%%%%%%%%%%%%% NYU Telemetry - Both SA & SZ %%%%%%%%%%%%%%%%%%%%%
%% NY007
subj = 'NY007';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
% contained in an eeg file, 541 trials, 176 channels
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY07/NY07_SA/eeg';
sa_file  = [sa_indir '/' 'NY07_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
% contained in an eeg file, 398 trials, 176 channels
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY07/NY07_SZ/eeg';
sz_file  = [sz_indir '/' 'NY07_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY008
subj = 'NY008';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
% contained in an eeg file, 387 trials, 176 channels
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY08/NY08_SA/eeg/';
sa_file  = [sa_indir '/' 'NY08_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
% contained in an eeg file, 365 trials, 176 channels
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY08/NY08_SZ/eeg';
sz_file  = [sz_indir '/' 'NY08_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

sz_epoch_data.label = sa_epoch_data.label; % - EJK sz had some weird append to it's chan names

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY011 - EJK removed 48 channels from SZ to make it match with SA
subj = 'NY011';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY11/NY11_SA/eeg';
sa_file  = [sa_indir '/' 'NY11_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY11/NY11_SZ/eeg';
sz_file  = [sz_indir '/' 'NY11_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

chans = sz_epoch_data.label(1:48);
chans = strcat('-',chans(:));
cfg = []; cfg.channel = {'all',chans{:}};
sz_epoch_data = ft_func(@ft_preprocessing,cfg,sz_epoch_data);


sz_epoch_data.label = sa_epoch_data.label; % - EJK sz had some weird append to it's chan names

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY012 - PROBLEM with missing SA data (contact Thomas) - EJK
subj = 'NY012';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '.cn';
sa_file  = [sa_indir '/' 'NY11_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY11/NY11_SZ/eeg';
sz_file  = [sz_indir '/' 'NY11_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

keep outdir

%% NY013
subj = 'NY013';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY13/NY13_SA/eeg';
sa_file  = [sa_indir '/' 'NY13_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY13/NY13_SZ/eeg';
sz_file  = [sz_indir '/' 'NY13_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY014
subj = 'NY014';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY14/NY14_SA/eeg';
sa_file  = [sa_indir '/' 'NY14_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY14/NY14_SZ/eeg';
sz_file  = [sz_indir '/' 'NY14_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY017
subj = 'NY017';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY17/NY17_SA/eeg';
sa_file  = [sa_indir '/' 'NY17_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY17/NY17_SZ/eeg';
sz_file  = [sz_indir '/' 'NY17_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY018
subj = 'NY018';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY18/NY18_SA/eeg';
sa_file  = [sa_indir '/' 'NY18_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY18/NY18_SZ/eeg';
sz_file  = [sz_indir '/' 'NY18_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY019 - PROBLEM with missing SA data (contact Thomas) - EJK
subj = 'NY';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '';
sa_file  = [sa_indir '/' ''];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '';
sz_file  = [sz_indir '/' ''];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY023
subj = 'NY023';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY23/NY23_SA/eeg';
sa_file  = [sa_indir '/' 'NY23_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY23/NY23_SZ/eeg';
sz_file  = [sz_indir '/' 'NY23_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY024
subj = 'NY024';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY24/NY24_SA/eeg';
sa_file  = [sa_indir '/' 'NY24_SA.eeg'];

cfg = [];
cfg.datafile = sa_file;
cfg.subject = subj;
cfg.task = 'SA';

sa_epoch_data = ft_func(@ft_preprocessing,cfg);
sa_epoch_data = ft_func(@ft_sasz_addevent,sa_epoch_data);

% SZ
sz_indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY24/NY24_SZ/eeg';
sz_file  = [sz_indir '/' 'NY24_SZ.eeg'];

cfg = [];
cfg.datafile = sz_file;
cfg.subject = subj;
cfg.task = 'SZ';

sz_epoch_data = ft_func(@ft_preprocessing,cfg);
sz_epoch_data = ft_func(@ft_sasz_addevent,sz_epoch_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir

%% NY349 - EJK this is slightly different than the other NYU data as it is EDF, and must be epoched
% The trigger channel in this data leads me to believe that the data is
% useless. Perhaps check other data from this subject to see if the trigger
% channel for all subjects is worthless
subj = 'NY349';
subj_outdir = [outdir '/' subj];
if ~exist(subj_outdir,'dir'); mkdir(subj_outdir);end

% SA
sa_indir = '/space/mdeh2/6/halgdev/incoming/iEEG_NYU/NY349/NY349/NY349_SZ_08022012';
sa_file  = [sa_indir '/' 'NY349_SZ_08022012.edf'];

cfg = [];
cfg.datafile = sa_file;
cfg.continuous = 'yes';
cfg.subject = subj;
cfg.task = 'SA';

sa_cont_data = ft_func(@ft_preprocessing,cfg);

% Epoch Data
cfg = [];
cfg.trialfun = 'ft_NYU_iEEG_trialfun';
cfg.dataset = sprintf(sa_file);
cfg.timelim = [1.9 2.1];
cfg.pre = 1.5;
cfg.pos = 2.5;
trl = ft_definetrial(cfg);
cfg = []; cfg.trl = trl.trl;

sa_epoch_data = ft_func(@ft_redefinetrial,cfg,sa_cont_data);

% SZ
sz_indir = '/space/mdeh2/6/halgdev/incoming/iEEG_NYU/NY349/NY349/NY349_SZ_08022012';
sz_file  = [sz_indir '/' 'NY349_SZ_08022012.edf'];

cfg = [];
cfg.datafile = sz_file;
cfg.continuous = 'yes';
cfg.subject = subj;
cfg.task = 'SZ';

sz_cont_data = ft_func(@ft_preprocessing,cfg);

% Epoch Data
cfg = [];
cfg.trialfun = 'ft_NYU_iEEG_trialfun';
cfg.dataset = sprintf(sz_file);
cfg.minduration = 1;
cfg.pre = 1.5;
cfg.pos = 2.5;
trl = ft_definetrial(cfg);
cfg = []; cfg.trl = trl.trl;

sz_epoch_data = ft_func(@ft_redefinetrial,cfg,sz_cont_data);

cfg = [];
epoch_data = ft_func(@ft_appenddata,cfg,sa_epoch_data,sz_epoch_data);
epoch_data.cfg.task = 'SA/SZ';

save([subj_outdir '/' subj '_epoch_data.mat'],'epoch_data')

keep outdir






