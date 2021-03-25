%% SZ NY17 Analysis Script - has 2 SZ
clear; clc;

plt_spc = 0;

%% Data Paths
subj = 'NY17_SZ';

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY17/NY17_SZ/eeg';
 
cln_dir = {'CLIN1'};

inpath_holder.(cln_dir{1}) = strsplit(ls([indir '/' '*.eeg']),'.eeg');
inpath = [inpath_holder.(cln_dir{1})];

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sem_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
inpath = cellfun(@(x) [x '.eeg'],inpath,'uni',0);

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
cfg.continuous = 'no';
cfg.dataset    = inpath;
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

cfg      = [];
cfg.task = 'SZ';
sem_dat  = ft_func(@ft_sasz_addevent,cfg,sem_dat);