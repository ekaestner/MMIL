%% SA NY20 Analysis Script
clear; clc;

plt_spc = 0;

%% Data Paths
subj = 'NY20_SA';

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

indir = '/space/mdeh5/1/halgdev/data/iEEG_NYU/NY20/NY2006/NY2006@20050907_1050_1133MAP';
 
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
cfg.task = 'SA';
sem_dat  = ft_func(@ft_sasz_addevent,cfg,sem_dat);