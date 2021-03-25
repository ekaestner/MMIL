%% SA MG39 Analysis Script

clear; clc;

trg_chn = 1;
plt_spc = 1;

%% Data Paths
subj = 'MG039_SZ';

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

indir   = '/home/halgdev/analysis/iEEG_MGH/MG39';
 
cln_dir = {'CLIN1'};

inpath_holder.(cln_dir{1}) = strsplit(ls([indir '/' '*.edf']),'.edf');
inpath = [inpath_holder.(cln_dir{1})];

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sem_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
inpath = cellfun(@(x) [x '.edf'],inpath,'uni',0);

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
cfg.continuous = 'yes';
cfg.dataset    = inpath;
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

if trg_chn
    dat = 'edf';
    cfg = [];
    cfg.channel = {'Pulse*'}; % - EJK might need to revisit depending on older NYU data
    sem_dat = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1}));  
    save(['/space/mdeh4/1/halgdev/projects/mmilanguage/TriggerChannel' '/' subj '_' dat '.mat'],'sem_dat');
else
    cfg      = [];
    cfg.task = 'SA';
    sem_dat  = ft_func(@ft_sasz_addevent,cfg,sem_dat);
end