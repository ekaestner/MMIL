%% SA NY011 Analysis Script
clear; clc;

plt_spc = 1;

%% Data Paths
subj = 'NY011_SA_SZ';

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

indir   = {'/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY11/NY11_SZ/eeg' ...
           '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY11/NY11_SA/eeg'};
 
cln_dir = {'SA' 'SZ'};

inpath_holder.(cln_dir{1}) = strsplit(ls([indir{1} '/' '*' cln_dir{1} '*.eeg']),'.eeg');
inpath_holder.(cln_dir{2}) = strsplit(ls([indir{2} '/' '*' cln_dir{2} '.eeg']),'.eeg');
inpath = [inpath_holder.(cln_dir{1}) inpath_holder.(cln_dir{2})];

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sem_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);