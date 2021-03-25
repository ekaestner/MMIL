%% Initial wrapper script for iSA/SZ project
% This script will perform basic pre-processing & look for electrodes with
% a response to either SA, SZ, or SA/SZ

clear; clc;
addpath /home/ekaestne/fieldtrip-20131031/;
ft_defaults

indir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/raw_data/';

subj = {'NY007' ...
    'NY008' ...
    'NY011' ...
    'NY013' ...
    'NY014' ...
    'NY017' ...
    'NY018' ...
    'NY023' ...
    'NY024' ...
    };

infile = strcat(subj,'_epoch_data.mat');

outdir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/proc_data';

cfg = [];
cfg.preprocess = 1;
cfg.rej = 1;
cfg.init_plot = 1;
cfg.baseline = 1;
cfg.stats = 1;
cfg.plot = 1;

for isub = 1:length(subj)
    Macro_Init_proc(indir,subj{isub},infile{isub},outdir,cfg)
end