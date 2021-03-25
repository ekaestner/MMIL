function SL_lng_eve(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial events work on %s \n'],sbj)

infile = [fcfg.dat_fld sbj '_overall_data.mat'];
outpath = fcfg.dat_fld;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Events
cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 3 4};
cfg.new_events  = {601 603 604 };
cfg.crt_alt_eve = 'lng_tot_nse';
cfg.use_alt_eve = 'trialinfo';
bcc_dat = ft_func(@ft_redefine_events,cfg,bcc_dat);

%% Size
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end