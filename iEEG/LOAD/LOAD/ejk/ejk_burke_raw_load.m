function [out_dat,trg_chn,fcfg] = ejk_burke_raw_load(fcfg,hld_dat)

fcfg.in_dir  = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir');
fcfg.cln_fld = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_fld');
fcfg.cln_dir = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir');
fcfg.end_dir = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir');

%% Burke Load
% BURKE DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idr_hld     = fcfg.in_dir{1};
cln_fld_hld = fcfg.cln_fld{1};
cln_dir_hld = fcfg.cln_dir{1};

tot_fle = dir([idr_hld '/' cln_fld_hld '/' '*' cln_dir_hld '*' fcfg.end_dir{1}]); tot_fle = {tot_fle(:).name};
hld_dta = load([idr_hld '/' cln_fld_hld '/' tot_fle{1}]);

%% IF MICRO DATA
% MICRO DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itn_dir = '/space/seh9/1/halgdev/projects/micromacro/3_25_19/';

itn_fle_nme = { 'SD018_Mon_190325_084425.rhd','SD018_Mon_190325_084525.rhd','SD018_Mon_190325_084625.rhd','SD018_Mon_190325_084725.rhd','SD018_Mon_190325_084825.rhd','SD018_Mon_190325_084925.rhd','SD018_Mon_190325_085025.rhd' ...
                'SD018_Mon_190325_085125.rhd','SD018_Mon_190325_085225.rhd','SD018_Mon_190325_085326.rhd','SD018_Mon_190325_085426.rhd','SD018_Mon_190325_085526.rhd','SD018_Mon_190325_085626.rhd','SD018_Mon_190325_085726.rhd','SD018_Mon_190325_085826.rhd','SD018_Mon_190325_085926.rhd','SD018_Mon_190325_090026.rhd' ...
                'SD018_Mon_190325_090126.rhd','SD018_Mon_190325_090226.rhd','SD018_Mon_190325_090326.rhd','SD018_Mon_190325_090426.rhd','SD018_Mon_190325_090526.rhd','SD018_Mon_190325_090626.rhd','SD018_Mon_190325_090726.rhd','SD018_Mon_190325_090826.rhd','SD018_Mon_190325_090926.rhd','SD018_Mon_190325_091026.rhd' ...
                'SD018_Mon_190325_091126.rhd','SD018_Mon_190325_091226.rhd','SD018_Mon_190325_091326.rhd','SD018_Mon_190325_091426.rhd','SD018_Mon_190325_091526.rhd','SD018_Mon_190325_091626.rhd','SD018_Mon_190325_091726.rhd','SD018_Mon_190325_091826.rhd' };
   
dta_itn = [];
fss_itn = [];
for file = 1:numel(itn_fle_nme)
    fprintf('%s %i complete\n',itn_fle_nme{file},file)
    dta_itn_tmp = read_Intan_RHD2000_file_nogui([itn_dir '/' itn_fle_nme{file}]);
    dta_itn = [dta_itn dta_itn_tmp.amplifier_data];
    fss_itn = [fss_itn dta_itn_tmp.frequency_parameters.amplifier_sample_rate];
end

clear dta_itn_tmp

dta_itn_dsm = resample(dta_itn', hld_dta.Fs, fss_itn(1))';

clear dta_itn

% ORIGINAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idr_hld     = fcfg.in_dir{2};
try cln_fld_hld = fcfg.cln_fld{2}; catch; cln_fld_hld = fcfg.cln_fld{1}; end
cln_dir_hld = fcfg.cln_dir{2};

tot_fle = dir([idr_hld '/' cln_fld_hld '/' '*' cln_dir_hld '*' fcfg.end_dir{2}]); tot_fle = {tot_fle(:).name};

targetSignals = {'TRIG'}; %{'Macro1', 'Macro2', 'Macro7', 'TRIG' };
[header, dta_cln] = edfread([idr_hld '/' cln_fld_hld '/' tot_fle{1}], 'targetSignals', targetSignals);

dta_cln = resample(dta_cln', hld_dta.Fs, header.frequency(1))';

% Trigger Channels
trg_chn = dta_cln(1,:);

% Align data
dta_ind_hld = [1710126 2760685];


trg_chn = trg_chn(1,dta_ind_hld(1):dta_ind_hld(2));
trg_chn = trg_chn * -1;
trg_chn(trg_chn<0) = 0;

% Append Micro & Macro Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hld_dta.dat = hld_dta.dat(dta_ind_hld(1):dta_ind_hld(2),:)';

% EC
chn_num = [1:13];

hld_dta.dat = [ hld_dta.dat ; dta_itn_dsm(chn_num,:) ];
hld_dta.lab = [ hld_dta.lab ; cellfun(@(x) strcat('EC_',num2str(x)),num2cell(chn_num),'uni',0)' ];

% SC
chn_num = [17:32];

hld_dta.dat = [ hld_dta.dat ; dta_itn_dsm(chn_num,:) ];
hld_dta.lab = [ hld_dta.lab ; cellfun(@(x) strcat('SC_',num2str(x)),num2cell(chn_num),'uni',0)' ];

% HC
chn_num = [33:48];

hld_dta.dat = [ hld_dta.dat ; dta_itn_dsm(chn_num,:) ];
hld_dta.lab = [ hld_dta.lab ; cellfun(@(x) strcat('HC_',num2str(x)),num2cell(chn_num),'uni',0)' ];

% STG
chn_num = [49:64];

hld_dta.dat = [ hld_dta.dat ; dta_itn_dsm(chn_num,:) ];
hld_dta.lab = [ hld_dta.lab ; cellfun(@(x) strcat('STG_',num2str(x)),num2cell(chn_num),'uni',0)' ];

%% FieldTrip Conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to fieldtrip format
cont_data.fsample    = hld_dta.Fs;

cont_data.sampleinfo = [ 1 dta_ind_hld(2)-dta_ind_hld(1)+1 ];
cont_data.time     = {linspace(0,cont_data.sampleinfo(2)*1/cont_data.fsample,cont_data.sampleinfo(2))};
cont_data.trial{1} = hld_dta.dat;

cont_data.label = hld_dta.lab(:);

% .hdr
cont_data.hdr.Fs = cont_data.fsample;
cont_data.hdr.nSamples = cont_data.sampleinfo(2);
cont_data.hdr.nSamplesPre = 0;
cont_data.hdr.nTrials = 1;
cont_data.hdr.nChans = size(cont_data.trial{1},1);

% .hdr.orig
cont_data.hdr.orig.file_id    = hld_dta.file;
cont_data.hdr.orig.DataPoints = numel(cont_data.time{1});
cont_data.hdr.orig.DataDurationSec = cont_data.time{1}(end);
cont_data.hdr.orig.isaverage  = 0;
cont_data.hdr.orig.iscontinuous = 1;
cont_data.hdr.orig.isepoched  = 0;

% .cfg
cont_data.cfg.dataset      = cont_data.hdr.orig.file_id;
cont_data.cfg.precision    = 'double';
cont_data.cfg.trl = [1 cont_data.sampleinfo(2) 0];

% .cfg.callinfo
cont_data.cfg.callinfo.usercfg = sprintf('dataset: %s',cont_data.cfg.dataset);
cont_data.cfg.callinfo.inputhash = {};
cont_data.cfg.callinfo.fieldtrip = ft_version;
cont_data.cfg.callinfo.matlab = version;
cont_data.cfg.callinfo.computer = computer;
[~,user] = system('whoami');
cont_data.cfg.callinfo.user = user;
cont_data.cfg.callinfo.outputhash = {};

out_dat.data_name{1} = ['n' mmil_spec_char(fcfg.cln_fld{1},{'-'}) '_' num2str(1)];
out_dat.(out_dat.data_name{1}) = cont_data; %clear cont_data;

%% Channel Naming Convention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/']); mkdir([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/']); end
save([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample_micro' '.mat'],'out_dat','-v7.3');

end