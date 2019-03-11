function [out_dat,trg_chn,fcfg] = mmil_pedot_raw_load(fcfg,hld_dat)

fcfg.in_dir  = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir');
fcfg.cln_fld = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_fld');
fcfg.cln_dir = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir');
fcfg.end_dir = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir');

%% Intan Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(fcfg.end_dir{1}(end-2:end),'rhd')
    
    tot_fle = dir([fcfg.in_dir{1} '/' fcfg.cln_fld{1} '/' '*' fcfg.end_dir{1}]); tot_fle = {tot_fle(:).name};
    
    for iTF = 1:numel(fcfg.cln_dir);
        
        fle_hld{iTF} = tot_fle(string_find(tot_fle,fcfg.cln_dir{iTF}));
        
        lod_dta = read_Intan_RHD2000_file(fle_hld{iTF}{1},[fcfg.in_dir{1} '/' fcfg.cln_fld{1} '/' ]);
        
        dta_hld{iTF}.data     = lod_dta.amplifier_data;
        dta_hld{iTF}.trg_chn  = lod_dta.board_adc_data;
        dta_hld{iTF}.idx      = zeros(size(lod_dta.amplifier_data,1),1);
        dta_hld{iTF}.file     = fle_hld{iTF}{1};
        dta_hld{iTF}.idx_name = {lod_dta.amplifier_channels(:).native_channel_name}';
        dta_hld{iTF}.fs       = lod_dta.frequency_parameters.amplifier_sample_rate;
        dta_hld{iTF}.time     = 1:size(lod_dta.amplifier_data,2);
        dta_hld{iTF}.imp      = [lod_dta.amplifier_channels.electrode_impedance_magnitude]';
        dta_hld{iTF}.cond     = fcfg.cln_dir{iTF};
        
        clear lod_dta
        
        for iD = 2:numel(fle_hld{iTF})
            
            lod_dta = read_Intan_RHD2000_file(fle_hld{iTF}{iD},[fcfg.in_dir{1} '/' fcfg.cln_fld{1} '/' ]);
            
            dta_hld{iTF}.data     = [dta_hld{iTF}.data lod_dta.amplifier_data];
            dta_hld{iTF}.trg_chn  = [dta_hld{iTF}.trg_chn lod_dta.board_adc_data];
            dta_hld{iTF}.time     = [dta_hld{iTF}.time dta_hld{iTF}.time(end)+1:dta_hld{iTF}.time(end)+size(lod_dta.amplifier_data,2)];
            
            clear lod_dta
            
        end
    end
    
end

%% Open EPhys Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(fcfg.end_dir{1}(end-9:end),'continuous')
    for iCL = 1:numel(fcfg.cln_fld)
     
        tot_fle = dir([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' '*' fcfg.end_dir{1}]); tot_fle = {tot_fle(:).name};
        tot_fle_chn = tot_fle(string_find(tot_fle,'_CH'));
        tot_fle_adc = tot_fle(string_find(tot_fle,'_ADC'));
        tot_fle_aux = tot_fle(string_find(tot_fle,'_AUX'));
    
        open_ephys_get_session_info([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/']);
        
        for iCH = 1:numel(tot_fle_chn);  lod_dta(iCH,:) = load_open_ephys_data([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' tot_fle_chn{iCH}]); end
        
    end
end

%% FieldTrip Conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iF = 1:numel(fcfg.cln_dir)
    
    trg_chn{iF} = hld_dat{iF}.trg_chn;
    
    dsfact = 20;
    anti_alias_factor = 2.5;
    filter_order = 4;
    
    %% Change to fieldtrip format
    cont_data.fsample    = hld_dat{iF}.fs;
    
    cont_data.sampleinfo = [1 size(hld_dat{iF}.data,2)];
    cont_data.time     = {linspace(0,cont_data.sampleinfo(2)*1/cont_data.fsample,cont_data.sampleinfo(2))};
    cont_data.trial{1} = hld_dat{iF}.data;
    cont_data.trial{1}(end+1,:) = hld_dat{iF}.trg_chn;
    
    cont_data.label = [hld_dat{1}.idx_name(:) ; {'trigger'}];
    
    % .hdr
    cont_data.hdr.Fs = cont_data.fsample;
    cont_data.hdr.nSamples = cont_data.sampleinfo(2);
    cont_data.hdr.nSamplesPre = 0;
    cont_data.hdr.nTrials = 1;
    cont_data.hdr.nChans = size(hld_dat{iF}.data,1);
    
    % .hdr.orig
    cont_data.hdr.orig.file_id    = hld_dat{iF}.file;
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
    
    out_dat.data_name{iF} = fcfg.cln_dir{iF};
    out_dat.(fcfg.cln_dir{iF}) = cont_data; %clear cont_data;
    
end

clear hld_dat

%% Downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Downsampling...\n')

for iF = 1:numel(fcfg.cln_dir)
    
    dat_lng = out_dat.(out_dat.data_name{iF}).sampleinfo(end);
    Fs = out_dat.(out_dat.data_name{iF}).fsample;
    [b,a] = butter( filter_order, ((Fs/(dsfact*anti_alias_factor))/(Fs/2)), 'low' );
    
    for iC = 1:size(out_dat.(out_dat.data_name{iF}).trial{1},1)
        fprintf('Downsampling Channel %i \n',iC)
        xf = filtfilt( b, a, out_dat.(out_dat.data_name{iF}).trial{1}(iC,:));
        chn_out(iC,:) = downsample( xf, dsfact );
    end
    
    out_dat.(out_dat.data_name{iF}).trial{1} = chn_out;
    out_dat.(out_dat.data_name{iF}).sampleinfo(end) = out_dat.(out_dat.data_name{iF}).sampleinfo(end) / dsfact;
    out_dat.(out_dat.data_name{iF}).fsample  = out_dat.(out_dat.data_name{iF}).fsample/dsfact;
    out_dat.(out_dat.data_name{iF}).time{1}  = downsample( out_dat.(out_dat.data_name{iF}).time{1}, dsfact );
    
    clear chn_out
    
end
    
save([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample' '.mat'],'out_dat','trg_chn','-v7.3');

end