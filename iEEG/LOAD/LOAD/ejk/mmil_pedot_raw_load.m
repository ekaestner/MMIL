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
     
        %%%%%%%%%%
        tot_fle = dir([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' '*' fcfg.end_dir{1}]); tot_fle = {tot_fle(:).name};
        tot_fle_chn = tot_fle(string_find(tot_fle,'_CH'));
        tot_fle_adc = tot_fle(string_find(tot_fle,'_ADC'));
        tot_fle_aux = tot_fle(string_find(tot_fle,'_AUX'));
        
        %%%%%%%%%%
        for iCH = 1:numel(tot_fle_chn); chn_num(iCH) = str2num(tot_fle_chn{iCH}(strfind(tot_fle_chn{iCH},'CH')+2:strfind(tot_fle_chn{iCH},'.continuous')-1)); end
        [~,srt_ind] = sort(chn_num);
        tot_fle_chn = tot_fle_chn(srt_ind); clear chn_num;
        
        for iCH = 1:numel(tot_fle_adc); chn_num(iCH) = str2num(tot_fle_adc{iCH}(strfind(tot_fle_adc{iCH},'ADC')+3:strfind(tot_fle_adc{iCH},'.continuous')-1)); end
        [~,srt_ind] = sort(chn_num);
        tot_fle_adc = tot_fle_adc(srt_ind); clear chn_num;
        
        for iCH = 1:numel(tot_fle_aux); chn_num(iCH) = str2num(tot_fle_aux{iCH}(strfind(tot_fle_aux{iCH},'AUX')+3:strfind(tot_fle_aux{iCH},'.continuous')-1)); end
        [~,srt_ind] = sort(chn_num);
        tot_fle_aux = tot_fle_aux(srt_ind); clear chn_num;
        
        %%%%%%%%%%        
        for iCH = 1:numel(tot_fle_chn); lod_dta(iCH,:) = load_open_ephys_data([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' tot_fle_chn{iCH}]); end
        for iCH = 1:numel(tot_fle_adc); adc_dta(iCH,:) = load_open_ephys_data([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' tot_fle_adc{iCH}]); end
        for iCH = 1:numel(tot_fle_aux); aux_dta(iCH,:) = load_open_ephys_data([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' tot_fle_aux{iCH}]); end
        
        mta_dta = mmil_readtext([fcfg.in_dir{1} '/' fcfg.cln_fld{iCL} '/' 'settings.xml'],['\t']);
        mta_dta = mta_dta{string_find(mta_dta,{'SampleRateString'})};
        mta_dta_str = strfind(mta_dta,'SampleRateString="')+numel('SampleRateString="'); mta_dta_end = strfind(mta_dta,' kS/s"');
        frq_smp = str2num(mta_dta(mta_dta_str:mta_dta_end))*1000;
        
        dta_hld{iCL}.data     = lod_dta;
        dta_hld{iCL}.trg_chn  = adc_dta;
        dta_hld{iCL}.idx      = zeros(size(lod_dta,1),1);
        dta_hld{iCL}.file     = tot_fle_chn{iCL};
        dta_hld{iCL}.idx_name = cellfun(@(x) x(5:end-11),tot_fle_chn,'uni',0)';
        dta_hld{iCL}.fs       = frq_smp; % HARDCODED NOT GREAT :(
        dta_hld{iCL}.time     = 1:size(lod_dta,2);
        dta_hld{iCL}.imp      = []; % ??????????????????????
        dta_hld{iCL}.cond     = fcfg.cln_fld{iCL};
         
        clear lod_dta adc_dta aux_dta

    end
end

%% FieldTrip Conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iF = 1:numel(fcfg.cln_fld)
    
    trg_chn{iF} = dta_hld{iF}.trg_chn;
    
    dsfact = dta_hld{iCL}.fs / 1000;
    anti_alias_factor = 2.5;
    filter_order = 4;
    
    cont_data.dsfact = dta_hld{iCL}.fs / 1000;
    
    %% Change to fieldtrip format
    cont_data.fsample    = dta_hld{iF}.fs;
    
    cont_data.sampleinfo = [1 size(dta_hld{iF}.data,2)];
    cont_data.time     = {linspace(0,cont_data.sampleinfo(2)*1/cont_data.fsample,cont_data.sampleinfo(2))};
    cont_data.trial{1} = dta_hld{iF}.data;
%     cont_data.trial{1}(end+1,:) = dta_hld{iF}.trg_chn;
    
    cont_data.label = dta_hld{1}.idx_name(:);
    
    % .hdr
    cont_data.hdr.Fs = cont_data.fsample;
    cont_data.hdr.nSamples = cont_data.sampleinfo(2);
    cont_data.hdr.nSamplesPre = 0;
    cont_data.hdr.nTrials = 1;
    cont_data.hdr.nChans = size(dta_hld{iF}.data,1);
    
    % .hdr.orig
    cont_data.hdr.orig.file_id    = dta_hld{iF}.file;
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
    
    out_dat.data_name{iF} = ['n' mmil_spec_char(fcfg.cln_fld{iF},{'-'}) '_' num2str(iF)];
    out_dat.(out_dat.data_name{iF}) = cont_data; %clear cont_data;
    
end

clear dta_hld

%% Downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Downsampling...\n')

for iF = 1:numel(out_dat.data_name)
    
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