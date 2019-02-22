function [out_dat,trg_chn] = mmil_pedot_raw_load(fcfg,hld_dat)

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

%% Downsample
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