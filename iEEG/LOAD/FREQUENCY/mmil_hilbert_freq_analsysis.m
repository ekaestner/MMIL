function [ft_dat] = ft_hilbert_freq_analsysis(fcfg,ft_dat)
%
% Outputs
% freq_ft_dat: Normal frequency structure (based on the structures created by ) 
% 
% Inputs
% fcfg.freq_band: 1 x 2 matrix of cutoff bands for frequency data, example: ([low_cf high_cf])
% fcfg.hilbert:   'amp' (default) or 'angle'
% 
% CURRENTLY USES 4TH ORDER 

if size(fcfg.freq_band,1) ~= 1 && size(fcfg.freq_band,2) ~= 2; error('Only run one band at a time'); end

for iFR = 1:numel(fcfg.freq_band)
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = fcfg.freq_band{iFR};
    tmp_dat      = ft_func_v2(@ft_preprocessing,cfg,ft_dat);
    
    if strcmpi(fcfg.hilbert,'amp')
        cfg.hilbert  = 'abs';
        tmp_dat_amp  = ft_func_v2(@ft_preprocessing,cfg,tmp_dat);
        tmp_dat      = cat(3,tmp_dat_amp.trial{:});
    elseif strcmpi(fcfg.hilbert,'phs')
        cfg.hilbert  = 'complex';
        tmp_dat_phs  = ft_func_v2(@ft_preprocessing,cfg,tmp_dat);
        tmp_dat      = cat(3,tmp_dat_phs.trial{:});
    end
               
    frq_dat(:,:,iFR,:) = permute(tmp_dat,[1 4 3 2]);
    
end

if numel(ft_dat.trial) == 1
    ft_dat.trial  = {frq_dat};
else
    ft_dat.trial  = cellfun(@squeeze,num2cell(frq_dat,[1 3 4]),'uni',0); %% - EJK need to test whether this adds to previously created frequencies
end

ft_dat.freq   = fcfg.freq_band; %% - EJK need to test whether this adds to previously created frequencies
ft_dat.type   = fcfg.hilbert;
ft_dat.dimord = '{trial}(chan,freq,time)';
 
end