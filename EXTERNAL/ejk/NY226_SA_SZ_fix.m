%% Fix the problems with NY226 and convert into FT
function sem_dat = NY226_SA_SZ_fix(tmp_dat)

for iDT = 1:numel(tmp_dat.data_name)
    
    switch tmp_dat.data_name{iDT}
        
        case 'NY226_091214_SA_03_500Hz_flt32_nspike_epoch_data'
            
            sem_dat.data_name{1} = 'NY226_091214_SA_03_500Hz_flt32_nspike_epoch_data';
            
            % convert data
            con_dat.trial = {};
            con_dat.trialinfo = [];
            con_dat.fsample   = tmp_dat.(tmp_dat.data_name{1}).sfreq;
            con_dat.label     = {tmp_dat.(tmp_dat.data_name{1}).sensor_info.label}';
            
            for iPR = 1:numel(tmp_dat.(tmp_dat.data_name{1}).epochs)
                con_dat.trial = [con_dat.trial(:) ; squeeze(mat2cell(tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).data,size(tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).data,1),size(tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).data,2),ones(1,size(tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).data,3))))];
                con_dat.trialinfo = [con_dat.trialinfo ; repmat(tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).event_code,tmp_dat.(tmp_dat.data_name{1}).epochs(iPR).num_trials,1)];
            end
            
            con_dat.sampleinfo = [1 size(con_dat.trial{1},2)];
            con_dat.time{1}    = tmp_dat.(tmp_dat.data_name{1}).epochs(1).time;
            for iPR = 2:numel(con_dat.trial)
                con_dat.sampleinfo(end+1,:) = [con_dat.sampleinfo(end,1)+1 con_dat.sampleinfo(end,1)+1+numel(con_dat.trial{1})];
                con_dat.time{iPR}           = tmp_dat.(tmp_dat.data_name{1}).epochs(1).time;
            end
            
            sem_dat.(sem_dat.data_name{1}) = con_dat;
            clear con_dat
            sem_dat.(sem_dat.data_name{1}).cfg.orig_lab = sem_dat.(sem_dat.data_name{1}).label;
            
            sem_dat.(sem_dat.data_name{1}).trialinfo = sem_dat.(sem_dat.data_name{1}).trialinfo + 10;
            
            cfg = [];
            cfg.latency = [-0.2 0.8];
            sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
            
        case 'NY226_091214_SZ_SAblocks_500to1000_500Hz_flt32_nspike_epoch_data'
            
            sem_dat.data_name{2} = 'NY226_091214_SZ_SAblocks_500to1000_500Hz_flt32_nspike_epoch_dat';
            
            % convert data
            con_dat.trial = {};
            con_dat.trialinfo = [];
            con_dat.fsample   = tmp_dat.(tmp_dat.data_name{2}).sfreq;
            con_dat.label     = {tmp_dat.(tmp_dat.data_name{2}).sensor_info.label}';
            
            for iPR = 1:numel(tmp_dat.(tmp_dat.data_name{2}).epochs)
                con_dat.trial = [con_dat.trial(:) ; squeeze(mat2cell(tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).data,size(tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).data,1),size(tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).data,2),ones(1,size(tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).data,3))))];
                con_dat.trialinfo = [con_dat.trialinfo ; repmat(tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).event_code,tmp_dat.(tmp_dat.data_name{2}).epochs(iPR).num_trials,1)];
            end
            
            con_dat.sampleinfo = [1 size(con_dat.trial{2},2)];
            con_dat.time{1}    = tmp_dat.(tmp_dat.data_name{2}).epochs(1).time;
            for iPR = 2:numel(con_dat.trial)
                con_dat.sampleinfo(end+1,:) = [con_dat.sampleinfo(end,1)+1 con_dat.sampleinfo(end,1)+1+numel(con_dat.trial{2})];
                con_dat.time{iPR}           = tmp_dat.(tmp_dat.data_name{2}).epochs(1).time;
            end
            
            sem_dat.(sem_dat.data_name{2}) = con_dat;
            clear con_dat
            sem_dat.(sem_dat.data_name{2}).cfg.orig_lab = sem_dat.(sem_dat.data_name{2}).label;
            
            sem_dat.(sem_dat.data_name{2}).trialinfo = sem_dat.(sem_dat.data_name{2}).trialinfo + 10;
            
        case 'NY226_091214_SZ_SZblocks_500to1000_500Hz_flt32_nspike_epoch_data'
            
            sem_dat.data_name{3} = 'NY226_091214_SZ_SZblocks_500to1000_500Hz_flt32_nspike_epoch_dat';
            
            % convert data
            con_dat.trial = {};
            con_dat.trialinfo = [];
            con_dat.fsample   = tmp_dat.(tmp_dat.data_name{3}).sfreq;
            con_dat.label     = {tmp_dat.(tmp_dat.data_name{3}).sensor_info.label}';
            
            end_ind = [50 50 49 50 60 40 40 0];
            for iPR = 1:numel(tmp_dat.(tmp_dat.data_name{3}).epochs)
                con_dat.trial = [con_dat.trial(:) ; squeeze(mat2cell(tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).data,size(tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).data,1),size(tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).data,2),ones(1,size(tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).data,3))))];
                
                tmp_trl           = repmat(tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).event_code,tmp_dat.(tmp_dat.data_name{3}).epochs(iPR).num_trials,1);
                tmp_trl(end_ind(iPR)+1:end) = tmp_trl(end_ind(iPR)+1:end)+10;
                con_dat.trialinfo = [con_dat.trialinfo ; tmp_trl];
                
            end
            
            con_dat.sampleinfo = [1 size(con_dat.trial{3},2)];
            con_dat.time{1}    = tmp_dat.(tmp_dat.data_name{3}).epochs(1).time;
            for iPR = 2:numel(con_dat.trial)
                con_dat.sampleinfo(end+1,:) = [con_dat.sampleinfo(end,1)+1 con_dat.sampleinfo(end,1)+1+numel(con_dat.trial{3})];
                con_dat.time{iPR}           = tmp_dat.(tmp_dat.data_name{3}).epochs(1).time;
            end
            
            sem_dat.(sem_dat.data_name{3}) = con_dat;
            clear con_dat
            sem_dat.(sem_dat.data_name{3}).cfg.orig_lab = sem_dat.(sem_dat.data_name{3}).label;
            
            cfg = [];
            cfg.latency = [-0.2 0.8];
            sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
    end
    
end

cfg           = [];
cfg.data_name = [1 ; 3];
cfg.data_new  = 'yes';
cfg.methapp   = 'trials';
sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);

cfg = [];
cfg.rmfield   = 'yes';
cfg.data_name = [2 3];
sem_dat = ft_func([],cfg,sem_dat);

end
