function iSASZ_trg_sem_eve(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile  = [fcfg.dat_fld '/' sbj '_overall_data.mat'];
outpath = [fcfg.dat_fld '/' ];

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Events
ttt_sem = mmil_readtext([fcfg.clr_fld '/' 'task' '/' 'semantic_database.csv']);

for iD = 1:numel(bcc_dat.data_name)
    for iPS = [6 8]
        for iT = 1:size(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn)
            ttt = find(strcmpi(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn{iT},ttt_sem(:,1)));
            if ~isempty(ttt) && ( bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==1 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==2 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==3 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==4 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==11 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==12 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==13 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==14 )
                bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.(ttt_sem{1,iPS})(iT,1) = cell2mat(ttt_sem(ttt,iPS));
            else
                bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.(ttt_sem{1,iPS})(iT,1) = nan;
            end
        end
    end
end

for iD = 1:numel(bcc_dat.data_name)
    for iT = 1:size(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.stm_idn)
        if ~isnan(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.motion(iT))
            if     bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==1 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==3
                if bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.motion(iT)>2 && bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.manipulation(iT)<2;
                    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = 901;
                end
            elseif bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==2 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==4
                if bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.motion(iT)<2 && bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.manipulation(iT)>2;
                    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = 902;
                end
            elseif bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==11 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==13
                if bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.motion(iT)>2 && bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.manipulation(iT)<2;
                    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = 911;
                end
            elseif bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==12 || bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo(iT)==14
                if bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.motion(iT)<2 && bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.manipulation(iT)>2;
                    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = 912;
                end
            else
                bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = nan;
            end
        else
            bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.sem_mot_mnp(iT) = nan;
        end
    end
end

%% Stats
eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

% Visual
if any(eve_typ < 9)
    cfg = [];
    cfg.stt_fnc   = { 'iSASZ_vis_sem_mot_mnp' };
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name;
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

% Auditory
if any(eve_typ > 9)
    cfg = [];
    cfg.stt_fnc   = { 'iSASZ_aud_sem_mot_mnp' };
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name;
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

%% Mask
eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ < 9)
    cfg     = [];
    cfg.stt     = { 'vis_sem_mot_mnp_stt' }; %
    cfg.stt_msk = { 'vis_new_ovr_stt' }; %
    cfg.pst_fix = '_msk_ovr';
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

if any(eve_typ > 9)
    cfg     = [];
    cfg.stt     = { 'aud_sem_mot_mnp_stt' }; %
    cfg.stt_msk = { 'aud_new_ovr_stt' }; %
    cfg.pst_fix = '_msk_ovr';
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end