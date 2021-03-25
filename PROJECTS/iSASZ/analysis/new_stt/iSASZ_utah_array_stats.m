function iSASZ_utah_array_stats(fcfg)

%% Load
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' 'MG049_SA_SZ' '_overall_data.mat'];
evt_dat  = ft_func([],cfg);

evt_key = mmil_readtext([ fcfg.clr_fld '/' 'events' '/' 'MG049_SA_SZ' '/' 'MG049_SA_SZ.csv' ]);
evt_key = cell2mat(evt_key(2:end,1));

%% Events
% Put together events
evt_lbl = fieldnames(evt_dat.(evt_dat.data_name{1}).cfg.alt_eve);

for iD = 1:numel(bcc_dat.data_name)
    for iE = 1:numel(evt_lbl)
        for iT = 1:numel(evt_dat.(evt_dat.data_name{1}).trialinfo)
            bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.(evt_lbl{iE})(iT) = evt_dat.(evt_dat.data_name{1}).cfg.alt_eve.(evt_lbl{iE})(iT);
        end
    end
end

% Events
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

%% Overall Channels
% Stats
cfg = [];
cfg.stt_fnc  = {'sasz_vis_ovr_new' 'sasz_aud_ovr_new'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% % Plot
% for iT = 1:numel(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse)
%     if ~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT)) & isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT))
%         ttt(iT) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT);
%     elseif isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT)) & ~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT))
%         ttt(iT) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT);
%     else
%         ttt(iT) = nan;
%     end
% end
% bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.plt_new_bse = ttt;
% 
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {bcc_dat.(bcc_dat.data_name{1})};
% 
% cfg.alt_eve   = 'plt_new_bse';
% cfg.eve       = [102 202];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('red') rgb('blue')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.stt_dat = repmat({{ 'vis_new_ovr_stt'      'aud_new_ovr_stt'}},10,10);
% cfg.stt_col = repmat({{ ft_stt_col(rgb('red')) ft_stt_col(rgb('blue'))}},10,10);
% cfg.stt_cmp = repmat({{ '0%3'                  '3%6'               }},10,10);
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('black')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot' '/'];
% cfg.prefix     = ['MG049_SASZ' '_' 'overall' '_' 'stt'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)

%% Repetition
% Stats
cfg = [];
cfg.stt_fnc  = {'sasz_vis_old_new' 'sasz_aud_old_new'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Mask

% % Plot
% for iT = 1:numel(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse)
%     if ~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT)) & isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT))
%         ttt(iT) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT);
%     elseif isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_new_bse(iT)) & ~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT))
%         ttt(iT) = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.aud_new_bse(iT);
%     else
%         ttt(iT) = nan;
%     end
% end
% bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.plt_new_bse = ttt;
% 
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {bcc_dat.(bcc_dat.data_name{1})};
% 
% cfg.alt_eve   = 'plt_new_bse';
% cfg.eve       = [102 202];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('red') rgb('blue')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.stt_dat = repmat({{ 'vis_new_ovr_stt'      'aud_new_ovr_stt'}},10,10);
% cfg.stt_col = repmat({{ ft_stt_col(rgb('red')) ft_stt_col(rgb('blue'))}},10,10);
% cfg.stt_cmp = repmat({{ '0%3'                  '3%6'               }},10,10);
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('black')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot' '/'];
% cfg.prefix     = ['MG049_SASZ' '_' 'overall' '_' 'stt'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)

%% Frequency
% Stats
cfg = [];
cfg.stt_fnc  = {'sasz_vis_frq' 'sasz_aud_frq'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% % Mask
% cfg     = [];
% cfg.stt     = { 'vis_new_old_stt' 'aud_new_old_stt'}; %
% cfg.stt_msk = { 'vis_new_ovr_stt' 'aud_new_ovr_stt' }; %
% cfg.pst_fix = '_msk_ovr';
% bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
%     
% % Vis Plot
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {bcc_dat.(bcc_dat.data_name{1})};
% 
% cfg.alt_eve   = 'vis_new_old';
% cfg.eve       = [111 112];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('red')  rgb('orangish red')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.stt_dat = repmat({{ 'vis_new_ovr_stt' 'vis_new_old_stt_msk_ovr'  }},10,10);
% cfg.stt_col = repmat({{ ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('bright red')) }},10,10);
% cfg.stt_cmp = repmat({{ '0%3'                  '3%6'               }},10,10);
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('black')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot' '/'];
% cfg.prefix     = ['MG049_SASZ' '_' 'vis' '_' 'repetition' '_' 'stt'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)
% 
% % Aud Plot
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {bcc_dat.(bcc_dat.data_name{1})};
% 
% cfg.alt_eve   = 'aud_new_old';
% cfg.eve       = [211 212];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('blue')  rgb('bluish purple')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.stt_dat = repmat({{ 'aud_new_ovr_stt' 'aud_new_old_stt_msk_ovr'  }},10,10);
% cfg.stt_col = repmat({{ ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('bright blue')) }},10,10);
% cfg.stt_cmp = repmat({{ '0%3'                  '3%6'               }},10,10);
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('black')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot' '/'];
% cfg.prefix     = ['MG049_SASZ' '_' 'aud' '_' 'repetition' '_' 'stt'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)

%% Semantic
% Stats
cfg = [];
cfg.stt_fnc  = {'sasz_vis_ani_obj' 'sasz_aud_ani_obj'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Mask

% Plot

%% Targetted Semantic
% Stats
    cfg = [];
    cfg.stt_fnc   = { 'iSASZ_vis_sem_mot_mnp' 'iSASZ_aud_sem_mot_mnp' };
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name;
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Mask

% Plot


%% Save Data
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end