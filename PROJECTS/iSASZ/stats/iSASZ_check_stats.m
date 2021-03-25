clear; clc;

chn_chs = 1;
mve_plt = 1;

for sbj_num = [1:17 19:numel(sbj_hld)];
    
    % Setting up Variables
    subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
    subj  = subj{sbj_num};
    
    fprintf('Starting work on %s \n',subj)
    
    infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/' subj '_overall_data.mat'];
    outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = infile;
    sem_dat  = ft_func([],cfg);
    
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    %% Choose Channels based on hypotheses
    if chn_chs == 1;
        
        if has_vis && has_aud
            
            cfg = [];
            %         cfg.ovr_wrt = 1;
            cfg.alt_stt = {'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_dif_rep_stt' 'vis_dif_rep_stt' 'vis_dif_sem_stt' 'vis_dif_sem_stt' ...
                'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_dif_rep_stt' 'aud_dif_rep_stt' 'aud_dif_sem_stt' 'aud_dif_sem_stt'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark orange')) ft_stt_col(rgb('dark orange')) ...
                [0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark orange')) ft_stt_col(rgb('dark orange'))};
            cfg.cmp_stt = [1 2 3 4  5  6 ...
                7 8 9 10 11 12];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'sem_all_eve' 'sem_all_eve' ...
                'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'sem_all_eve' 'sem_all_eve' };
            cfg.cmp     = {'101!999' '101!999' '111!112' '111!112' '121!122' '121!122'...
                '102!999' '102!999' '113!114' '113!114' '123!124' '123!124'};
            cfg.cmp_nme = {'vis_ovr_all_erl' 'vis_ovr_all_lte' 'vis_dif_rep_erl' 'vis_dif_rep_lte' 'vis_dif_sem_erl' 'vis_dif_sem_lte' ...
                'aud_ovr_all_erl' 'aud_ovr_all_lte' 'aud_dif_rep_erl' 'aud_dif_rep_lte' 'aud_dif_sem_erl' 'aud_dif_sem_lte'};
            cfg.tme_win = {[0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200] ...
                [0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200]};
            cfg.clr_fld = outpath; %'/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        elseif has_vis
            
            cfg = [];
            %         cfg.ovr_wrt = 1;
            cfg.alt_stt = {'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_dif_rep_stt' 'vis_dif_rep_stt' 'vis_dif_sem_stt' 'vis_dif_sem_stt'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark orange')) ft_stt_col(rgb('dark orange'))};
            cfg.cmp_stt = [1 2 3 4  5  6];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'sem_all_eve' 'sem_all_eve'};
            cfg.cmp     = {'101!999' '101!999' '111!112' '111!112' '121!122' '121!122'};
            cfg.cmp_nme = {'vis_ovr_all_erl' 'vis_ovr_all_lte' 'vis_dif_rep_erl' 'vis_dif_rep_lte' 'vis_dif_sem_erl' 'vis_dif_sem_lte'};
            cfg.tme_win = {[0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200]};
            cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        elseif has_aud
            
            cfg = [];
            %         cfg.ovr_wrt = 1;
            cfg.alt_stt = {'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_dif_rep_stt' 'aud_dif_rep_stt' 'aud_dif_sem_stt' 'aud_dif_sem_stt'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark orange')) ft_stt_col(rgb('dark orange'))};
            cfg.cmp_stt = [1 2 3 4  5  6];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'sem_all_eve' 'sem_all_eve' };
            cfg.cmp     = {'102!999' '102!999' '113!114' '113!114' '123!124' '123!124'};
            cfg.cmp_nme = {'aud_ovr_all_erl' 'aud_ovr_all_lte' 'aud_dif_rep_erl' 'aud_dif_rep_lte' 'aud_dif_sem_erl' 'aud_dif_sem_lte'};
            cfg.tme_win = {[0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200] [0.100 0.500] [0.500 1.200]};
            cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        end
    end
    %% Correct stats based on SNR
    % tme_int = [0.100 0.500];  [tmp,tme_int(1)] = min(abs(tme_int(1) - sem_dat.(sem_dat.data_name{iMV}).time{1})); [tmp,tme_int(2)] = min(abs(tme_int(2) - sem_dat.(sem_dat.data_name{iMV}).time{1})); tme_int = sem_dat.(sem_dat.data_name{iMV}).time{1}(tme_int(1):tme_int(2));
    % tme_bse = [-0.400 0.000]; [tmp,tme_bse(1)] = min(abs(tme_bse(1) - sem_dat.(sem_dat.data_name{iMV}).time{1})); [tmp,tme_bse(2)] = min(abs(tme_bse(2) - sem_dat.(sem_dat.data_name{iMV}).time{1}));
    %
    % for iMV = 1:numel(sem_dat.data_name)
    %
    %     sig_fle = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}]);
    %
    %     ovr_dat_hld = cat(3,sem_dat.(sem_dat.data_name{iMV}).trial{:});
    %
    %     if has_vis && has_aud
    %         sig_col = [1 7];
    %     elseif has_vis
    %         sig_col = 1;
    %     elseif has_vis
    %         sig_col = 1;
    %     end
    %
    %     for iCL = 1:numel(sig_col)
    %         sig_ind     = find(cell2mat(sig_fle(2:end,sig_col(iCL)+1)))+1;
    %         for iCH = 1:numel(sig_ind); sig_chn_ind(iCH) = str2num(sig_fle{sig_ind(iCH),1}(1:strfind(sig_fle{sig_ind(iCH),1},' ')-1)) ; end
    %         for iCH = 1:numel(sig_ind); sig_prb_ind(iCH) = find(strcmpi(sem_dat.(sem_dat.data_name{iMV}).label(sig_chn_ind(iCH)),sem_dat.(sem_dat.data_name{iMV}).cfg.alt_stt.vis_ovr_all_stt.label)) ; end
    %
    %         for iSG = 1:numel(sig_chn_ind)
    %
    %             % calculate snr over the timecourse
    %             dat_hld = squeeze(ovr_dat_hld(sig_chn_ind(iSG),:,:))';
    %             for iTM = 1:size(dat_hld,2)
    %                 snr_hld(iTM) = nanmean(dat_hld(:,iTM)) / nanvar(dat_hld(:,iTM));
    %             end
    %
    %             % calculate snr for baseline
    %             bse_snr = mean(abs(snr_hld(tme_bse(1):tme_bse(2))));
    %
    %             sig_inc(iSG,1) = sem_dat.(sem_dat.data_name{2}).label(sig_chn_ind(iSG));
    %
    %             % find sig periods of interest
    %             sig_hld = squeeze(sem_dat.(sem_dat.data_name{iMV}).cfg.alt_stt.vis_ovr_all_stt.prob(sig_prb_ind(iSG),:)) < 0.05;
    %             sig_hld(1) = 0;
    %             sig_hld(end) = 0;
    %             sig_prd = crossing(sig_hld(:),[],0.5);
    %             sig_1st = sig_prd(1:2:end);
    %             sig_2nd = sig_prd(2:2:end);
    %
    %             cnt = 1;
    %             for iSP = 1:numel(sig_1st)
    %                 tme_sig = sem_dat.(sem_dat.data_name{iMV}).time{1}(sig_1st(iSP):sig_2nd(iSP));
    %                 if ~isempty(intersect(tme_sig,tme_int))
    %                     sig_inc{iSG,2}(cnt) = mean(snr_hld(sig_1st(iSP):sig_2nd(iSP))) / bse_snr;
    %                     cnt = cnt + 1;
    %                 end
    %             end
    %
    %             % calculate snr for each relevant period of interest
    % %             for iSP = 1:numel(sig_1st)
    % %                 sig_snr(iSP) = ;
    % %             end
    %
    %             % make decision if snr is high enough
    %
    %         end
    %     end
    % end
    
    %% Move files around
    if mve_plt == 1
        
        outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
        
        for iMV = 1:numel(sem_dat.data_name)
            
            if has_vis && has_aud
                cfg = [];
                cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
                cfg.cmb     = 1;
                cfg.cmp_ind = [1 3 7 9];
                cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
                [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
            elseif has_vis
                cfg = [];
                cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
                cfg.cmb     = 1;
                cfg.cmp_ind = [1 3];
                cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
                [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
            elseif has_aud
                cfg = [];
                cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
                cfg.cmb     = 1;
                cfg.cmp_ind = [1 3];
                cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
                [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
            end
            
            % Copyfiles over
            fld = [outpath '/' 'initial_plot' '/' 'ThreeStat' '/' subj '/' sem_dat.data_name{iMV}(end-2:end)];
            
            for iCP = 1:numel(pre_fix)
                if exist([fld '/' pre_fix{iCP}],'dir'); rmdir([fld '/' pre_fix{iCP}],'s'); mkdir([fld '/' pre_fix{iCP}]); else mkdir([fld '/' pre_fix{iCP}]); end
                for iFL = 1:numel(grp_chn{iCP});
                    copyfile([fld '/' subj '__average_chan_' sem_dat.(sem_dat.data_name{iMV}).cfg.alt_lab.label{grp_chn{iCP}(iFL)} '.png'], ...
                        [fld '/' pre_fix{iCP} '/' subj '__average_chan_' sem_dat.(sem_dat.data_name{iMV}).cfg.alt_lab.label{grp_chn{iCP}(iFL)} '.png']);
                end
            end
            
        end
    end
    
    %% Count Channels
    outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
    
    for iMV = 1:numel(sem_dat.data_name)
        
        % Make list of overall numbers of electrodes
        sig_chn     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}]);
        sig_chn_int = sig_chn(2:end,[1 2 4 8 10]);
        
        if has_vis && has_aud
            ele_hld = {'subject','v_act_ely','v_rep_ely','a_act_ely', 'a_rep_ely','ovr_lap_ovr','ovr_lap_rep'};
        elseif has_vis
            ele_hld = {'subject','v_act_ely','v_rep_ely'};
        elseif has_aud
            ele_hld = {'subject','a_act_ely', 'a_rep_ely'};
        end
        
        ele_hld{2,1} = subj;
        
        ind_icd = 1:size(sig_chn_int,1);
        dpt_ind = find(~cellfun(@isempty,strfind(sig_chn_int(:,1),'D')) | ~cellfun(@isempty,strfind(sig_chn_int(:,1),'d')));
        if any(~cellfun(@isempty,strfind(sig_chn_int(dpt_ind,1),'>')))
            str_ind = strfind(sig_chn_int(dpt_ind,1),'>');
        else
            str_ind = strfind(sig_chn_int(dpt_ind,1),'-');
        end
        if ~isempty(dpt_ind); dpt_ind = dpt_ind(cellfun(@(x,y) strcmpi(x(y+1),'d'),sig_chn_int(dpt_ind,1),str_ind)); end
        ind_icd(dpt_ind-1) = [];
        
        if has_vis && has_aud
            
            ele_hld{2,2} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ele_hld{2,3} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            ele_hld{2,4} = sum(cell2mat(sig_chn_int(ind_icd,4)));
            ele_hld{2,5} = sum(cell2mat(sig_chn_int(ind_icd,5)));
            ele_hld{2,6} = sum(cell2mat(sig_chn_int(ind_icd,4)) & cell2mat(sig_chn_int(ind_icd,2)));
            ele_hld{2,7} = sum(cell2mat(sig_chn_int(ind_icd,5)) & cell2mat(sig_chn_int(ind_icd,3)));
            
        elseif has_vis
            
            ele_hld{2,2} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ele_hld{2,3} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            
        elseif has_aud
            
            ele_hld{2,2} = sum(cell2mat(sig_chn_int(ind_icd,2)));
            ele_hld{2,3} = sum(cell2mat(sig_chn_int(ind_icd,3)));
            
        end
        
        cell2csv(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/paper_plot/ThreeStat/' subj '/' sem_dat.data_name{iMV}(end-2:end) '/' sem_dat.data_name{iMV} '.csv'],ele_hld)
        
    end
    
end

% for i = 1:40
%    fprintf('%s : \n', sig_inc{i,1})
%     sig_inc{i,2}
% end
% ejk = 1;
% sig_inc{ejk,1}
% sig_inc{ejk,2}
%
% figure()
% plot(snr_hld,'r')
%
% figure()
% sig_ind     = find(cell2mat(sig_fle(2:end,sig_col(iCL)+1)))+1;
% for iCH = 1:numel(sig_ind); sig_chn_ind(iCH) = str2num(sig_fle{sig_ind(iCH),1}(1:strfind(sig_fle{sig_ind(iCH),1},' ')-1)) ; end
% for iCH = 1:numel(sig_ind); sig_prb_ind(iCH) = find(strcmpi(sem_dat.(sem_dat.data_name{iMV}).label(sig_chn_ind(iCH)),sem_dat.(sem_dat.data_name{iMV}).cfg.alt_stt.vis_ovr_all_stt.label)) ; end
% for iSG = 1:9
% subplot(3,3,iSG)
% dat_hld = squeeze(ovr_dat_hld(sig_chn_ind(iSG),:,:))';
% plot(nanmean(dat_hld),'k')
% for iEV = 1:400
% low_bnd(iEV) = nanmean(dat_hld(iEV,tme_bse(1):tme_bse(2))) + (nanstd(dat_hld(iEV,tme_bse(1):tme_bse(2))) / sqrt(numel(tme_bse(1):tme_bse(2)))) * 1.96;
% hgh_bnd(iEV) = nanmean(dat_hld(iEV,tme_bse(1):tme_bse(2))) - (nanstd(dat_hld(iEV,tme_bse(1):tme_bse(2))) / sqrt(numel(tme_bse(1):tme_bse(2)))) * 1.96;
% end
% hline(nanmean(hgh_bnd));
% hline(nanmean(low_bnd));
% end
%
% figure()
% sig_ind     = find(~cell2mat(sig_fle(2:end,sig_col(iCL)+1)))+1;
% for iCH = 1:numel(sig_ind); sig_chn_ind(iCH) = str2num(sig_fle{sig_ind(iCH),1}(1:strfind(sig_fle{sig_ind(iCH),1},' ')-1)) ; end
% for iCH = 1:numel(sig_ind); sig_prb_ind(iCH) = find(strcmpi(sem_dat.(sem_dat.data_name{iMV}).label(sig_chn_ind(iCH)),sem_dat.(sem_dat.data_name{iMV}).cfg.alt_stt.vis_ovr_all_stt.label)) ; end
% for iSG = 1:9
% subplot(3,3,iSG)
% dat_hld = squeeze(ovr_dat_hld(sig_chn_ind(iSG),:,:))';
% plot(nanmean(dat_hld),'k')
% for iEV = 1:400
% low_bnd(iEV) = nanmean(dat_hld(iEV,tme_bse(1):tme_bse(2))) + (nanstd(dat_hld(iEV,tme_bse(1):tme_bse(2))) / sqrt(numel(tme_bse(1):tme_bse(2)))) * 1.96;
% hgh_bnd(iEV) = nanmean(dat_hld(iEV,tme_bse(1):tme_bse(2))) - (nanstd(dat_hld(iEV,tme_bse(1):tme_bse(2))) / sqrt(numel(tme_bse(1):tme_bse(2)))) * 1.96;
% end
% hline(nanmean(hgh_bnd));
% hline(nanmean(low_bnd));
% end
%
%
%
% cfg = [];
% cfg.data_name = 2;
% cfg.events   = [101];
% cfg.add_stt  = 'vis_ovr_all_stt';
% cfg.alt_eve  = 'ovr_all_evt';
% cfg.alpha    = .05;
% sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
%
% cfg = [];
% cfg.data_name = 2;
% cfg.basetime = [-0.400 0];
% cfg.events   = [101];
% cfg.add_stt  = 'vis_ovr_all_stt_v3';
% cfg.alt_eve  = 'ovr_all_evt';
% cfg.alpha    = .05;
% sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
%
% iSG = 1;
% alpha = 0.05;
% dat_hld = squeeze(ovr_dat_hld(sig_chn_ind(iSG),:,:))';
% for iT = 1:size(dat_hld,2)
%    [tmp,pval_vec(iT)] = ttest( dat_hld(:,iT));
% end
% pval_sorted = sort(pval_vec);
% fdrpvec=((1:length(pval_vec))./length(pval_vec))*alpha;
% fdrtest = pval_sorted-fdrpvec;
% fdrp = find(fdrtest < 0,1,'last');
% fdrmask = zeros(1,length(pval_vec));
% if ~isempty(fdrp),
%     fdrthresh = pval_sorted(fdrp);
%     fdrmask(pval_vec < fdrthresh)=1;
% else
%     fdrthresh=0;
% end
% pval_vec(~fdrmask) = 1;
%
% figure()
% plot(pval_vec); hold on;
% plot(sem_dat.(sem_dat.data_name{2}).cfg.alt_stt.vis_ovr_all_stt.prob(sig_prb_ind(iSG),:),'--g');
% plot(sem_dat.(sem_dat.data_name{2}).cfg.alt_stt.vis_ovr_all_stt_3.prob(sig_prb_ind(iSG),:),'-.c');
%
% figure();
% plot(Pmtx(ichan,:)); hold on;
% plot(Pmtx(ichan,:),'--r');
