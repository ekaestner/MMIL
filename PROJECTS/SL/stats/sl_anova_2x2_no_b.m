%% Auditory 2-x-2 Anova
vow_ind = {1:4 5:8 9:11};
con_ind = {[1 3 4] 5:8 9:12};

tme_pnt = numel(stt_dat.time{1});

for iBL = 1:numel(vow_ind)
    for iVW = 1:numel(vow_ind{iBL})
        for iCN = 1:numel(con_ind{iBL})
            trl_ind{iBL}{iVW,iCN} = find(stt_dat.cfg.alt_eve.a_vow==vow_ind{iBL}(iVW) & stt_dat.cfg.alt_eve.a_con==con_ind{iBL}(iCN) & stt_dat.cfg.alt_eve.trialinfo~=4);
            vow_idn{iBL}{iVW,iCN} = repmat({num2str(iVW)},1,numel(trl_ind{iBL}{iVW,iCN}));
            con_idn{iBL}{iVW,iCN} = repmat({num2str(iCN)},1,numel(trl_ind{iBL}{iVW,iCN}));
        end
    end
end

ovr_dat_hld = cat(3,stt_dat.trial{:});

aud_vow_anv_stt = zeros(numel(vow_ind),numel(stt_dat.label),tme_pnt);
aud_con_anv_stt = zeros(numel(vow_ind),numel(stt_dat.label),tme_pnt);
aud_int_anv_stt = zeros(numel(vow_ind),numel(stt_dat.label),tme_pnt);

for iCH = 1:numel(stt_dat.label)
    for iBL = 1:numel(vow_ind)
                
        % Create dat
        p_hld = zeros(3,tme_pnt);
        dat_hld = squeeze(ovr_dat_hld(iCH,:,cat(1,trl_ind{iBL}{:})'));
        g_vow   = cat(2,vow_idn{iBL}{:});
        c_vow   = cat(2,con_idn{iBL}{:});
        
        % Run Tests
        for iTM = 1:tme_pnt
            p = anovan(dat_hld(iTM,:),{g_vow,c_vow},'model','interaction','display','off');
            aud_vow_anv_stt(iBL,iCH,iTM) = p(1);
            aud_con_anv_stt(iBL,iCH,iTM) = p(2);
            aud_int_anv_stt(iBL,iCH,iTM) = p(3);
        end
        
    end
end

stt_dat.cfg.alt_stt.aud_vow_blk1.prob = squeeze(aud_vow_anv_stt(1,:,:))'; stt_dat.cfg.alt_stt.aud_vow_blk1.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_vow_blk1.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_vow_blk2.prob = squeeze(aud_vow_anv_stt(2,:,:))'; stt_dat.cfg.alt_stt.aud_vow_blk2.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_vow_blk2.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_vow_blk3.prob = squeeze(aud_vow_anv_stt(3,:,:))'; stt_dat.cfg.alt_stt.aud_vow_blk3.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_vow_blk3.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_con_blk1.prob = squeeze(aud_con_anv_stt(1,:,:))'; stt_dat.cfg.alt_stt.aud_con_blk1.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_con_blk1.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_con_blk2.prob = squeeze(aud_con_anv_stt(2,:,:))'; stt_dat.cfg.alt_stt.aud_con_blk2.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_con_blk2.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_con_blk3.prob = squeeze(aud_con_anv_stt(3,:,:))'; stt_dat.cfg.alt_stt.aud_con_blk3.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_con_blk3.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_int_blk1.prob = squeeze(aud_int_anv_stt(1,:,:))'; stt_dat.cfg.alt_stt.aud_int_blk1.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_int_blk1.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_int_blk2.prob = squeeze(aud_int_anv_stt(2,:,:))'; stt_dat.cfg.alt_stt.aud_int_blk2.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_int_blk2.time = stt_dat.time{1};
stt_dat.cfg.alt_stt.aud_int_blk3.prob = squeeze(aud_int_anv_stt(3,:,:))'; stt_dat.cfg.alt_stt.aud_int_blk3.label = stt_dat.label; stt_dat.cfg.alt_stt.aud_int_blk3.time = stt_dat.time{1};

% Add FDR correction
ttt = fieldnames(stt_dat.cfg.alt_stt);
ttt = ttt(~cellfun(@isempty,regexpi(ttt,'aud.+blk')));
for iFD = 1:numel(ttt)
    stt_dat.cfg.alt_stt.([ttt{iFD} '_fdr']) = stt_dat.cfg.alt_stt.(ttt{iFD});
    for iCH = 1:size(stt_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob,1)
        stt_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:) = fdr_correction(stt_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:),0.05);
    end
end