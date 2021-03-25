% 
% %% Load Data
% 
% % HGP data
% flt_dat_hld = cat(3,flt_vow_alg_dat.trial{:});
% smt_dat_hld = cat(3,smt_vow_alg_dat.trial{:});
% 
% % Timing
% [~,ddd.tme_000_050] = min(abs(0.000 - flt_vow_alg_dat.time{1})); [~,ddd.tme_000_050(2)] = min(abs(0.050 - flt_vow_alg_dat.time{1}));
% [~,ddd.tme_000_100] = min(abs(0.000 - flt_vow_alg_dat.time{1})); [~,ddd.tme_000_100(2)] = min(abs(0.100 - flt_vow_alg_dat.time{1}));
% [~,ddd.tme_200_250] = min(abs(0.200 - flt_vow_alg_dat.time{1})); [~,ddd.tme_200_250(2)] = min(abs(0.250 - flt_vow_alg_dat.time{1}));
% [~,ddd.tme_200_300] = min(abs(0.200 - flt_vow_alg_dat.time{1})); [~,ddd.tme_200_300(2)] = min(abs(0.300 - flt_vow_alg_dat.time{1}));
% 
% tme_lab{1} = 'tme_000_050';
% tme_lab{2} = 'tme_000_100';
% tme_lab{3} = 'tme_200_250';
% tme_lab{4} = 'tme_200_300';
% 
% % Trial Inds
% vow_ind = {101:104 105:108 109:111};
% con_ind = {1:4 5:8 9:12};
% 
% vow_idn = {{'IH' 'UH' 'OH' 'AY'} {'EE' 'OO' 'AH' 'UR'} {'EH' 'AW' 'OY'}};
% con_idn = {{'N'  'B'  'S'  'W'}  {'M'  'T'  'H'  'F'}  {'L'  'Y'  'R' 'K'}};
% 
% for iBL = 1:numel(vow_ind)
%     for iVW = 1:numel(vow_ind{iBL})
%         for iCN = 1:numel(con_ind{iBL})
%             trl_ind{iBL}{iVW,iCN} = find(flt_vow_alg_dat.cfg.alt_eve.a_vow==vow_ind{iBL}(iVW) & flt_vow_alg_dat.cfg.alt_eve.a_con==con_ind{iBL}(iCN) & flt_vow_alg_dat.cfg.alt_eve.trialinfo~=4);
%             vow_idn{iBL}{iVW,iCN} = repmat({num2str(iVW)},1,numel(trl_ind{iBL}{iVW,iCN}));
%             con_idn{iBL}{iVW,iCN} = repmat({num2str(iCN)},1,numel(trl_ind{iBL}{iVW,iCN}));
%         end
%     end
% end
% 
% %% Make stat table
% clear p1 p2 p3 p4
% for iC = 1:size(flt_dat_hld,1)
%     for iBL = 1:numel(vow_ind)
%         
%         dat_hld = squeeze(flt_dat_hld(iC,:,cat(1,trl_ind{iBL}{:})'));
%         g_vow   = cat(2,vow_idn{iBL}{:});
%         c_vow   = cat(2,con_idn{iBL}{:});
%         
%         p = anovan(mean(dat_hld(ddd.tme_000_050(1):ddd.tme_000_050(2),:)),{g_vow,c_vow},'model','interaction','display','off'); p1(iC,iBL) = p(3); clear p;
%         p = anovan(mean(dat_hld(ddd.tme_000_100(1):ddd.tme_000_100(2),:)),{g_vow,c_vow},'model','interaction','display','off'); p2(iC,iBL) = p(3); clear p;
%         p = anovan(mean(dat_hld(ddd.tme_200_250(1):ddd.tme_200_250(2),:)),{g_vow,c_vow},'model','interaction','display','off'); p3(iC,iBL) = p(3); clear p;
%         p = anovan(mean(dat_hld(ddd.tme_200_300(1):ddd.tme_200_300(2),:)),{g_vow,c_vow},'model','interaction','display','off'); p4(iC,iBL) = p(3); clear p;
%         
%     end
% end
% 
% pvl_tab = cat(3,p1,p2,p3,p4);
% 
% [ind_chn{1},ind_tme{1}] = find(squeeze(pvl_tab(:,1,:))<0.05);
% [ind_chn{2},ind_tme{2}] = find(squeeze(pvl_tab(:,2,:))<0.05);
% [ind_chn{3},ind_tme{3}] = find(squeeze(pvl_tab(:,3,:))<0.05);

%% Make Figures
for iBL = 1:numel(ind_chn)
         
    for iCH = 1:size(flt_dat_hld,1)
   
        if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/InteractionPlots/' flt_vow_alg_dat.label{iCH}],'dir'); mkdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/InteractionPlots/' flt_vow_alg_dat.label{iCH}]); end
   
        if any(ind_chn{iBL}==iCH)
            
            tme_typ = ind_tme{iBL}(ind_chn{iBL}==iCH);
            
            for iTM = 1:numel(tme_typ)
                
                figure('Visible','off'); hold on;
                for iV = 1:size(trl_ind{iBL},1);
                    for iC = 1:size(trl_ind{iBL},2);
                        dat_for_plt(iV,iC) = nanmean(nanmean(squeeze(flt_dat_hld(iCH,ddd.(tme_lab{tme_typ(iTM)})(1):ddd.(tme_lab{tme_typ(iTM)})(2),trl_ind{iBL}{iV,iC})),2),1);
                        plot(iC,dat_for_plt(iV,iC),'o','MarkerFaceColor',cmb_col(int_num(con_ind{iBL}(iC),vow_ind{iBL}(iV)-100),:),'MarkerSize',15)
                    end
                    lne(iV) = plot(dat_for_plt(iV,:),'Color',vow_plt_col(vow_ind{iBL}(iV)-100,:));
                end
                xlim([0.75 4.25])
                set(gca,'XTick',[1 2 3 4])
                set(gca,'XTickLabel',con_idn{iBL})
                legend(lne,vow_idn{iBL})
                
                print(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/InteractionPlots/' flt_vow_alg_dat.label{iCH} '/' flt_vow_alg_dat.label{iCH} '_' 'Block' num2str(iBL) '_' tme_lab{tme_typ(iTM)} '.png'],'-dpng','-r300');
                close all
                
            end
                       
        end
        
        % Copy over waves
%         copyfile(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/2x2FltAlg' '/' 'Filtered_Aud_2x2_Phoneme_average_chan_' flt_vow_alg_dat.label{iCH} '.png'],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/InteractionPlots/' flt_vow_alg_dat.label{iCH}]);
        
        % Copy over Individual Phonemes
%         copyfile(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/ChanFigAlgFlt' '/' flt_vow_alg_dat.label{iCH}],['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot/NY439_SL/TimingCheck/InteractionPlots/' '/' flt_vow_alg_dat.label{iCH} '/' 'phonemes']);
        
    end
end



%% Make Figures for significant periods

%% Biphoneme Probability 


