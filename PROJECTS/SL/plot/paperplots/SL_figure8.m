function SL_figure8

fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';

fcfg.sbj_nme = {   'NY439_SL'               'NY591_SL' };                % {   'NY439_SL'               'NY591_SL' };
fcfg.chn_nme = { { 'LPTr05'  'G20'  'G29' 'GR33' 'GR35' 'G44' 'GR53' 'GR60'} { 'PT01' 'G20' 'G54' } }; % { { 'LPTr15' 'G37' 'G46' } { 'G20' 'G54' 'O03' 'O04' } };
fcfg.loc_nme = { { '' '' '' }             { '' '' '' } };                % { { '' '' '' }             { '' '' '' } };

% 'G21' 'G28' 'GR25' 'GR26' 'GR34' 'GR36' 'GR37' 'GR44' 'GR45' 'GR56' 'G43'  
% 'LPTr05' 'G13' 'G20' 'G21' 'G28' 'G29' 'G36' 'GR19' 'GR21' 'GR24' 'GR25' 'GR26' 'GR33' 'GR34' 'GR35' 'GR36' 'GR37' 'GR44' 'GR45' 'G22' 'G46' 'G55' 'G56' 'GR56' 'LO08' 'G43' 'G44' 'GR52' 'GR53' 'GR60'

fcfg.eff_lbl = {'Word'};
fcfg.eff_col = {rgb('red')};

fcfg.sub_plt = [3 1];
fcfg.eff_typ = [1 3 4];
fcfg.eff_nme = {'Match' 'False-Font' 'Noise-Vocoding'};
fcfg.eff_xlm = [-0.1 1.250];

for iS = 1:numel(fcfg.sbj_nme)
    
    %%
    mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv'])
    mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' '3d'])
    mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' 'line'])
    mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' 'lineHGP'])
    
    load([fcfg.clr_fld '/' 'sig_chn' '/' 'connectivity' '/' 'plv' '/' fcfg.sbj_nme{iS} '/' fcfg.sbj_nme{iS} '_sig_plv.mat'])
    
    % Get highest value
    plv_hld = [plt(1).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(2).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(3).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(4).plv(:,:,1:dsearchn(plt(1).tme',0))];
    plv_hld = plv_hld(:);
    plv_max = quantile(plv_hld(:),.999);
    
    % Smooth
    HGP_smooth_msec = 0.050; %0.025;
    w               = round(HGP_smooth_msec* (1/(plt(1).tme(2)-plt(1).tme(1))) );
    gauss_x         = -w:1:w;
    gauss_y         = normpdf(gauss_x,0,round(0.032* (1/(plt(1).tme(2)-plt(1).tme(1))) ));
    pcfg.window      = gauss_y/sum(gauss_y);
    
    for iE = 1:numel(plt)
        plt(iE).plv(plt(iE).plv<plv_max) = 0;
        plt_dat(iE).plv = zeros(size(plt(iE).plv));
        for iC = 1:size(plt_dat(iE).plv,1)
            for iF = 1:size(plt_dat(iE).plv,2)
                plt_dat(iE).plv(iC,iF,:) = conv(squeeze(plt(iE).plv(iC,iF,:)),pcfg.window,'same');
            end
        end
        plt(iE).plv = plt_dat(iE).plv;
    end
    
    %
    scr_sze = get(0,'screensize');
    eff_pos = 1;
    for iP1 = 1:numel(fcfg.chn_nme{iS})-1
        
        for iP2 = iP1+1:numel(fcfg.chn_nme{iS})
            
            %% 3D Plot
            figure('units','pixels','Visible','off','Position',[1 1 min(scr_sze(3:4)) min(scr_sze(3:4))]);
            
            for iSP = 1:numel(fcfg.eff_typ)
                
                cmb_num = find(strcmpi(plt(fcfg.eff_typ(iSP)).labelcmb,[fcfg.chn_nme{iS}{iP1} '--' fcfg.chn_nme{iS}{iP2}]));
                if isempty(cmb_num); cmb_num = find(strcmpi(plt(fcfg.eff_typ(iSP)).labelcmb,[fcfg.chn_nme{iS}{iP2} '--' fcfg.chn_nme{iS}{iP1}])); end
                
                subplot(fcfg.sub_plt(1),fcfg.sub_plt(2),iSP)
                
                colormap(hot(200));
                
                plt_clr = squeeze(double(plt(fcfg.eff_typ(iSP)).plv(cmb_num,:,:)));
                
                surf(plt(fcfg.eff_typ(iSP)).tme,plt(fcfg.eff_typ(iSP)).frq,plt_clr,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
                
                set(gca,'View',[0 90]); axis tight;
                set(gca,'YGrid','off','XGrid','off','clim',[0 0.5]);
                %             title(fcfg.eff_nme{iSP},'FontSize',23)
                
                set(gca,'FontSize',1)
                set(gca,'XColor',[1 1 1]);
                set(gca,'YColor',[1 1 1]);
                xlim(fcfg.eff_xlm)
                
                zlm = max(get(gca,'zlim'));
                
                line([0.000 0.000],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('reddish grey'),'LineWidth',7)
                line([0.450 0.450],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('bluish grey'),'LineWidth',7)
                line([0.900 0.900],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('grey'),'LineWidth',7)
            end
            
            %         set(gca,'XTickLabels',num2str(str2num(get(gca,'XTickLabels'))*1000));
            %         set(gca,'XColor',[0 0 0]);
            %         ylabel('FREQUENCY (hz)','FontSize',23);
            %         xlabel('TIME (ms)','FontSize',23);
            
            set(gcf,'color','w');
            tightfig();
            
            pcfg = [];
            pcfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' '3d' '/' fcfg.sbj_nme{iS} '_' plt(fcfg.eff_typ(iSP)).labelcmb{cmb_num}];
            pcfg.prn_typ = 'eps';
            pcfg.jpg = 1;
            pcfg.cst_col=1;
            mmil_print_plot(pcfg);
            
            pcfg.prn_typ = 'png';
            mmil_print_plot(pcfg);
            close all;
            
        end
    end
    
    %% Line Plot
    tot_sbj = fcfg.sbj_nme{iS};
    
    % Parameters
    alt_eve = { 'lng_tot_nse' };
    eve     = { [ 601                 603                 604 ] };
    col_ord = { { rgb('light purple') rgb('reddish grey') rgb('bluish grey') } };
    
    stt_dat = { { 'vis_nse_stt_msk_anv'  'aud_nse_stt_msk_anv'} };
    stt_col = { { ft_stt_col(rgb('red')) ft_stt_col(rgb('bright blue')) }  };
    stt_cmp = { { '0%4'                  '4%8' } };
    
    % Put together plot
    % Load
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.dat_fld '/' tot_sbj '_overall_data.mat'];
    bcc_dat  = ft_func([],cfg);
    
    for iP1 = 1:numel(fcfg.chn_nme{iS})
        
        % Plot
        cfg = [];
        
        cfg.type      = 'chan';
        cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,fcfg.chn_nme{iS}{iP1}))};
        cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
        cfg.dat_loc   = 1;
        
        cfg.plt_dim   = [1 1];
        
        cfg.y_lim     = 'maxmin'; %ylm{1};
        
        cfg.lgd       = 0;
        cfg.std_err   = 1;
        
        cfg.alt_eve = alt_eve{1};
        cfg.eve     = eve{1};
        cfg.lnstyle.col_ord = col_ord{1};
        
        cfg.stt_dat = stt_dat(1);
        cfg.stt_col = stt_col(1);
        cfg.stt_lab = 'stt_lab';
        cfg.stt_cmp = stt_cmp(1);
        
        cfg.v_lne       = {[0         0.450       0.900]};
        cfg.v_lne_wdt   = {[10        10          2.5]};
        cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};
        
        cfg.x_lim = [-0.1 1.2];
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'png';
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' 'lineHGP'];
        cfg.prefix     = [fcfg.sbj_nme{iS} '_' fcfg.chn_nme{iS}{iP1}];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
        %% Iterate
        tme_hld{eff_pos,1} = plt(fcfg.eff_typ(1)).labelcmb{cmb_num};
        
        eff_pos = eff_pos + 1;
        
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% 2D Plot
%         plt_frq = [6 12];
%
%         cmb_num = find(strcmpi(plt(fcfg.eff_typ(1)).labelcmb,[fcfg.chn_nme{iP1} '--' fcfg.chn_nme{iP2}]));
%         if isempty(cmb_num); cmb_num = find(strcmpi(plt(fcfg.eff_typ(1)).labelcmb,[fcfg.chn_nme{iP2} '--' fcfg.chn_nme{iP1}])); end
%
%         figure('units','pixels','Visible','off','Position',[1 1 min(scr_sze(3:4)) min(scr_sze(3:4))]);
%         hold on
%
%         line([0 0],[-0.1 1],'Color',rgb('red'),'LineWidth',3)
%
%         dat_plt_one = squeeze(double(plt(fcfg.eff_typ(1)).plv(cmb_num,:,:)));
%         frq_smt_msc = 0.025;
%         win         = round(frq_smt_msc * 500);
%         gss_xxx     = -win:1:win;
%         gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
%         wdw  = gss_yyy/sum(gss_yyy);
%         for iF = 1:size(dat_plt_one,1)
%             dat_plt_one(iF,:) = conv(dat_plt_one(iF,:),wdw,'same');
%         end
%
%         dat_plt_two = squeeze(double(plt(fcfg.eff_typ(2)).plv(cmb_num,:,:)));
%         frq_smt_msc = 0.025;
%         win         = round(frq_smt_msc * 500);
%         gss_xxx     = -win:1:win;
%         gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
%         wdw  = gss_yyy/sum(gss_yyy);
%         for iF = 1:size(dat_plt_two,1)
%             dat_plt_two(iF,:) = conv(dat_plt_two(iF,:),wdw,'same');
%         end
%
%         dat_plt_thr = squeeze(double(plt(fcfg.eff_typ(3)).plv(cmb_num,:,:)));
%         frq_smt_msc = 0.025;
%         win         = round(frq_smt_msc * 500);
%         gss_xxx     = -win:1:win;
%         gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
%         wdw  = gss_yyy/sum(gss_yyy);
%         for iF = 1:size(dat_plt_thr,1)
%             dat_plt_thr(iF,:) = conv(dat_plt_thr(iF,:),wdw,'same');
%         end
%
%         dat_plt_for = squeeze(double(plt(fcfg.eff_typ(4)).plv(cmb_num,:,:)));
%         frq_smt_msc = 0.025;
%         win         = round(frq_smt_msc * 500);
%         gss_xxx     = -win:1:win;
%         gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
%         wdw  = gss_yyy/sum(gss_yyy);
%         for iF = 1:size(dat_plt_for,1)
%             dat_plt_for(iF,:) = conv(dat_plt_for(iF,:),wdw,'same');
%         end
%
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('green'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('green'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('green'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('green'),'LineWidth',15);
%         stt.one = find(mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)) > plv_max);
%
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('yellow'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('yellow'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('yellow'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('yellow'),'LineWidth',15);
%         stt.two = find(mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)) > plv_max);
%
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_thr(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_thr(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_thr(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_thr(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('reddish grey'),'LineWidth',15);
%         stt.thr = find(mean(dat_plt_thr(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(2)),:)) > plv_max);
%
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_for(dsearchn(plt(fcfg.eff_typ(4)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('bluish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_for(dsearchn(plt(fcfg.eff_typ(4)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('bluish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_for(dsearchn(plt(fcfg.eff_typ(4)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('bluish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_for(dsearchn(plt(fcfg.eff_typ(4)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('bluish grey'),'LineWidth',15);
%         stt.for = find(mean(dat_plt_for(dsearchn(plt(fcfg.eff_typ(4)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)) > plv_max);
%
%         set(gcf,'color','w');
%
%         line([0.000 0.000],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('reddish grey'),'LineWidth',7)
%         line([0.450 0.450],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('bluish grey'),'LineWidth',7)
%         line([0.900 0.900],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('grey'),'LineWidth',7)
%
%         line([0.000 0.000],[-0.1 1],'Color',rgb('dark red'),'LineWidth',5)
%         line([0.450 0.450],[-0.1 1],'Color',rgb('dark blue'),'LineWidth',5)
%         line([0.900 0.900],[-0.1 1],'Color',rgb('black'),'LineWidth',5)
%
%         line([-0.2 1.250],[0 0],'Color',rgb('black'),'LineWidth',2)
%         line([-0.2 1.250],[plv_max plv_max],'Color',rgb('black'),'LineWidth',5,'LineStyle','--')
%
%         set(gca,'Xlim',[-0.1 1.250])
%         set(gca,'Ylim',[-0.1 0.8])
%         set(gca,'FontSize',1)
%         set(gca,'XColor',[1 1 1]);
%         set(gca,'YColor',[1 1 1]);
%
%         ylm_shd = {[0.60 0.65] [0.65 0.70] [0.70 0.75] [0.75 0.80]};
%         clr_shd = {rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')};
%         nme_shd = {'one' 'two' 'thr' 'for'};
%         for iST = [1 2 3 4]
%             stats = stt.(nme_shd{iST});
%             stats(1:2) = 0; stats(end-1:end) = 0;
%             sig_periods = crossing(stats(:),[],0.5); sig_1st = sig_periods(1:2:end); sig_2nd = sig_periods(2:2:end)-1;
%             stats(1) = []; stats(end) = [];
%
%             for iSG = 1:numel(sig_1st)
%                 dat_stt_one = repmat(ylm_shd{iST}(1),1,(sig_2nd(iSG)-sig_1st(iSG))+1);
%                 dat_stt_two = repmat(ylm_shd{iST}(2),1,(sig_2nd(iSG)-sig_1st(iSG))+1);
%                 patch([plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG):sig_2nd(iSG))) fliplr(plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG):sig_2nd(iSG))))],[dat_stt_one dat_stt_two],clr_shd{iST},'EdgeColor','none','FaceAlpha',0.45);
%                 tme_hld{eff_pos,iST+1} = [sprintf('%.0f',plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG)))*1000) ' - ' sprintf('%.0f',plt(fcfg.eff_typ(iST)).tme(stats(sig_2nd(iSG)))*1000)];
%             end
%         end
%
%         tightfig();
%
%         pcfg = [];
%         pcfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure8' '/' 'plv' '/' 'line' '/' plt(fcfg.eff_typ(1)).labelcmb{cmb_num} '_line_smt'];
%         pcfg.prn_typ = 'eps';
%         pcfg.jpg = 0;
%         pcfg.cst_col=1;
%         mmil_print_plot(pcfg);
%
%                 pcfg.prn_typ = 'png';
%                 mmil_print_plot(pcfg);
%         close all;