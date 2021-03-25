function FW_figure7

%% Set-up Places to look
fcfg = [];

fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.tsk     = 'FW';

fcfg.sbj_nme = 'NY226_FW';
fcfg.chn_nme = {'LPT1' 'LPT3' 'LMT2' 'LTO4' 'G37' 'G06' 'G21' 'G12' };
fcfg.eff_lbl = {'Word'};
fcfg.eff_col = {rgb('red')};

fcfg.sub_plt = [3 1];
fcfg.eff_typ = [1 3 4];
fcfg.eff_nme = {'Word' 'Letters' 'False-Font'};
fcfg.eff_xlm = [-0.1 0.550];

%% Central Picture
% Overall
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/']);

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' fcfg.sbj_nme '/' 'surf' '/' 'lh.pial']                  [fcfg.clr_fld '/' 'electrode_location_files' '/' fcfg.sbj_nme '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' fcfg.sbj_nme '/' 'output' '/' fcfg.sbj_nme '_lhs_ecog'] [fcfg.clr_fld '/' 'electrode_location_files' '/' fcfg.sbj_nme '/'  'output' '/' fcfg.sbj_nme '_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.chn_nme;
cfg.all_ele   = fcfg.chn_nme;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = {};

cfg.sep_str         = ',';

cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v3(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = {};
mmil_loc_dot(pcfg)

%% Connectivity Plots
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'plv' '/' 'lines_6_12_maxlne'])

load([fcfg.clr_fld '/' 'sig_chn' '/' 'connectivity' '/' 'plv' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_sig_plv.mat'])

% Smooth
HGP_smooth_msec = 0.050; %0.025;
w               = round(HGP_smooth_msec* (1/(plt(1).tme(2)-plt(1).tme(1))) );
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,round(0.032* (1/(plt(1).tme(2)-plt(1).tme(1))) ));
pcfg.window      = gauss_y/sum(gauss_y);

for iE = 1:numel(plt)
    plt_dat(iE).plv = zeros(size(plt(iE).plv));
    for iC = 1:size(plt_dat(iE).plv,1)
        for iF = 1:size(plt_dat(iE).plv,2)
            plt_dat(iE).plv(iC,iF,:) = conv(squeeze(plt(iE).plv(iC,iF,:)),pcfg.window,'same');
        end
    end
    plt(iE).plv = plt_dat(iE).plv;
end

% Get highest value
plv_hld = [plt(1).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(2).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(3).plv(:,:,1:dsearchn(plt(1).tme',0)) plt(4).plv(:,:,1:dsearchn(plt(1).tme',0))];
plv_hld = plv_hld(:);
plv_max = quantile(plv_hld(:),.999);

% 
scr_sze = get(0,'screensize');
eff_pos = 1;
for iP1 = 1:numel(fcfg.chn_nme)-1
    for iP2 = iP1+1:numel(fcfg.chn_nme)
          
        %% 3D Plot
%         figure('units','pixels','Visible','off','Position',[1 1 min(scr_sze(3:4)) min(scr_sze(3:4))]);
%         
%         for iSP = 1:numel(fcfg.eff_typ)
%             
%             cmb_num = find(strcmpi(plt(fcfg.eff_typ(iSP)).labelcmb,[fcfg.chn_nme{iP1} '--' fcfg.chn_nme{iP2}]));
%             if isempty(cmb_num); cmb_num = find(strcmpi(plt(fcfg.eff_typ(iSP)).labelcmb,[fcfg.chn_nme{iP2} '--' fcfg.chn_nme{iP1}])); end
%             
%             subplot(fcfg.sub_plt(1),fcfg.sub_plt(2),iSP)
%             
%             colormap(hot(200));
%             
%             surf(plt(fcfg.eff_typ(iSP)).tme,plt(fcfg.eff_typ(iSP)).frq,squeeze(double(plt(fcfg.eff_typ(iSP)).plv(cmb_num,:,:))),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
%             
%             set(gca,'View',[0 90]); axis tight;
%             set(gca,'YGrid','off','XGrid','off','clim',[0 0.5]);
% %             title(fcfg.eff_nme{iSP},'FontSize',23)
%             
%             set(gca,'FontSize',1)
%             set(gca,'XColor',[1 1 1]);
%             set(gca,'YColor',[1 1 1]);
%             xlim(fcfg.eff_xlm)
%             
%             zlm = max(get(gca,'zlim'));
%             
%             line([0 0],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('reddish grey'),'LineWidth',7)
%             line([0.2 0.2],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('white'),'LineWidth',7)
%             line([0.4 0.4],[4 24],[zlm+zlm*0.1 zlm+zlm*0.1],'Color',rgb('white'),'LineWidth',7)
%         end
%         
% %         set(gca,'XTickLabels',num2str(str2num(get(gca,'XTickLabels'))*1000));
% %         set(gca,'XColor',[0 0 0]);
% %         ylabel('FREQUENCY (hz)','FontSize',23);
% %         xlabel('TIME (ms)','FontSize',23);
%         
%         set(gcf,'color','w');
%         tightfig();
%         
%         pcfg = [];
%         pcfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'plv' '/' plt(fcfg.eff_typ(iSP)).labelcmb{cmb_num}];
%         pcfg.prn_typ = 'eps';
%         pcfg.jpg = 1;
%         pcfg.cst_col=1;
%         mmil_print_plot(pcfg);
%         close all;
        
        %% 2D Plot
        plt_frq = [6 12];
        
        cmb_num = find(strcmpi(plt(fcfg.eff_typ(1)).labelcmb,[fcfg.chn_nme{iP1} '--' fcfg.chn_nme{iP2}]));
        if isempty(cmb_num); cmb_num = find(strcmpi(plt(fcfg.eff_typ(1)).labelcmb,[fcfg.chn_nme{iP2} '--' fcfg.chn_nme{iP1}])); end
        
        figure('units','pixels','Visible','off','Position',[1 1 min(scr_sze(3:4)) min(scr_sze(3:4))]);
        hold on
        
        line([0 0],[-0.1 1],'Color',rgb('red'),'LineWidth',3)
        
        dat_plt_one = squeeze(double(plt(fcfg.eff_typ(1)).plv(cmb_num,:,:)));
        frq_smt_msc = 0.025;
        win         = round(frq_smt_msc * 500);
        gss_xxx     = -win:1:win;
        gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
        wdw  = gss_yyy/sum(gss_yyy);
        for iF = 1:size(dat_plt_one,1)
            dat_plt_one(iF,:) = conv(dat_plt_one(iF,:),wdw,'same');
        end
        
        dat_plt_two = squeeze(double(plt(fcfg.eff_typ(3)).plv(cmb_num,:,:)));
        frq_smt_msc = 0.025;
        win         = round(frq_smt_msc * 500);
        gss_xxx     = -win:1:win;
        gss_yyy     = normpdf(gss_xxx,0,round(0.016* 500));
        wdw  = gss_yyy/sum(gss_yyy);
        for iF = 1:size(dat_plt_two,1)
            dat_plt_two(iF,:) = conv(dat_plt_two(iF,:),wdw,'same');
        end
        
        plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('bright red'),'LineWidth',15);
        plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('bright red'),'LineWidth',15);
        plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('bright red'),'LineWidth',15);
        plot(plt(fcfg.eff_typ(1)).tme,mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('bright red'),'LineWidth',15);
        stt.one = find(mean(dat_plt_one(dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)) > plv_max);
        %         plot(plt(fcfg.eff_typ(2)).tme,mean(squeeze(double(plt(fcfg.eff_typ(2)).plv(cmb_num,dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(2)).frq',plt_frq(2)),:)))),'Color',rgb('purple'),'LineWidth',5)
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:)),'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0035,'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0070,'Color',rgb('reddish grey'),'LineWidth',15);
%         plot(plt(fcfg.eff_typ(3)).tme,mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(1)).frq',plt_frq(2)),:))+0.0105,'Color',rgb('reddish grey'),'LineWidth',15);
%         stt.thr = find(mean(dat_plt_two(dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(1)):dsearchn(plt(fcfg.eff_typ(3)).frq',plt_frq(2)),:)) > plv_max);
        
        set(gcf,'color','w');
           
        line([0 0],[-0.1 1],'Color',rgb('dark red'),'LineWidth',5)
        line([0.2 0.2],[-0.1 0.6],'Color',rgb('black'),'LineWidth',5)
        line([0.4 0.4],[-0.1 0.6],'Color',rgb('black'),'LineWidth',5)
        line([-0.2 0.7],[0 0],'Color',rgb('black'),'LineWidth',2)
        line([-0.2 0.7],[plv_max plv_max],'Color',rgb('black'),'LineWidth',5,'LineStyle','--')
        
        set(gca,'Xlim',[-0.1 0.550])
        set(gca,'Ylim',[-0.1 0.6])
        set(gca,'FontSize',1)
        set(gca,'XColor',[1 1 1]);
        set(gca,'YColor',[1 1 1]);
        
        ylm_shd = {[0.50 0.55] [] [0.55 0.60]};
        clr_shd = {rgb('red') [] rgb('reddish grey')};
        nme_shd = {'one' '' 'thr'};
        for iST = [1]% 3]
            stats = stt.(nme_shd{iST});
            stats(1:2) = 0; stats(end-1:end) = 0;
            sig_periods = crossing(stats(:),[],0.5); sig_1st = sig_periods(1:2:end); sig_2nd = sig_periods(2:2:end)-1;
            stats(1) = []; stats(end) = [];
            
            for iSG = 1:numel(sig_1st)
                dat_stt_one = repmat(ylm_shd{iST}(1),1,(sig_2nd(iSG)-sig_1st(iSG))+1); 
                dat_stt_two = repmat(ylm_shd{iST}(2),1,(sig_2nd(iSG)-sig_1st(iSG))+1); 
                patch([plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG):sig_2nd(iSG))) fliplr(plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG):sig_2nd(iSG))))],[dat_stt_one dat_stt_two],clr_shd{iST},'EdgeColor','none','FaceAlpha',0.45);
                tme_hld{eff_pos,iST+1} = [sprintf('%.0f',plt(fcfg.eff_typ(iST)).tme(stats(sig_1st(iSG)))*1000) ' - ' sprintf('%.0f',plt(fcfg.eff_typ(iST)).tme(stats(sig_2nd(iSG)))*1000)];
            end
        end
        
        tightfig();
        
        pcfg = [];
        pcfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'plv' '/' 'lines_6_12_maxlne' '/' plt(fcfg.eff_typ(1)).labelcmb{cmb_num} '_line_smt'];
        pcfg.prn_typ = 'eps';
        pcfg.jpg = 0;
        pcfg.cst_col=1;
        mmil_print_plot(pcfg);
        close all;
        
        tme_hld{eff_pos,1} = plt(fcfg.eff_typ(1)).labelcmb{cmb_num};
        
        eff_pos = eff_pos + 1;
        
    end
end
cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'plv' '/' 'plv_timing.csv'],tme_hld);

% Colorbar
pcfg = [];
pcfg.col_map = 'hot(200)';
pcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'plv' '/'];
pcfg.col_bar = [0 0.5];
mmil_color_bar(pcfg)

%% Channel Plots

% Line Plots
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'line_plots']);

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

for iT = 1:numel(fcfg.chn_nme)
    
    cfg = [];
    
    cfg.type      = 'chan';
    cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,fcfg.chn_nme{iT}))};
    cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
    cfg.dat_loc   = 1;
    
    cfg.plt_dim   = [1 1];
    
    cfg.y_lim     = [-6*10^9 10*10^10];
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    
    cfg.alt_eve = 'trialinfo';
    cfg.eve     = [3 6]; %[3 5 6];
    cfg.lnstyle.col_ord = {rgb('bright red') rgb('reddish grey')}; %{rgb('bright red') rgb('purple') rgb('reddish grey')};
    
    cfg.stt_dat = { 'vis_wrd_ffn_msk_01'};
    cfg.stt_col = { { rgb('red') } };
    cfg.stt_lab = 'stt_lab';
    cfg.stt_cmp = { { '90%99' } };
    
    cfg.v_lne       = [0 0.2 0.4];
    cfg.v_lne_wdt    = [3 1 1];
    cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
    
    cfg.x_lim = [-0.1 0.550];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure7' '/' 'line_plots_two' '/' ];
    cfg.prefix     = [fcfg.sbj_nme '_' fcfg.chn_nme{iT} '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
    cfg.print_type = 'eps';
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

end



