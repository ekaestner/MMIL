%% Load Data
% cfg = [];
% cfg.channel   = sll_dat.(sll_dat.data_name{1}).label(iCH);
% sll_dat       = ft_func(@ft_preprocessing,cfg,sll_dat);
% 
% foi =  [2:12,14:2:24,25,30:5:55,70:10:170];    %frequency of interest
% sf  =   [1,1,repmat(2,1,9),repmat(3,1,6),repmat(5,1,7),repmat(10,1,11)]; %specific frequency
% gwidth = ones(size(foi))*pi; %wavelet
% 
% cfg = [];
% cfg.data_name  = 1;
% cfg.method     = 'wavelet';
% cfg.output     = 'fourier';
% cfg.foi        = foi;	      % frequency of interest
% cfg.sf         = sf;        % specific frequency
% cfg.width      = foi./sf;
% cfg.gwidth     = gwidth;
% cfg.toi        = sll_dat.(sll_dat.data_name{1}).time{1}; %cellfun(@(x) sll_dat.(x).time{1},sll_dat.data_name,'uni',0);
% cfg.keeptrials = 'yes';
% sll_dat        = ft_func(@ft_freqanalysis,cfg,sll_dat);
% 
% % Trim Timing
% cfg = [];
% cfg.data_name = 1;
% cfg.latency   = [-0.5 1.5];
% sll_dat       = ft_func(@ft_selectdata,cfg,sll_dat);
% 
% sll_dat       = sll_dat.NY439_SL_Day3_Block1_1_Clin1;
% sll_dat.power = abs(sll_dat.fourierspctrm);
% sll_dat       = rmfield(sll_dat,'fourierspctrm');
% 
% cfg = [];
% cfg.return_events = 0;
% cfg.old_events  = {[1 2] 4};
% cfg.new_events  = {111 112};
% cfg.crt_alt_eve = 'aud_nse';
% sll_dat = ft_redefine_events(cfg,sll_dat);

sbj_num = 1;

subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting work on %s \n',subj)

infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_freq/' subj '_overall_power_data.mat'];
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_freq/initial_plot/NY439_SL/power_plot';

cfg = [];
cfg.load = 'yes';
cfg.file = infile;
sll_dat  = ft_func([],cfg);
sll_dat  = sll_dat.(sll_dat.data_name{1});

%% Initial specifications
eve_cmp     = {[101 102] [111 112] [1 2]};
eve_cmp_lab = {'vis_nse' 'aud_nse' 'trialinfo'};
eve_lab     = {{'Visual' 'Noise'} {'Auditory' 'Noise'} {'Match' 'Mismatch'}};
alt_eve_nme = 'trialinfo';
sve_lab     = 'repair_macro_ejk1st_meso_label';
bse_tme     = [-0.5 0];

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {[1 2] 4};
cfg.new_events  = {111 112};
cfg.crt_alt_eve = 'aud_nse';
sll_dat = ft_redefine_events(cfg,sll_dat);

%% Stats & Plot
% bonferroni
for iCH = 1:numel(sll_dat.label)
    
    figure('Visible','off');
    
    for iCM = 1:numel(eve_cmp)
        
        % bad trials
        bad_trl = find(squeeze(sll_dat.cfg.badtrl(end,iCH,:)));
        
        ev1_inc = setdiff(find(sll_dat.cfg.alt_eve.(eve_cmp_lab{iCM})==eve_cmp{iCM}(1)),bad_trl);
        ev2_inc = setdiff(find(sll_dat.cfg.alt_eve.(eve_cmp_lab{iCM})==eve_cmp{iCM}(2)),bad_trl);
        
        %
        ev1 = squeeze(mean(squeeze(sll_dat.power(ev1_inc,iCH,:,:)),1));
        ev2 = squeeze(mean(squeeze(sll_dat.power(ev2_inc,iCH,:,:)),1));
        
        % z-score
        [~,btm(1)] = min(abs(bse_tme(1) - sll_dat.time));
        [~,btm(2)] = min(abs(bse_tme(2) - sll_dat.time));
        
        ev1_men = nanmean(ev1(:,btm(1):btm(2)),2);
        ev1_std = nanstd(ev1(:,btm(1):btm(2)),[],2);
        
        ev2_men = nanmean(ev2(:,btm(1):btm(2)),2);
        ev2_std = nanstd(ev2(:,btm(1):btm(2)),[],2);
        
        ev1     = bsxfun(@rdivide,bsxfun(@minus,ev1,ev1_men),ev1_std);
        ev2     = bsxfun(@rdivide,bsxfun(@minus,ev2,ev2_men),ev1_std);
        
        % diff
        eve_dff = ev1 - ev2;
        
        % Bonferroni
        ev1(ev1<6 & ev1>-6) = 0;
        ev2(ev2<6 & ev2>-6) = 0;
        eve_dff(eve_dff<6 & eve_dff>-6) = 0;
                
        % plot
        subplot(numel(eve_cmp),3,1+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,ev1,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title(eve_lab{iCM}{1});
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        subplot(numel(eve_cmp),3,2+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,ev2,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title(eve_lab{iCM}{2});
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        subplot(numel(eve_cmp),3,3+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,eve_dff,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title([eve_lab{iCM}{1} ' - ' eve_lab{iCM}{2}]);
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
    end
    
    ppos = get(gcf,'PaperPosition');
    su = get(gcf,'Units');
    pu = get(gcf,'PaperUnits');
    set(gcf,'Units',pu);
    spos = get(gcf,'Position');
    set(gcf,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
    set(gcf,'Units',su)
    print(gcf,sprintf('%s/bonf/%s_%s.png',outpath,subj,sll_dat.cfg.alt_lab.(sve_lab){iCH}),'-dpng','-r300') % - EJK implement chan-by-chan
    close all
end

% fdr
for iCH = 1:numel(sll_dat.label)
    
    figure('Visible','off');
    
    for iCM = 1:numel(eve_cmp)
        
        % bad trials
        bad_trl = find(squeeze(sll_dat.cfg.badtrl(end,iCH,:)));
        
        ev1_inc = setdiff(find(sll_dat.cfg.alt_eve.(eve_cmp_lab{iCM})==eve_cmp{iCM}(1)),bad_trl);
        ev2_inc = setdiff(find(sll_dat.cfg.alt_eve.(eve_cmp_lab{iCM})==eve_cmp{iCM}(2)),bad_trl);
        
        %
        ev1 = squeeze(mean(squeeze(sll_dat.power(ev1_inc,iCH,:,:)),1));
        ev2 = squeeze(mean(squeeze(sll_dat.power(ev2_inc,iCH,:,:)),1));
        
        % z-score
        [~,btm(1)] = min(abs(bse_tme(1) - sll_dat.time));
        [~,btm(2)] = min(abs(bse_tme(2) - sll_dat.time));
        
        ev1_men = nanmean(ev1(:,btm(1):btm(2)),2);
        ev1_std = nanstd(ev1(:,btm(1):btm(2)),[],2);
        
        ev2_men = nanmean(ev2(:,btm(1):btm(2)),2);
        ev2_std = nanstd(ev2(:,btm(1):btm(2)),[],2);
        
        ev1     = bsxfun(@rdivide,bsxfun(@minus,ev1,ev1_men),ev1_std);
        ev2     = bsxfun(@rdivide,bsxfun(@minus,ev2,ev2_men),ev1_std);
        
        % diff
        eve_dff = ev1 - ev2;
               
        % FDR
        ev1_pvl     = arrayfun(@(x) 2 * (1 - normcdf(x,0,1)),ev1);
        ev2_pvl     = arrayfun(@(x) 2 * (1 - normcdf(x,0,1)),ev2);
        eve_dff_pvl = arrayfun(@(x) 2 * (1 - normcdf(x,0,1)),eve_dff);
        
        [~,ev1_fdrmask]     = fdr_cor([ev1_pvl(:)]',0.01);
        [~,ev2_fdrmask]     = fdr_cor([ev2_pvl(:)]',0.01);
        [~,eve_dff_fdrmask] = fdr_cor([eve_dff_pvl(:)]',0.01);
        ev1_fdrmask     = vec2mat(ev1_fdrmask,size(ev1_pvl,1))';
        ev2_fdrmask     = vec2mat(ev2_fdrmask,size(ev2_pvl,1))';
        eve_dff_fdrmask = vec2mat(eve_dff_fdrmask,size(eve_dff_pvl,1))';

        ev1(~ev1_fdrmask)     = 0;
        ev2(~ev1_fdrmask)     = 0;
        eve_dff(~ev1_fdrmask) = 0;
        
        % plot
        subplot(numel(eve_cmp),3,1+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,ev1,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title(eve_lab{iCM}{1});
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        subplot(numel(eve_cmp),3,2+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,ev2,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title(eve_lab{iCM}{2});
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
        subplot(numel(eve_cmp),3,3+(3*(iCM-1)))
        surf(sll_dat.time, sll_dat.freq,eve_dff,'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
        set(gca,'View',[0 90]); axis tight;
        set(gca, 'XColor', [1 1 1],'YColor', [1 1 1],'Color', [1 1 1],'YGrid','off','XGrid','off','clim',[-15 15],'xlim',[sll_dat.time(1) sll_dat.time(end)]);
        title([eve_lab{iCM}{1} ' - ' eve_lab{iCM}{2}]);
        line('XData',[0 0],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        line('XData',[0.45 0.45],'YData',ylim,'ZData',[1 1],'LineWidth',0.4,'Color',[.6 .6 .6]); hold off;
        
    end
    
    ppos = get(gcf,'PaperPosition');
    su = get(gcf,'Units');
    pu = get(gcf,'PaperUnits');
    set(gcf,'Units',pu);
    spos = get(gcf,'Position');
    set(gcf,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
    set(gcf,'Units',su)
    print(gcf,sprintf('%s/fdr/%s_%s.png',outpath,subj,sll_dat.cfg.alt_lab.(sve_lab){iCH}),'-dpng','-r300') % - EJK implement chan-by-chan
    close all
end







