function iSASZ_word_char(fcfg)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.sbj_dat_hld '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

%% Test interrelatedness of variables
eve_loc = { 'vis_lng_med' 'vis_big_med' 'vis_frq_med' 'vis_ngh_med' };
num_loc = { 'vis_lng'     'vis_big'     'vis_frq'     'vis_ngh' };

cmp = nchoosek(1:4,2);
fix = size(cmp,1);

tabulate(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.vis_lng_med)

cfg = [];
cfg.eve_rul     = {'eve_sql' };
cfg.old_events  = { {131 132} };
cfg.use_alt_eve = {'vis_lng_med' };
cfg.fld_nme     = {'vis_big' };
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

figure()
for iC = 1:size(cmp,1)
    
    dat1 = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)});
    dat1 = dat1(~isnan(dat1));
    
    dat2 = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)});
    dat2 = dat2(~isnan(dat2));
    
    [~,p] = corrcoef(dat1,dat2); p = p(1,2);
    subplot(3,3,iC);
    scatter(dat1,dat2);
    title([num2str(p) ' ' num2str(numel(dat1))])
    xlabel(mmil_spec_char(eve_loc{cmp(iC,1)},{'_'}));
    ylabel(mmil_spec_char(eve_loc{cmp(iC,2)},{'_'}));
    
end
tightfig()


eve_loc = { 'aud_lng_med' 'aud_big_med' 'aud_frq_med' 'aud_ngh_med' };
num_loc = { 'aud_lng'     'aud_big'     'aud_frq'     'aud_ngh' };

cmp = nchoosek(1:4,2);
fix = size(cmp,1);

figure()
for iC = 1:size(cmp,1)
    
    dat1 = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)});
    dat1 = dat1(~isnan(dat1));
    
    dat2 = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)});
    dat2 = dat2(~isnan(dat2));
    
    [~,p] = corrcoef(dat1,dat2); p = p(1,2);
    subplot(3,3,iC);
    scatter(dat1,dat2);
    title([num2str(p) ' ' num2str(numel())])
    xlabel(mmil_spec_char(eve_loc{cmp(iC,1)},{'_'}));
    ylabel(mmil_spec_char(eve_loc{cmp(iC,2)},{'_'}));
end
tightfig()

%%
% figure()
% 
% iC = 1;
% subplot(2,2,iC)
% hold on;
% eve     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC}))));
% col_hld = distinguishable_colors(3);

for iCP = 1:numel(eve)
    eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})==eve(iCP);
    val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{iC})(eve_hld);
    scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5, ...
        val_hld,'MarkerFaceColor',col_hld(iCP,:));
    line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
end
xlim([1.5 numel(eve)+1.5])

subplot(2,2,iC+1)
hold on;
eve     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC}))));
col_hld = distinguishable_colors(3);

for iCP = 1:numel(eve)
    eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})==eve(iCP);
    val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{iC+1})(eve_hld);
    scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5, ...
        val_hld,'MarkerFaceColor',col_hld(iCP,:));
    line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
end
xlim([1.5 numel(eve)+1.5])

cfg = [];
cfg.eve_rul     = {'eve_sql' };
cfg.old_events  = { {131 132} };
cfg.use_alt_eve = {'vis_lng_med' };
cfg.fld_nme     = {'vis_big' };
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

iC = 1;
subplot(2,2,iC+2)
hold on;
eve     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC}))));
col_hld = distinguishable_colors(3);

for iCP = 1:numel(eve)
    eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})==eve(iCP);
    val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{iC})(eve_hld);
    scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5, ...
        val_hld,'MarkerFaceColor',col_hld(iCP,:));
    line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
end
xlim([1.5 numel(eve)+1.5])

subplot(2,2,iC+3)
hold on;
eve     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC}))));
col_hld = distinguishable_colors(3);

for iCP = 1:numel(eve)
    eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{iC})==eve(iCP);
    val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{iC+1})(eve_hld);
    scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5, ...
        val_hld,'MarkerFaceColor',col_hld(iCP,:));
    line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
end
xlim([1.5 numel(eve)+1.5])

%% Normalize relations - SZ
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.sbj_dat_hld '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152] [161 161 0 162 162]};
cfg.crt_alt_eve = {'vis_lng_med'       'vis_big_med'       'vis_frq_med'       'vis_ngh_med'};
cfg.fld_nme     = {'vis_lng'           'vis_big'           'vis_frq'           'vis_ngh'};
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

%
mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr']);

eve_loc = { 'vis_lng_med' 'vis_big_med' 'vis_frq_med' 'vis_ngh_med' };
num_loc = { 'vis_lng'     'vis_big'     'vis_frq'     'vis_ngh' };

cmp = nchoosek(1:4,2);
cmp = [cmp ; fliplr(cmp)];
fix = size(cmp,1);

%
for iC = 1:fix
   
%     figure()
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    eve_two     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)}))));
    col_hld = distinguishable_colors(numel(eve_one));

    %
    num_hld{iC,1} = [num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)}];
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,2} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
%     % NaturalPlot - one/two
%     subplot(2,5,1); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,2); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % NaturalPlot - two/one
%     subplot(2,5,4); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,5); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%         
    % Fix Numbers
    cfg = [];
    cfg.eve_rul     = { 'eve_sql' };
    cfg.old_events  = { { eve_one(2) eve_one(3) } };
    cfg.use_alt_eve = { eve_loc{cmp(iC,1)} };
    cfg.fld_nme     = { num_loc{cmp(iC,2)} };
    tst_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
        
%     % FixedPlot - one/two
%     subplot(2,5,6); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,7); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % FixedPlot - two/one
%     subplot(2,5,9); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,10); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     tightfig()
%     
%     print(gcf,[fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr' '/' num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)} '.png'],'-dpng')
%     close all
    
    % Save numbers
    eve_one     = unique(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,3} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
end

%% Normalize relations - SA
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.sbj_dat_hld '/' 'epoch_data' '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152] [161 161 0 162 162]};
cfg.crt_alt_eve = {'aud_lng_med'       'aud_big_med'       'aud_frq_med'       'aud_ngh_med'};
cfg.fld_nme     = {'aud_lng'           'aud_big'           'aud_frq'           'aud_ngh'};
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

%
mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr']);

eve_loc = { 'aud_lng_med' 'aud_big_med' 'aud_frq_med' 'aud_ngh_med' };
num_loc = { 'aud_lng'     'aud_big'     'aud_frq'     'aud_ngh' };

cmp = nchoosek(1:4,2);
cmp = [cmp ; fliplr(cmp)];
fix = size(cmp,1);

%
for iC = 1:fix
   
%     figure()
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    eve_two     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)}))));
    col_hld = distinguishable_colors(numel(eve_one));

    %
    num_hld{iC,1} = [num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)}];
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,2} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
%     % NaturalPlot - one/two
%     subplot(2,5,1); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,2); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % NaturalPlot - two/one
%     subplot(2,5,4); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,5); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%         
    % Fix Numbers
    cfg = [];
    cfg.eve_rul     = { 'eve_sql' };
    cfg.old_events  = { { eve_one(2) eve_one(3) } };
    cfg.use_alt_eve = { eve_loc{cmp(iC,1)} };
    cfg.fld_nme     = { num_loc{cmp(iC,2)} };
    tst_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
        
%     % FixedPlot - one/two
%     subplot(2,5,6); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,7); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % FixedPlot - two/one
%     subplot(2,5,9); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,10); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     tightfig()
%     
%     print(gcf,[fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr' '/' num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)} '.png'],'-dpng')
%     close all
    
    % Save numbers
    eve_one     = unique(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,3} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
end

%% Normalize relations - FW
cfg = [];
cfg.load    = 'yes';
cfg.file    = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data/NY035_FW_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152] [161 161 0 162 162]};
cfg.crt_alt_eve = {'wrd_lng_med'       'wrd_big_med'       'wrd_frq_med'       'wrd_ngh_med'};
cfg.fld_nme     = {'wrd_lng'           'wrd_big'           'wrd_frq'           'wrd_ngh'};
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

%
mkdir([fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr']);

eve_loc = { 'wrd_lng_med' 'wrd_big_med' 'wrd_frq_med' 'wrd_ngh_med' };
num_loc = { 'wrd_lng'     'wrd_big'     'wrd_frq'     'wrd_ngh' };

cmp = nchoosek(1:4,2);
cmp = [cmp ; fliplr(cmp)];
fix = size(cmp,1);

%
for iC = 1:fix
   
%     figure()
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    eve_two     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)}))));
    col_hld = distinguishable_colors(numel(eve_one));

    %
    num_hld{iC,1} = [num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)}];
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,2} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
%     % NaturalPlot - one/two
%     subplot(2,5,1); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,2); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % NaturalPlot - two/one
%     subplot(2,5,4); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,5); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%         
    % Fix Numbers
    cfg = [];
    cfg.eve_rul     = { 'eve_sql' };
    cfg.old_events  = { { eve_one(2) eve_one(3) } };
    cfg.use_alt_eve = { eve_loc{cmp(iC,1)} };
    cfg.fld_nme     = { num_loc{cmp(iC,2)} };
    tst_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
        
%     % FixedPlot - one/two
%     subplot(2,5,6); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     subplot(2,5,7); hold on;
%     for iCP = 1:numel(eve_one)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,1)},{'_'})]);
%     
%     % FixedPlot - two/one
%     subplot(2,5,9); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,2)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,2)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,2)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     subplot(2,5,10); hold on;
%     for iCP = 1:numel(eve_two)
%         eve_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,2)})==eve_two(iCP);                 val_hld = bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(num_loc{cmp(iC,1)})(eve_hld);
%         scatter(ones(numel(val_hld),1)+iCP+rand(numel(val_hld),1)-0.5,val_hld,'MarkerFaceColor',col_hld(iCP,:)); line([iCP+0.75 iCP+1.25],[mean(val_hld) mean(val_hld)],'Color',col_hld(iCP,:))
%     end; xlim([1.5 numel(eve_one)+1.5]); ylabel(mmil_spec_char(num_loc{cmp(iC,1)},{'_'})); title([mmil_spec_char(num_loc{cmp(iC,1)},{'_'}) ' BY ' mmil_spec_char(num_loc{cmp(iC,2)},{'_'})]);
%     
%     tightfig()
%     
%     print(gcf,[fcfg.clr_fld '/' 'events' '/' fcfg.sbj_nme '/' 'lng_chr' '/' num_loc{cmp(iC,1)} '_BY_' num_loc{cmp(iC,2)} '.png'],'-dpng')
%     close all
    
    % Save numbers
    eve_one     = unique(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    num_eve_one = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(2));
    num_eve_two = sum(tst_dat.(tst_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})==eve_one(3));
    
    num_hld{iC,3} = [num2str(num_eve_one) '/' num2str(num_eve_two)];
    
end


end
