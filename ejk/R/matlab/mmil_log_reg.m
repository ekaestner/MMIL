
function mmil_log_reg(fcfg)

if ~isfield(fcfg,'acc_cut_off'); fcfg.acc_cut_off = .50; end

%%
clear trg_dta fol_dta rmv_mtx_trg rmv_mtx_fol

%% Put together
trg_grp = cell(0);
fol_grp = cell(0);

for iS = 2:size(fcfg.sbj_grp,1)
    if     strcmpi(fcfg.sbj_grp{iS,2},fcfg.lbl_ord{1})
        trg_grp{end+1,1} = fcfg.sbj_grp{iS,1};
    elseif strcmpi(fcfg.sbj_grp{iS,2},fcfg.lbl_ord{2})
        fol_grp{end+1,1} = fcfg.sbj_grp{iS,1};
    end
end

%%
for iM = 1:numel(fcfg.mdl)
        
    trg_dta{iM} = nan(size(trg_grp,1),numel(fcfg.mdl{iM}));
    fol_dta{iM} = nan(size(fol_grp,1),numel(fcfg.mdl{iM}));
    
    for iC = 1:numel(fcfg.mdl{iM})
        if ischar(fcfg.dta{1,strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})})
            chr_lvl{iC} = unique(fcfg.dta(~cellfun(@isempty,fcfg.dta(:,strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC}))),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})));
        else            
            chr_lvl{iC} = '';
        end
    end
    
    % Target Data
    for iS = 1:size(trg_dta{iM},1)
        for iC = 1:numel(fcfg.mdl{iM})
            if ischar(fcfg.dta{strcmpi(fcfg.dta(:,1),trg_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})})
                try trg_dta{iM}(iS,iC) = find(strcmpi(chr_lvl{iC},fcfg.dta{strcmpi(fcfg.dta(:,1),trg_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})})); catch; end
            else
                try trg_dta{iM}(iS,iC) = fcfg.dta{strcmpi(fcfg.dta(:,1),trg_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})}; catch; end
            end
        end
    end
    
    % Foil Data
    for iS = 1:size(fol_dta{iM},1)
        for iC = 1:numel(fcfg.mdl{iM})
            if ischar(fcfg.dta{strcmpi(fcfg.dta(:,1),fol_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})})
                try fol_dta{iM}(iS,iC) = find(strcmpi(chr_lvl{iC},fcfg.dta{strcmpi(fcfg.dta(:,1),fol_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})})); catch; end
            else
                try fol_dta{iM}(iS,iC) = fcfg.dta{strcmpi(fcfg.dta(:,1),fol_grp{iS}),strcmpi(fcfg.dta_lbl(1,:),fcfg.mdl{iM}{iC})}; catch; end
            end
        end
    end
        
end

%% Systematize
if isfield(fcfg,'nrm_grp') && fcfg.nrm_grp==1
    for iM = 1:numel(fcfg.mdl)
        rmv_mtx_trg(:,iM) = sum(isnan(trg_dta{iM}),2);
        rmv_mtx_fol(:,iM) = sum(isnan(fol_dta{iM}),2);
    end
    
    rmv_mtx_trg = sum(rmv_mtx_trg,2);
    rmv_mtx_fol = sum(rmv_mtx_fol,2);
    
    for iM = 1:numel(fcfg.mdl)
        trg_dta{iM}(find(rmv_mtx_trg),:) = [];
        fol_dta{iM}(find(rmv_mtx_fol),:) = [];
    end
    
    trg_grp(find(rmv_mtx_trg),:) = [];
    fol_grp(find(rmv_mtx_fol),:) = [];
else
    for iM = 1:numel(fcfg.mdl)
        rmv_mtx_trg(:,iM) = sum(isnan(trg_dta{iM}),2);
        rmv_mtx_fol(:,iM) = sum(isnan(fol_dta{iM}),2);
                        
        trg_dta{iM}(find(rmv_mtx_trg(:,iM)),:) = [];
        fol_dta{iM}(find(rmv_mtx_fol(:,iM)),:) = [];

    end
end

%% Calculte Model
xax = cell(0); yax = cell(0); auc = cell(0); opt = cell(0);
for iM = 1:numel(fcfg.mdl)
    ydt{iM} = [ones(size(trg_dta{iM},1),1) ; ones(size(fol_dta{iM},1),1)*2 ];
    xdt{iM} = [trg_dta{iM}                 ; fol_dta{iM}];
    
    [bbb{iM},~,stt{iM}] = mnrfit(xdt{iM},ydt{iM});
    pih{iM} = mnrval(bbb{iM},xdt{iM});
    
    [xax{iM},yax{iM},~,auc{iM},opt{iM}] = perfcurve(ydt{iM},pih{iM}(:,1),1,'NBoot',10000);
    
    acc(iM) = ( sum((pih{iM}(:,1)>fcfg.acc_cut_off) & ydt{iM}==1) + sum((pih{iM}(:,1)<fcfg.acc_cut_off) & ydt{iM}==2)) / numel(ydt{iM});
    
end

%% Plot
for iP = 1:numel(fcfg.mdl_cmp_plt)

    if ~exist([ fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP}]); mkdir([ fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP}]); end
   
    % ROC Curves
    figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
    hold on;
    
    for iSP = 1:numel(fcfg.mdl_cmp_plt{iP})
        plot(xax{fcfg.mdl_cmp_plt{iP}(iSP)}(:,1),yax{fcfg.mdl_cmp_plt{iP}(iSP)}(:,1),'Color',fcfg.mld_cmp_col{iP}{iSP},'LineWidth',5.5)
    end
    
    set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
    set(gca,'LineWidth',4)
    set(gca,'FontSize',26)
    [leg_obj,leg_lne] = legend(cellfun(@(x) mmil_spec_char(x,{'_'}),fcfg.mdl_nme(fcfg.mdl_cmp_plt{iP}),'uni',0),'Position',[0.5 0.15 0.4 0.3],'FontSize',20,'Box','off');
    leg_lne = findobj(leg_lne,'type','line');
    set(leg_lne,'LineWidth',8)
    
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'ROC.png'],'-dpng')
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    [leg_obj,leg_lne] = legend(repmat({''},1,numel(fcfg.mdl_cmp_plt{iP})),'Position',[0.6 0.15 0.1 0.3],'FontSize',20,'Box','off');
    leg_lne = findobj(leg_lne,'type','line');
    set(leg_lne,'LineWidth',8)
    
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'ROC_FIGURE.png'],'-dpng')
    close all
    
    % AUC Bar Plot
    figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
    hold on;
    
    for iSP = 1:numel(fcfg.mdl_cmp_plt{iP})
        bar_hld(iC) = barh(numel(fcfg.mdl_cmp_plt{iP}) - (iSP-1),auc{fcfg.mdl_cmp_plt{iP}(iSP)}(1));
            bar_hld(iC).FaceColor = fcfg.mld_cmp_col{iP}{iSP};
            bar_hld(iC).EdgeColor = rgb('black');
            bar_hld(iC).LineWidth = 3.5;
    end
    set(gca,'LineWidth',4)
    set(gca,'FontSize',26)
    set(gca,'Xlim',[0.5 1])
    set(gca,'YTick',1:numel(fcfg.mdl_cmp_plt{iP}))
    set(gca,'YTickLabel',fliplr(cellfun(@(x) mmil_spec_char(x,{'_'}),fcfg.mdl_nme(fcfg.mdl_cmp_plt{iP}),'uni',0)))
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'BarPlot.png'],'-dpng')
    
    
    figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
    hold on;
    
    for iSP = 1:numel(fcfg.mdl_cmp_plt{iP})
        bar_hld(iC) = barh(numel(fcfg.mdl_cmp_plt{iP}) - (iSP-1),auc{fcfg.mdl_cmp_plt{iP}(iSP)}(1));
            bar_hld(iC).FaceColor = fcfg.mld_cmp_col{iP}{iSP};
            bar_hld(iC).EdgeColor = rgb('black');
            bar_hld(iC).LineWidth = 3.5;
    end
    set(gca,'LineWidth',4)
    set(gca,'FontSize',26)
    set(gca,'Xlim',[0.5 1])
    set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9 1.0])
    set(gca,'YTick',1:numel(fcfg.mdl_cmp_plt{iP}))
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'BarPlot_FIGURE.png'],'-dpng')
    close all
    
    % Accuracy Plot
    figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
    hold on;
    
    for iSP = 1:numel(fcfg.mdl_cmp_plt{iP})
        bar_hld(iC) = bar(iSP,acc(fcfg.mdl_cmp_plt{iP}(iSP)));
            bar_hld(iC).FaceColor = fcfg.mld_cmp_col{iP}{iSP};
            bar_hld(iC).EdgeColor = rgb('black');
            bar_hld(iC).LineWidth = 3.5;
    end
    set(gca,'LineWidth',4)
    set(gca,'FontSize',26)
    set(gca,'Ylim',[0.5 1])
    set(gca,'YTick',[0.5 0.6 0.7 0.8 0.9 1])
    set(gca,'XTick',1:numel(fcfg.mdl_cmp_plt{iP}))
    set(gca,'XTickLabel',cellfun(@(x) mmil_spec_char(x,{'_'}),fcfg.mdl_nme(fcfg.mdl_cmp_plt{iP}),'uni',0))
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'ACC' '_' 'BarPlot.png'],'-dpng')
    close all
        
    figure('units','pixels','outerposition',[0 0 1080 1080],'Visible','off');
    hold on;
    
    for iSP = 1:numel(fcfg.mdl_cmp_plt{iP})
        bar_hld(iC) = bar(iSP,acc(fcfg.mdl_cmp_plt{iP}(iSP)));
            bar_hld(iC).FaceColor = fcfg.mld_cmp_col{iP}{iSP};
            bar_hld(iC).EdgeColor = rgb('black');
            bar_hld(iC).LineWidth = 3.5;
    end
    set(gca,'LineWidth',4)
    set(gca,'FontSize',26)
    set(gca,'Ylim',[0.5 1])
    set(gca,'YTick',[0.5 0.6 0.7 0.8 0.9 1])
    set(gca,'XTick',1:numel(fcfg.mdl_cmp_plt{iP}))
    
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    tightfig();
    set(gcf,'units','pixels','outerposition',[0 0 1080 1080]);
    
    print([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'ACC' '_' 'BarPlot_FIGURE.png'],'-dpng')
    close all
    
    % Save AUC & Accuracy
    out_txt = [cellfun(@(x) mmil_spec_char(x,{'_'}),fcfg.mdl_nme(fcfg.mdl_cmp_plt{iP}),'uni',0)' ...
               num2cell(roundsd(cat(1,auc{fcfg.mdl_cmp_plt{iP}}),3)) ...
               repmat({''},numel(fcfg.mdl_cmp_plt{iP}),1) ...
               num2cell(roundsd(acc(fcfg.mdl_cmp_plt{iP})',3)) ...
               repmat({''},numel(fcfg.mdl_cmp_plt{iP}),1) ...
               num2cell(cellfun(@(x) size(x,2),trg_dta(fcfg.mdl_cmp_plt{iP}))') ...
               num2cell(cellfun(@(x) size(x,1),trg_dta(fcfg.mdl_cmp_plt{iP}))') ...
               num2cell(cellfun(@(x,y) size(x,1)+size(y,1),trg_dta(fcfg.mdl_cmp_plt{iP}),fol_dta(fcfg.mdl_cmp_plt{iP}))') ];
    cell2csv([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Plot' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'AUC.csv'],out_txt)
    
    % Save Stats
    for iSV = 1:numel(fcfg.mdl_cmp_plt{iP})
        
        stt_out = num2cell([stt{fcfg.mdl_cmp_plt{iP}(iSV)}.beta stt{fcfg.mdl_cmp_plt{iP}(iSV)}.p ]);
        cell2csv([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'Stats' '_' num2str(iP) '_' fcfg.mdl_cmp_nme{iP} '_' 'ModelNum_' num2str(fcfg.mdl_cmp_plt{iP}(iSV)) '.csv'],stt_out)
        
    end
    
end

%% Save Data
var_fct = {'one' 'two' 'three' 'four' 'five' 'six' 'seven'};
for iP = 1:numel(fcfg.mdl_cmp_plt)
    for iM = 1:numel(fcfg.mdl_cmp_plt{iP})
       
    nme_fct = { 'handedness' 'MTS' 'side' };    
        
    raw_dta = [  [{''} {''} fcfg.mdl{fcfg.mdl_cmp_plt{iP}(iM)}] ; ...
                 [ [trg_grp ; fol_grp ] ...
                   num2cell(ydt{fcfg.mdl_cmp_plt{iP}(iM)}) ...
                   num2cell(xdt{fcfg.mdl_cmp_plt{iP}(iM)}) ] ];
    
    col_fct = find(ismember(raw_dta(1,:),nme_fct));
    for iC = 1:numel(col_fct) 
        lvl_hld = unique(cell2mat(raw_dta(2:end,col_fct(iC))));
        for iL = 1:numel(lvl_hld)
            rep_ind{iL} = find(cell2mat(cellfun(@(x) x==lvl_hld(iL),raw_dta(2:end,col_fct(iC)),'uni',0)))+1;
        end
        for iL = 1:numel(lvl_hld)    
            raw_dta(rep_ind{iL},col_fct(iC)) = var_fct(iL);
        end  
        clear rep_ind
    end
               
    cell2csv([fcfg.out_dir '/' 'Models' '/' 'LogReg' '/' fcfg.mdl_cmp_nme{iP} '/' 'RawData' '_' mmil_spec_char(fcfg.mdl_nme{fcfg.mdl_cmp_plt{iP}(iM)},{';'}) '_' fcfg.mdl_cmp_nme{iP} '_' 'ModelNum_' num2str(fcfg.mdl_cmp_plt{iP}(iSV)) '.csv'],raw_dta)
    
    end
end

end