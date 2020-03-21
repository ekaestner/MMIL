function [gdd_grd_lab,bdd_grd_lab,gdd_tal_lab,dat] = mmil_pedot_chn(fcfg,dat)

imp_num = mmil_readtext([fcfg.in_dir '/' 'chn_inf' '/' fcfg.sbj_nme '_impedance.csv']);
imp_map = cell2mat(mmil_readtext([fcfg.in_dir '/' 'chn_inf' '/' fcfg.sbj_nme '_impedance_map.csv']));
grd_map = cell2mat(mmil_readtext([fcfg.in_dir '/' 'chn_inf' '/' fcfg.sbj_nme '_grid_map.csv']));

num_ele = floor((size(imp_num,1)-1)/128);
imp_num = imp_num(2:128*num_ele+1,:);

%% Get map layout & make figure
if strcmpi(fcfg.sbj_nme,'bw02_pedot')
    sde_ltr = {{'A' 'C'}};
else
    sde_ltr = {{'A' 'B'} {'C' 'D'}};
end

if ~exist([fcfg.clr_fld '/' 'electrode_info' '/' fcfg.sbj_nme '/']); mkdir([fcfg.clr_fld '/' 'electrode_info' '/' fcfg.sbj_nme '/']); end

for iE = 1:num_ele
    grd_idn{iE} = cell(size(imp_map));
    grd_idn{iE}(imp_map>0) = {sde_ltr{iE}{1}};
    grd_idn{iE}(imp_map<0) = {sde_ltr{iE}{2}};
    grd_lab{iE} = cell(size(imp_map));
    grd_imp{iE} = zeros(size(imp_map));
    
    tmp_lay_out = abs(imp_map);
    
    for iR = 1:size(grd_idn{iE},1);
        for iC = 1:size(grd_idn{iE},2);
            
            if ~isnan(tmp_lay_out(iR,iC))
                grd_lab{iE}{iR,iC} = imp_num{find(strcmpi(imp_num(:,1),[grd_idn{iE}{iR,iC} '-0' num2str(tmp_lay_out(iR,iC))])),1};
                grd_imp{iE}(iR,iC) = imp_num{find(strcmpi(imp_num(:,1),[grd_idn{iE}{iR,iC} '-0' num2str(tmp_lay_out(iR,iC))])),5};
                
                grd_lab{iE}{iR,iC} = imp_num{find(strcmpi(imp_num(:,1),[grd_idn{iE}{iR,iC} '-0' num2str(tmp_lay_out(iR,iC))])),1};
                grd_imp{iE}(iR,iC) = imp_num{find(strcmpi(imp_num(:,1),[grd_idn{iE}{iR,iC} '-0' num2str(tmp_lay_out(iR,iC))])),5};
            end
        end
    end
    
    if ~exist([fcfg.clr_fld '/' 'electrode_info' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_electrode_' num2str(iE) '_impedance.png'])
        figure('units','normalized','outerposition',[0 0 1 1],'Visible','off');
        heatmap_advance(grd_imp{iE}/1000,[],[],'%0.2f','ColorMap',@hot,'colorbar','true','MinColorValue', 0, 'MaxColorValue', 50,'NaNColor', [0 0 0],'FontSize', 26)
        title('Grid 1 Impedance Magnitude/1000 at 1000 Hz (ohms)','FontSize',20)
        hold on
        for iR = 1:size(grd_idn{iE},1);
            for iC = 1:size(grd_idn{iE},2);
                if ~isempty(grd_lab{iE}{iR,iC})
                    text(iC-0.3,iR,grd_lab{iE}{iR,iC},'FontSize',12)
                end
            end
        end
        print(gcf,[fcfg.clr_fld '/' 'electrode_info' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '_electrode_' num2str(iE) '_impedance.png'],'-dpng');
    end
    
end

%% Get channel names
for iE = 1:num_ele
    for iR = 1:size(grd_lab{iE},1);
        for iC = 1:size(grd_lab{iE},2);
            if ~isempty(grd_lab{iE}{iR,iC})
                dat.(fcfg.sbj_nme).label{grd_map(iR,iC)+64*(iE-1)} = grd_lab{iE}{iR,iC};
            else
                
            end
        end
    end
end

emp_chn = find(cellfun(@isempty,dat.(fcfg.sbj_nme).label(:)));
for iFX = 1:numel(emp_chn)
    dat.(fcfg.sbj_nme).label{emp_chn(iFX)} = ['hold_chan_' num2str(iFX)];
end

%% Get good channels
tal_grd_map = grd_map(:,9);
grd_grd_map = grd_map(:,1:8);

for iE = 1:num_ele
    tal_imp{iE} = grd_imp{iE}(:,9);
    grd_imp{iE} = grd_imp{iE}(:,1:8);
    
    gdd_grd_lab{iE} = grd_grd_map((grd_imp{iE}/1000)<60 & ~(grd_imp{iE}==0));
    bdd_grd_lab{iE} = grd_grd_map((grd_imp{iE}/1000)>60 & ~(grd_imp{iE}==0));
    gdd_tal_lab{iE} = tal_grd_map((tal_imp{iE}/1000)<60 & ~(tal_imp{iE}==0));
end

%% Chan Labels
dat.(dat.data_name{1}).label = dat.(dat.data_name{1}).label';
dat.(dat.data_name{1}).cfg.alt_lab.label = dat.(dat.data_name{1}).label;

%% Save information
tal_lab = [];
for iE = 1:numel(grd_lab)
    tal_hld = grd_lab{iE}(:,end);
    tal_hld = tal_hld(~cellfun(@isempty,tal_hld));
    tal_lab = [tal_lab ; tal_hld];
end

cfg         = [];
cfg.dat_nme = dat.data_name;
cfg.clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.alt_lab = 'label';
cfg.typ     = 'pedot';
cfg.tal_lab = tal_lab;
cfg.specific = {'dat_nme' ; 1:numel(dat.data_name)};
ft_func(@mmil_create_depth2,cfg,dat);

% GoodChan
gdd_lab = [];
for iE = 1:numel(grd_lab)
    gdd_hld = dat.(dat.data_name{1}).cfg.alt_lab.label(gdd_grd_lab{iE});
    gdd_hld = [gdd_hld ; dat.(dat.data_name{1}).cfg.alt_lab.label(gdd_tal_lab{iE})];
    gdd_lab = [gdd_lab ; gdd_hld];
end

cfg = [];
cfg.clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.alt_lab = 'label';
cfg.msc_fld = 'gdd_chn';
cfg.msc_dat = ones(numel(gdd_lab),1);
cfg.msc_ord = gdd_lab;
mmil_misc_chan(cfg,dat)

% Impedance
imp_num = [];
imp_lab = [];
for iE = 1:numel(grd_lab)
    grd_imp_hld = grd_imp{iE}(:,1:8);
    grd_imp_hld = grd_imp_hld(:);
    imp_num = [imp_num ; grd_imp_hld];
    
    grd_lab_hld = grd_lab{iE}(:,1:8);
    grd_lab_hld = grd_lab_hld(:);
    imp_lab = [imp_lab ; grd_lab_hld];
end


cfg = [];
cfg.clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;
cfg.alt_lab = 'label';
cfg.msc_fld = 'grd_imp';
cfg.msc_dat = imp_num;
cfg.msc_ord = imp_lab;
mmil_misc_chan(cfg,dat)

end