function SL_plot_depth(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial plots work on %s \n'],sbj)

infile = [fcfg.sbj_dat_hld '/' sbj '_overall_data_depth.mat'];

if exist([fcfg.sbj_dat_hld '/' sbj '_overall_data_depth.mat'],'file')
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.sbj_dat_hld '/' sbj '_overall_data_depth.mat'];
    fwv_dat  = ft_func([],cfg);
    
    for iD = 1:numel(fwv_dat.data_name);
        fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.stt_lab = fwv_dat.(fwv_dat.data_name{iD}).label;
        fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.label = fwv_dat.(fwv_dat.data_name{iD}).label;
    end
    
    %% Put together depths
    % Find Probes
    dpt_hld = cellfun(@(x,y) y(1:x(1)-1),regexp(fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.stt_lab,['\d']),fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.stt_lab,'uni',0);
    dpt_tot = unique(dpt_hld);
    
    for iD = 1:numel(dpt_tot)
        dpt{iD} = find(strcmpi(dpt_hld,dpt_tot{iD}));
    end
    
    % Find Location
   
    %% Plot
    for iD = 1:numel(dpt_tot)
        
        % Setup
        chn_grp = zeros(numel(dpt{iD}),3);
        chn_grp(1:numel(dpt{iD}),[1 3]) = repmat(dpt{iD}',1,2);
        
        dat_loc = zeros(numel(dpt{iD}),3);
        dat_loc(1:numel(dpt{iD}),1) = ones(numel(dpt{iD}),1);
        dat_loc(1:numel(dpt{iD}),3) = ones(numel(dpt{iD}),1)*2;
        
        col_ord = cell(numel(dpt{iD}),3);
        col_ord(1:numel(dpt{iD}),[1 3]) = [repmat({{rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')}},numel(dpt{iD}),1) repmat({{rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')}},numel(dpt{iD}),1)];
        
        eve = cell(numel(dpt{iD}),3);
        eve(1:numel(dpt{iD}),[1 3]) = [repmat({[1 2 3 4]},numel(dpt{iD}),1) repmat({[1 2 3 4]},numel(dpt{iD}),1)];
        
        alt_eve = cell(numel(dpt{iD}),3);
        alt_eve(1:numel(dpt{iD}),[1 3]) = [repmat({'trialinfo'},numel(dpt{iD}),1) repmat({'trialinfo'},numel(dpt{iD}),1)];
        
        y_lnk = zeros(numel(dpt{iD}),3);
        y_lnk(1:numel(dpt{iD}),1) = ones(numel(dpt{iD}),1);
        y_lnk(1:numel(dpt{iD}),3) = ones(numel(dpt{iD}),1)*2;
        
        plt_lbl = {};
        
        stt_dat = {};
        stt_col = {};
        
        % Plot
        cfg = [];
        
        cfg.dat = {fwv_dat.(fwv_dat.data_name{1}) fwv_dat.(fwv_dat.data_name{2})};
        
        cfg.type      = 'macro';
        cfg.plt_dim   = [numel(dpt{iD}) 3];
        
        cfg.dat_loc = dat_loc;
        cfg.chn_grp = {chn_grp};
        
        cfg.alt_eve = alt_eve;
        cfg.eve     = eve;
        cfg.lnstyle.col_ord = col_ord;
        cfg.std_err   = 1;
        
        cfg.ttl_num   = zeros(numel(dpt{iD}),3);
        
        cfg.alt_lbl = 'label';
        
        cfg.y_lnk = y_lnk;
        
        cfg.x_lim       = [-0.25 1.25];
        
        cfg.v_lne       = {repmat({[0.000 0.450 0.900]},numel(dpt{iD}),3)};
        cfg.v_lne_col   = {repmat({{rgb('red') rgb('blue') rgb('black')}},numel(dpt{iD}),3)};
        
        cfg.v_lne_wdt = repmat(0.9,numel(dpt{iD}),3);
        cfg.axe_fnt_sze = repmat(5,numel(dpt{iD}),3);
        cfg.axe_lne_sze = repmat(0.9,numel(dpt{iD}),3);
        cfg.ttl_lne_sze = repmat(6,numel(dpt{iD}),3);
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'png';
        cfg.outdir     = [fcfg.sbj_dat_hld '/' 'initial_plot_depth' '/' ];
        cfg.prefix     = [fcfg.sbj_nme '_' dpt_tot{iD}];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end


end