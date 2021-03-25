function iSASZ_plot_depth(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial plots work on %s \n'],sbj)

infile = [fcfg.dat_fld '/' 'epoch_data' '/' 'depth' '/' sbj '_overall_data_depth.mat'];

if exist([fcfg.dat_fld '/' 'epoch_data' '/' 'depth' '/' sbj '_overall_data_depth.mat'],'file')
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.dat_fld '/' 'epoch_data' '/' 'depth' '/' sbj '_overall_data_depth.mat'];
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
        if size(dpt{iD},1) > size(dpt{iD},2); dpt{iD} = dpt{iD}'; end
    end
    
    % Find Location
    
    %% Plot
    for iD = 1:numel(dpt_tot)
      
        % Setup
        chn_grp = zeros(numel(dpt{iD}),9);
        chn_grp(1:numel(dpt{iD}),[1:4 6:9]) = repmat(dpt{iD}',1,8);
        
        dat_loc = zeros(numel(dpt{iD}),9);
        dat_loc(1:numel(dpt{iD}),1:4) = ones(numel(dpt{iD}),4);
        dat_loc(1:numel(dpt{iD}),6:9) = ones(numel(dpt{iD}),4)*2;
        
        col_ord = cell(numel(dpt{iD}),9);
        col_ord(1:numel(dpt{iD}),[1:4 6:9]) = [repmat({{rgb('red')  rgb('orange')}},numel(dpt{iD}),1) repmat({{rgb('lime') rgb('reddish grey')}},numel(dpt{iD}),1) repmat({{rgb('blue')  rgb('orange')}},numel(dpt{iD}),1) repmat({{rgb('aqua') rgb('bluish grey')}},numel(dpt{iD}),1) repmat({{rgb('red')  rgb('orange')}},numel(dpt{iD}),1) repmat({{rgb('lime') rgb('reddish grey')}},numel(dpt{iD}),1) repmat({{rgb('blue')  rgb('orange')}},numel(dpt{iD}),1) repmat({{rgb('aqua') rgb('bluish grey')}},numel(dpt{iD}),1) ];
        
        eve = cell(numel(dpt{iD}),9);
        eve(1:numel(dpt{iD}),[1:4 6:9]) = [repmat({[111 112]},numel(dpt{iD}),1) repmat({[121 122]},numel(dpt{iD}),1) repmat({[211 212]},numel(dpt{iD}),1) repmat({[221 222]},numel(dpt{iD}),1) repmat({[111 112]},numel(dpt{iD}),1) repmat({[121 122]},numel(dpt{iD}),1) repmat({[211 212]},numel(dpt{iD}),1) repmat({[221 222]},numel(dpt{iD}),1) ];
        
        alt_eve = cell(numel(dpt{iD}),9);
        alt_eve(1:numel(dpt{iD}),[1:4 6:9]) = [repmat({'vis_new_old'},numel(dpt{iD}),1) repmat({'vis_ani_obj'},numel(dpt{iD}),1) repmat({'aud_new_old'},numel(dpt{iD}),1) repmat({'aud_ani_obj'},numel(dpt{iD}),1) repmat({'vis_new_old'},numel(dpt{iD}),1) repmat({'vis_ani_obj'},numel(dpt{iD}),1) repmat({'aud_new_old'},numel(dpt{iD}),1) repmat({'aud_ani_obj'},numel(dpt{iD}),1)];
        
        alt_stt = cell(numel(dpt{iD}),9);
        alt_stt(1:numel(dpt{iD}),[1:4 6:9]) = [repmat({{'vis_new_ovr_stt' 'vis_new_old_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'vis_new_ovr_stt' 'vis_ani_obj_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'aud_new_ovr_stt' 'aud_new_old_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'aud_new_ovr_stt' 'aud_ani_obj_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'vis_new_ovr_stt' 'vis_new_old_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'vis_new_ovr_stt' 'vis_ani_obj_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'vis_new_ovr_stt' 'vis_new_old_stt'}},numel(dpt{iD}),1) ...
                                               repmat({{'aud_new_ovr_stt' 'aud_ani_obj_stt'}},numel(dpt{iD}),1) ];
        
        stt_col = cell(numel(dpt{iD}),9);
        stt_col(1:numel(dpt{iD}),[1:4 6:9]) = [repmat({{rgb('black') rgb('red')}},numel(dpt{iD}),1)  ...
                                               repmat({{rgb('black') rgb('lime')}},numel(dpt{iD}),1) ...
                                               repmat({{rgb('black') rgb('blue')}},numel(dpt{iD}),1) ...
                                               repmat({{rgb('black') rgb('aqua')}},numel(dpt{iD}),1) ...
                                               repmat({{rgb('black') rgb('red')}},numel(dpt{iD}),1)  ...
                                               repmat({{rgb('black') rgb('lime')}},numel(dpt{iD}),1) ...
                                               repmat({{rgb('black') rgb('blue')}},numel(dpt{iD}),1) ...
                                               repmat({{rgb('black') rgb('aqua')}},numel(dpt{iD}),1) ];
                
        stt_cmp = cell(numel(dpt{iD}),9);
        stt_cmp(1:numel(dpt{iD}),[1:4 6:9]) =  [repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...    
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ...
                                                repmat({{'0%7' '7%14'}},numel(dpt{iD}),1) ];
        
        y_lnk = zeros(numel(dpt{iD}),9);
        y_lnk(1:numel(dpt{iD}),1:4) = ones(numel(dpt{iD}),4);
        y_lnk(1:numel(dpt{iD}),6:9) = ones(numel(dpt{iD}),4)*2;
        
        plt_lbl = {};
                
        % Fix for events that are present
        if ~isfield(fwv_dat.(fwv_dat.data_name{1}).cfg.alt_eve,'vis_new_old')
            chn_grp(1:numel(dpt{iD}),[1:2 6:7]) = zeros(numel(dpt{iD}),4);
            dat_loc(1:numel(dpt{iD}),[1:2 6:7]) = zeros(numel(dpt{iD}),4);
            col_ord(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            eve(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            alt_eve(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            y_lnk(1:numel(dpt{iD}),[1:2 6:7]) = zeros(numel(dpt{iD}),4);
            alt_stt(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            stt_col(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            stt_cmp(1:numel(dpt{iD}),[1:2 6:7]) = cell(numel(dpt{iD}),4);
            
            chn_grp = chn_grp(1:numel(dpt{iD}),[3:9]);
            dat_loc = dat_loc(1:numel(dpt{iD}),[3:9]);
            col_ord = col_ord(1:numel(dpt{iD}),[3:9]);
            eve = eve(1:numel(dpt{iD}),[3:9]);
            alt_eve = alt_eve(1:numel(dpt{iD}),[3:9]);
            y_lnk = y_lnk(1:numel(dpt{iD}),[3:9]);
            alt_stt = alt_stt(1:numel(dpt{iD}),[3:9]);
            stt_col = stt_col(1:numel(dpt{iD}),[3:9]);
            stt_cmp = stt_cmp(1:numel(dpt{iD}),[3:9]);
            
            num = 7;
        end
        
        if ~isfield(fwv_dat.(fwv_dat.data_name{1}).cfg.alt_eve,'aud_new_old')
            chn_grp(1:numel(dpt{iD}),[3:4 8:9]) = zeros(numel(dpt{iD}),4);
            dat_loc(1:numel(dpt{iD}),[3:4 8:9]) = zeros(numel(dpt{iD}),4);
            col_ord(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            eve(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            alt_eve(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            y_lnk(1:numel(dpt{iD}),[3:4 8:9]) = zeros(numel(dpt{iD}),4);

            alt_stt(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            stt_col(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            stt_cmp(1:numel(dpt{iD}),[3:4 8:9]) = cell(numel(dpt{iD}),4);
            
            num = 9;
        end
        
        if isfield(fwv_dat.(fwv_dat.data_name{1}).cfg.alt_eve,'aud_new_old') && isfield(fwv_dat.(fwv_dat.data_name{1}).cfg.alt_eve,'vis_new_old')
            num = 9;
        end
        
        dat_loc = [dat_loc ; zeros(max([numel(dpt{iD}) num])-size(dat_loc,1),size(dat_loc,2))];
        
        % Plot
        cfg = [];
        
        cfg.dat = {fwv_dat.(fwv_dat.data_name{1}) fwv_dat.(fwv_dat.data_name{2})};
        
        cfg.type      = 'macro';
        cfg.plt_dim   = [max([numel(dpt{iD}) num]) num];
        
        cfg.dat_loc = dat_loc;
        cfg.chn_grp = {[ chn_grp ; zeros(max([numel(dpt{iD}) num])-size(chn_grp,1),size(chn_grp,2))]};
        
        cfg.alt_eve = [ alt_eve ; cell(max([numel(dpt{iD}) num])-size(alt_eve,1),size(alt_eve,2)) ];
        cfg.eve     = [ eve ; cell(max([numel(dpt{iD}) num])-size(eve,1),size(eve,2))];
        cfg.lnstyle.col_ord = [col_ord ; cell(max([numel(dpt{iD}) num])-size(col_ord,1),size(col_ord,2))];
        cfg.std_err   = 1;
        
        cfg.stt_dat = {[alt_stt ; cell(max([numel(dpt{iD}) num])-size(alt_stt,1),size(alt_stt,2))]};
        cfg.stt_col = {[stt_col ; cell(max([numel(dpt{iD}) num])-size(stt_col,1),size(stt_col,2))]};
        cfg.stt_cmp = {[stt_cmp ; cell(max([numel(dpt{iD}) num])-size(stt_cmp,1),size(stt_cmp,2))]};
        
        cfg.ttl_num   = zeros(max([numel(dpt{iD}) num]),num);
        
        cfg.alt_lbl = 'label';
        
        cfg.y_lnk = y_lnk;
        
        cfg.x_lim       = [-0.15 0.75];
        
        cfg.v_lne       = {repmat({[0.000 0.200 0.400]},max([numel(dpt{iD}) num]),num)};
        cfg.v_lne_col   = {repmat({{rgb('red') rgb('black') rgb('black')}},numel(dpt{iD}),num)};
        
        cfg.v_lne_wdt = repmat(0.9,max([numel(dpt{iD}) num]),num);
        cfg.axe_fnt_sze = repmat(5,max([numel(dpt{iD}) num]),num);
        cfg.axe_lne_sze = repmat(0.9,max([numel(dpt{iD}) num]),num);
        cfg.ttl_lne_sze = repmat(6,max([numel(dpt{iD}) num]),num);
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'png';
        cfg.outdir     = [fcfg.sbj_dat_hld '/' 'initial_plot_depth' '/' ];
        cfg.prefix     = [fcfg.sbj_nme '_' dpt_tot{iD}];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end

end