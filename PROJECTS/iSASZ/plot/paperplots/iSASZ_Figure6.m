function iSASZ_Figure6

clear; clc;

dta_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
out_dir = [clr_fld '/' 'sig_chn' '/' 'amp' '/' ];

%% AUDITORY AMPLITUDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.dta_fld = dta_fld;
fcfg.clr_fld = clr_fld;

fcfg.typ     = {'hgp'};
fcfg.ele_typ = {'ecog'};
fcfg.loc_typ = {'split'};

fcfg.alt_eve = { 'aud_new_bse' }; % 'aud_new_bse' % 'vis_new_bse' 'aud_new_bse'
fcfg.eve     = [ 202       ];     % 202           % 102           202
fcfg.eve_nme = { 'aud_act' };     % 'aud_act'     % 

fcfg.mdl     = 'pap_aud_act'; %
fcfg.out_nme = 'pap_aud_act';

fcfg.dta_typ = 2;

fcfg.zsc_bse = [-0.300 0.000];
fcfg.zsc_epc = [0.050 0.600] ;

out_dir = [ clr_fld '/' 'manuscript' '/' 'figure6' '/' 'voice_amplitude' '/' ];

%% Find Amplitude
sbj_nme = mmil_readtext([fcfg.clr_fld '/' 'subjects']);
mdl_ele = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.typ{1} '/' fcfg.ele_typ{1} '/' fcfg.loc_typ{1} '/' fcfg.mdl '/' 'subjects' '/' 'total' '/' fcfg.mdl '_plt']);
mdl_ele_lab = mdl_ele(1,1:2);
mdl_ele = mdl_ele(2:end,1:2);

for iS = 1:numel(sbj_nme)
    
    cfg = [];
    cfg.load    = 'yes';
    cfg.file    = [fcfg.dta_fld '/' sbj_nme{iS} '_overall_data.mat'];
    bcc_dat     = ft_func([],cfg);
    
    dta_hld = cat(3,bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).trial{:});
%     for iC = 1:size(dta_hld,1)
%         for iT = 1:size(dta_hld,3)
%             dta_zsc_hld(iC,:,iT) = ( dta_hld(iC,:,iT) - nanmean(dta_hld(iC,tme_bse_beg:tme_bse_end,iT)) ) ./ std(dta_hld(iC,tme_bse_beg:tme_bse_end,iT));
%         end
%     end
    
    tme_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).time{1};
    
    tme_bse_beg = dsearchn(tme_hld',fcfg.zsc_bse(1)); tme_bse_end = dsearchn(tme_hld',fcfg.zsc_bse(2));
    tme_epc_beg = dsearchn(tme_hld',fcfg.zsc_epc(1)); tme_epc_end = dsearchn(tme_hld',fcfg.zsc_epc(2));
    
    for iA = 1:numel(fcfg.alt_eve)
        
        if isfield(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve,fcfg.alt_eve{iA})
            
            eve_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve.(fcfg.alt_eve{iA});
            eve_loc = find(eve_hld==fcfg.eve(iA));
            
            for iC = 1:numel(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label)
                
                ele_num = find(strcmpi(mdl_ele(:,2),[sbj_nme{iS} '_' bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label{iC}]));
                
                if ~isempty(ele_num)
                    
                    dta_men = nanmean(dta_hld(iC,:,eve_loc),3);
                    zsc_hld = ( dta_men - nanmean(dta_men(tme_bse_beg:tme_bse_end)) ) ./ std(dta_men(tme_bse_beg:tme_bse_end));
                    
                    mdl_ele{ele_num,iA+2} = max(zsc_hld(tme_epc_beg:tme_epc_end)); % max(nanmean(dta_zsc_hld(iC,tme_epc_beg:tme_epc_end,eve_loc),3));
                    
                end
            end
        end
    end
    
end

%
mmil_chk_dir([fcfg.clr_fld '/' 'sig_chn' '/' 'amp']);

cell2csv([fcfg.clr_fld '/' 'sig_chn' '/' 'amp' '/' fcfg.out_nme ],[mdl_ele_lab fcfg.eve_nme{:} ; mdl_ele])

%% Plot
% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_aud_act'   };
fcfg.ana_clm = [ 1               ];
fcfg.eff_clm = { [1]             };
fcfg.ana_col = { 'dark blue'    };
fcfg.ana_lbl = { 'AuditoryActive' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
                 {'lateraloccipital'} ...
                 {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
                 {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
                 {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
                 {'supramarginal'} ...
                 {'inferior-precentral' 'middle-precentral'} ...
                 {'parstriangularis'} ...
                 {'parsopercularis'} ...
                 {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
                'Lateral Occipital' ...
                'ITG' ...
                'MTG' ...
                'STG' ...
                'Supramarginal' ...
                'Precentral' ...
                'Pars Triangularis' ...
                'Pars Operculatris' ...
                'middle-middlefrontal' };

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sbj_nme = 'total';
    [amp_tab,amp_plt] = mmil_amp_tab(cfg);
    
    for iH = 1:numel(fcfg.hms)
        % Make plot
        cfg = [];
        cfg.typ = fcfg.typ{iT};
        cfg.hms = fcfg.hms{iH};
        cfg.ana_col = fcfg.ana_col;
        cfg.ana_lbl = fcfg.ana_lbl;
        cfg.out_dir = out_dir;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_amp_plt(cfg,amp_plt{iH})
    end
    
end

%% VISUAL AMPLITUDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

fcfg.dta_fld = dta_fld;
fcfg.clr_fld = clr_fld;

fcfg.typ     = {'hgp'};
fcfg.ele_typ = {'ecog'};
fcfg.loc_typ = {'split'};

fcfg.alt_eve = { 'vis_new_bse' }; % 'aud_new_bse' % 'vis_new_bse' 'aud_new_bse'
fcfg.eve     = [ 102       ];     % 202           % 102           202
fcfg.eve_nme = { 'vis_act' };     % 'aud_act'     % 

fcfg.mdl     = 'pap_vis_act'; %
fcfg.out_nme = 'pap_vis_act';

fcfg.dta_typ = 2;

fcfg.zsc_bse = [-0.300 0.000];
fcfg.zsc_epc = [0.050 0.600] ;

out_dir = [ clr_fld '/' 'manuscript' '/' 'figure6' '/' 'text_amplitude' '/' ];

%% Find Amplitude
sbj_nme = mmil_readtext([fcfg.clr_fld '/' 'subjects']);
mdl_ele = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.typ{1} '/' fcfg.ele_typ{1} '/' fcfg.loc_typ{1} '/' fcfg.mdl '/' 'subjects' '/' 'total' '/' fcfg.mdl '_plt']);
mdl_ele_lab = mdl_ele(1,1:2);
mdl_ele = mdl_ele(2:end,1:2);

for iS = 1:numel(sbj_nme)
    
    cfg = [];
    cfg.load    = 'yes';
    cfg.file    = [fcfg.dta_fld '/' sbj_nme{iS} '_overall_data.mat'];
    bcc_dat     = ft_func([],cfg);
    
    dta_hld = cat(3,bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).trial{:});
%     for iC = 1:size(dta_hld,1)
%         for iT = 1:size(dta_hld,3)
%             dta_zsc_hld(iC,:,iT) = ( dta_hld(iC,:,iT) - nanmean(dta_hld(iC,tme_bse_beg:tme_bse_end,iT)) ) ./ std(dta_hld(iC,tme_bse_beg:tme_bse_end,iT));
%         end
%     end
    
    tme_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).time{1};
    
    tme_bse_beg = dsearchn(tme_hld',fcfg.zsc_bse(1)); tme_bse_end = dsearchn(tme_hld',fcfg.zsc_bse(2));
    tme_epc_beg = dsearchn(tme_hld',fcfg.zsc_epc(1)); tme_epc_end = dsearchn(tme_hld',fcfg.zsc_epc(2));
    
    for iA = 1:numel(fcfg.alt_eve)
        
        if isfield(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve,fcfg.alt_eve{iA})
            
            eve_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve.(fcfg.alt_eve{iA});
            eve_loc = find(eve_hld==fcfg.eve(iA));
            
            for iC = 1:numel(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label)
                
                ele_num = find(strcmpi(mdl_ele(:,2),[sbj_nme{iS} '_' bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label{iC}]));
                
                if ~isempty(ele_num)
                    
                    dta_men = nanmean(dta_hld(iC,:,eve_loc),3);
                    zsc_hld = ( dta_men - nanmean(dta_men(tme_bse_beg:tme_bse_end)) ) ./ std(dta_men(tme_bse_beg:tme_bse_end));
                    
                    mdl_ele{ele_num,iA+2} = max(zsc_hld(tme_epc_beg:tme_epc_end)); % max(nanmean(dta_zsc_hld(iC,tme_epc_beg:tme_epc_end,eve_loc),3));
                    
                end
            end
        end
    end
    
end

%
mmil_chk_dir([fcfg.clr_fld '/' 'sig_chn' '/' 'amp']);

cell2csv([fcfg.clr_fld '/' 'sig_chn' '/' 'amp' '/' fcfg.out_nme ],[mdl_ele_lab fcfg.eve_nme{:} ; mdl_ele])

%% Plot
% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'   };
fcfg.ana_clm = [ 1               ];
fcfg.eff_clm = { [1]             };
fcfg.ana_col = { 'dark red'    };
fcfg.ana_lbl = { 'VisualActive' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
                 {'lateraloccipital'} ...
                 {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
                 {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
                 {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
                 {'supramarginal'} ...
                 {'inferior-precentral' 'middle-precentral'} ...
                 {'parstriangularis'} ...
                 {'parsopercularis'} ...
                 {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
                'Lateral Occipital' ...
                'ITG' ...
                'MTG' ...
                'STG' ...
                'Supramarginal' ...
                'Precentral' ...
                'Pars Triangularis' ...
                'Pars Operculatris' ...
                'middle-middlefrontal' };

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sbj_nme = 'total';
    [amp_tab,amp_plt] = mmil_amp_tab(cfg);
    
    for iH = 1:numel(fcfg.hms)
        % Make plot
        cfg = [];
        cfg.typ = fcfg.typ{iT};
        cfg.hms = fcfg.hms{iH};
        cfg.ana_col = fcfg.ana_col;
        cfg.ana_lbl = fcfg.ana_lbl;
        cfg.out_dir = out_dir;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_amp_plt(cfg,amp_plt{iH})
    end
    
end

%% PLOT PIECES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end




