function SL_figure6

%% VISUAL
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_lng_950'   'pap_con_950'        };
fcfg.ana_clm = [ 1               1                    ];
fcfg.eff_clm = { [1]             [1]                  };
fcfg.ana_col = { 'bright red'    'dark magenta'       };
fcfg.ana_lbl = { 'TextSelective' 'FalseFontSelective' };

fcfg.sub_tme = [ 0.000           0.000                ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
    { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
    { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

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

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Visual' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% AUDITORY
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_lng_950'   'pap_con_950'             'pap_mtc_1450' };
fcfg.ana_clm = [ 2               2                         1            ];
fcfg.eff_clm = { [2]             [2]                       [1]            };
fcfg.ana_col = { 'blue'          'cyan'                    'dark yellow'       };
fcfg.ana_lbl = { 'VoiceSelective' 'NoiseVocodingSelective' 'Mismatch'     };

fcfg.sub_tme = [ 0.450            0.450                    0.450          ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
    { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
    { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

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

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Auditory' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [4 1 ; 5 2 ; 6 3 ]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Table Construction
fcfg.ana_nme = { 'pap_lng_950'   'pap_con_950'        'pap_phn_950'      'pap_lng_950'    'pap_con_950'            'pap_phn_950'       'pap_mtc_1450' };
fcfg.ana_clm = [ 1               1                    1                  2                2                        2                   1              ];
fcfg.eff_clm = { [1]             [1]                  [1]                [2]              [2]                      [2]                 [1]            };
fcfg.ana_col = { 'red'           'bright magenta'     'reddish orange'   'blue'           'cyan'                   'bluish purple'     'yellow'       };
fcfg.ana_lbl = { 'TextSelective' 'FalseFontSelective' 'Letter-Sensitive' 'VoiceSelective' 'NoiseVocodingSelective' 'Phoneme-Sensitive' 'Mismatch'     };
fcfg.sub_tme = [ 0.000           0.000                0.000              0.450            0.450                    0.450               0.450          ];

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Combined' '/'];
mmil_chk_dir(fcfg.out_dir)

for iT = 1:numel(fcfg.typ)
    
    %% Load Data
     cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        %         mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        tme_hld = 0:0.020:2;
        for iR = 1:size(tme_tab{iH},1)
            for iC = 2:size(tme_tab{iH},2)
                if ~isempty(tme_tab{iH}{iR,iC})
                    for iTM = 1:numel(tme_tab{iH}{iR,iC})
                        tme_tab{iH}{iR,iC}(iTM) = tme_hld(dsearchn(tme_hld',tme_tab{iH}{iR,iC}(iTM)));
                    end
                end
            end
        end
        
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_inclusive.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [ 2 1 ; 3 2 ; 4 3 ; 5 4 ; 6 5 ; 7 6 ; 8 7 ]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

% lhs
tme_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Combined' '/' 'eff_tme_tbl_' fcfg.hms{1} '_combined.csv']);
ttl_col = 2;

col = {'' 'Text Selective' 'False-Font Selective' 'Letter-Sensitive' '' 'Voice Selective' 'Noise-Vocoded Selective' 'Phoneme-Sensitive' '' 'Incongruent Effects'};

out_tbl = tme_tbl(:,1);
for iC = 1:numel(col)
    if ~strcmpi(col{iC},'')
        out_tbl(:,iC+1) = tme_tbl(:,string_find(tme_tbl(ttl_col,:),col(iC)));
    end
end

cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table4.csv'],out_tbl);

% rhs
tme_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Combined' '/' 'eff_tme_tbl_' fcfg.hms{2} '_combined.csv']);
ttl_col = 2;

col = {'' 'Text Selective' 'False-Font Selective' 'Letter-Sensitive' '' 'Voice Selective' 'Noise-Vocoded Selective' 'Phoneme-Sensitive' '' 'Incongruent Effects'};

out_tbl = tme_tbl(:,1);
for iC = 1:numel(col)
    if ~strcmpi(col{iC},'')
        out_tbl(:,iC+1) = tme_tbl(:,string_find(tme_tbl(ttl_col,:),col(iC)));
    end
end

cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table5.csv'],out_tbl);

%% Overall Activity
% fcfg.ana_nme = {'pap_anv_1500' };
% fcfg.ana_clm = [ 4            ];
% fcfg.eff_clm = { 4            };
% fcfg.ana_col = { 'black'     };
% fcfg.ana_lbl = { 'Specific' };
% 
% fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'specific'];
% mmil_chk_dir(fcfg.out_dir)
% 
% for iT = 1:numel(fcfg.typ)
%     
%     %% Load Data
%     cfg = [];
%     cfg.clr_fld = fcfg.clr_fld;
%     cfg.typ     = fcfg.typ{iT};
%     cfg.ele_typ = fcfg.ele_typ{1};
%     cfg.ana_nme = fcfg.ana_nme;
%     cfg.ana_clm = fcfg.ana_clm;
%     cfg.eff_clm = fcfg.eff_clm;
%     cfg.sbj_nme = 'total';
%     [tme_tab,tme_plt] = mmil_tme_tab(cfg);
%     
%     save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
%     
%     %% Plot Figures
%     for iH = 1:numel(fcfg.hms)
%         for iP = 1%:numel(fcfg.nme_reg)
%             
%             cfg = [];
%             cfg.typ = fcfg.typ{iT};
%             cfg.reg = fcfg.inc_reg{iP};
%             cfg.nme = fcfg.nme_reg{iP};
%             cfg.hms = fcfg.hms{iH};
%             cfg.ana_col = fcfg.ana_col;
%             cfg.ana_lbl = fcfg.ana_lbl;
%             cfg.out_dir = fcfg.out_dir;
%             %mmil_tme_plt(cfg,tme_plt{iH})
%             
%         end
%         
%         cfg.reg_nme = fcfg.reg_nme;
%         cfg.cmb_reg = fcfg.cmb_reg;
%         mmil_tme_plt(cfg,tme_plt{iH})
%         
%     end
%     
% end

%% Overall Activity
% fcfg.ana_nme = {'pap_phn_950' 'pap_phn_950' };
% fcfg.ana_clm = [ 1            2 ];
% fcfg.eff_clm = { 1            2 };
% fcfg.ana_col = { 'red'        'blue' };
% fcfg.sub_tme = [ 0.000        0.450  ];
% fcfg.ana_lbl = { 'Letter'     'Phoneme' };
% 
% fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme'];
% mmil_chk_dir(fcfg.out_dir)
% 
% for iT = 1:numel(fcfg.typ)
%     
%     %% Load Data
%     cfg = [];
%     cfg.clr_fld = fcfg.clr_fld;
%     cfg.typ     = fcfg.typ{iT};
%     cfg.ele_typ = fcfg.ele_typ{1};
%     cfg.ana_nme = fcfg.ana_nme;
%     cfg.ana_clm = fcfg.ana_clm;
%     cfg.eff_clm = fcfg.eff_clm;
%     cfg.sub_tme = fcfg.sub_tme;
%     cfg.sbj_nme = 'total';
%     [tme_tab,tme_plt] = mmil_tme_tab(cfg);
%     
%     save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
%     
%     %% Plot Figures
%     for iH = 1:numel(fcfg.hms)
%         for iP = 1%:numel(fcfg.nme_reg)
%             
%             cfg = [];
%             cfg.typ = fcfg.typ{iT};
%             cfg.reg = fcfg.inc_reg{iP};
%             cfg.nme = fcfg.nme_reg{iP};
%             cfg.hms = fcfg.hms{iH};
%             cfg.ana_col = fcfg.ana_col;
%             cfg.ana_lbl = fcfg.ana_lbl;
%             cfg.out_dir = fcfg.out_dir;
%             %mmil_tme_plt(cfg,tme_plt{iH})
%             
%         end
%         
%         cfg.reg_nme = fcfg.reg_nme;
%         cfg.cmb_reg = fcfg.cmb_reg;
%         mmil_tme_plt(cfg,tme_plt{iH})
%         
%     end
%     
% end

%% ALL EFFECTS (FOR STAT PURPOSES)
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_lng_950'   'pap_phn_950'     'pap_con_950'        'pap_lng_950'    'pap_phn_950'      'pap_con_950'            'pap_mtc_1450'  };
fcfg.ana_clm = [ 1               1                 1                    2                2                  2                        1               ];
fcfg.eff_clm = { [1]             [1]               [1]                  [2]              [2]                [2]                      [1]             };
fcfg.ana_col = { 'bright red'    'red'             'dark magenta'       'blue'           'blue'             'cyan'                   'dark yellow' 	};
fcfg.ana_lbl = { 'TextSelective' 'Letter-Specific' 'FalseFontSelective' 'VoiceSelective' 'Phoneme-specific' 'NoiseVocodingSelective' 'Mismatch'      };

fcfg.sub_tme = [ 0.000           0.000             0.000                0.450            0.450              0.450                     0.450           ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
    { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
    { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

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

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [4 1 ; 5 2 ; 6 3 ]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

end




