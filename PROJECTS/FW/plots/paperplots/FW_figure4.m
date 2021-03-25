function FW_figure4

%% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';

fcfg.typ     = {'hgp'};
fcfg.ele_typ = {'ecog'};

fcfg.ana_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.ana_clm = [ 1            2             1 ];
fcfg.eff_clm = { [1 3]        3             1 };

fcfg.ana_col = { 'purple' 'red'  'reddish grey'};
fcfg.ana_lbl = { 'letter' 'word' 'false-font'};

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { { 'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
                 { 'lateraloccipital' } ...
                 { 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'caudal-MTG' 'middle-MTG' 'rostral-MTG' } ...
                 { 'caudal-STG' 'middle-STG' 'rostral-STG' } ...
                 { 'inferior-precentral' 'middle-precentral' } ...
                 { 'parstriangularis' } ...
                 { 'parsopercularis' } };
fcfg.reg_nme = { 'Fusiform' ...
                 'Lateral Occipital' ...
                 'ITG' ...
                 'MTG' ...
                 'STG' ...
                 'Precentral' ...
                 'Pars Triangularis' ...
                 'Pars Operculatris' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize'];
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
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ; 7 3]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% 
fcfg.ana_nme = {'pap_lex_600' 'pap_lex_600' };
fcfg.ana_clm = [ 1            3          ];
fcfg.eff_clm = { 1            3          };
fcfg.ana_col = { 'orange'     'dark yellow'};
fcfg.ana_lbl = { 'repetition' 'lexical'};

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical'];
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
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%:numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
                mmil_tme_plt(cfg,tme_plt{iH})
        
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
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [6 1 ; 4 2 ; 5 3]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Table Construction
% lhs
wrd_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/' 'eff_tme_tbl_' fcfg.hms{1} '_combined.csv']);
lex_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/' 'eff_tme_tbl_' fcfg.hms{1} '_combined.csv']);
tme_tbl = [wrd_tbl lex_tbl];
ttl_col = 2;

col = {'' 'Orthographic-Selective' 'Word-Selective' '' 'Word-Frequency' 'Repetition' '' 'False-Font'};

out_tbl = tme_tbl(:,1);
for iC = 1:numel(col)
    if ~strcmpi(col{iC},'')
        out_tbl(:,iC+1) = tme_tbl(:,string_find(tme_tbl(ttl_col,:),col(iC)));
    end
end

cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table4.csv'],out_tbl);

% rhs
wrd_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/' 'eff_tme_tbl_' fcfg.hms{2} '_combined.csv']);
lex_tbl = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/' 'eff_tme_tbl_' fcfg.hms{2} '_combined.csv']);
tme_tbl = [wrd_tbl lex_tbl];
ttl_col = 2;

col = {'' 'Orthographic-Selective' 'Word-Selective' '' 'Word-Frequency' 'Repetition' '' 'False-Font'};

out_tbl = tme_tbl(:,1);
for iC = 1:numel(col)
    if ~strcmpi(col{iC},'')
        out_tbl(:,iC+1) = tme_tbl(:,string_find(tme_tbl(ttl_col,:),col(iC)));
    end
end

cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table5.csv'],out_tbl);

%% Overall Activity
fcfg.ana_nme = {'pap_rsp_600' };
fcfg.ana_clm = [ 1            ];
fcfg.eff_clm = { 1            };
fcfg.ana_col = { 'black'     };
fcfg.ana_lbl = { 'Specific' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific'];
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
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%:numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
    end
    
end

end




