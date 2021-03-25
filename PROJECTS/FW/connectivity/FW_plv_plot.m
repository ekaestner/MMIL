clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

sbj_tme = [];

for sbj_num = 1:numel(sbj_nme_hld);
    
    tic
        
    %% All Plot
    pcfg = [];
    pcfg.clr_fld = sbj_clr_hld{sbj_num};
    
    pcfg.eve = [1 3 4];
    pcfg.ttl = {'word' 'consonantstring' 'false-font'};
    
    pcfg.foi = [4 8];
    pcfg.mve_tme = [0.050 0.600];
    
    pcfg.mve_tle = 'WordFocus';
    
    pcfg.lbl     = {'Visual'      ''};
    pcfg.lbl_tme = {[0.000 0.400] [0.400 0.600]};
    
    pcfg.str_inc = 'focus'; % 'all','focus'
    pcfg.nde_nme = { { 'caudal-fusiform_L'     'caudal-fusiform_R'     'lateraloccipital_L'  'lateraloccipital_R' } ...
                     { 'inferior-precentral_L' 'inferior-precentral_R' 'middle-precentral_L' 'middle-precentral_R' 'parsopercularis_L' 'parsopercularis_R' } ...
                     { 'caudal-STG_L'          'caudal-STG_R'          'middle-STG_L'        'middle-STG_R' 'supramarginal_L' 'supramarginal_R'} };
    pcfg.nde_col = { 'red'    ...
                     'yellow' ...
                     'blue'   };
    
    pcfg.edg_nme = { [1 2]    [1 3]    [2 3]        [1 1] [2 2]    [3 3] };
    pcfg.edg_col = { 'orange' 'purple' 'aquamarine' 'red' 'yellow' 'blue'};
    
    pcfg.lne_mlt = 5;
    
    pcfg.sbj_num = sbj_num;
    NetworkPlotting2(pcfg)
    
    %% Focused Plot
    pcfg = [];
    pcfg.clr_fld = sbj_clr_hld{sbj_num};
    
    pcfg.eve = [1 3 4];
    pcfg.ttl = {'word' 'consonantstring' 'false-font'};
    
    pcfg.foi = [4 8];
    pcfg.mve_tme = [0.050 0.600];
    
    pcfg.mve_tle = 'Word';
    
    pcfg.lbl     = {'Visual'      ''};
    pcfg.lbl_tme = {[0.000 0.400] [0.400 0.600]};
    
    pcfg.str_inc = 'all'; % 'all','focus'
    pcfg.nde_nme = { { 'caudal-fusiform_L'     'caudal-fusiform_R'     'lateraloccipital_L'  'lateraloccipital_R' } ...
                     { 'inferior-precentral_L' 'inferior-precentral_R' 'middle-precentral_L' 'middle-precentral_R' 'parsopercularis_L' 'parsopercularis_R' } ...
                     { 'caudal-STG_L'          'caudal-STG_R'          'middle-STG_L'        'middle-STG_R' 'supramarginal_L' 'supramarginal_R'} };
    pcfg.nde_col = { 'red'    ...
                     'yellow' ...
                     'blue'   };
    
    pcfg.edg_nme = { [1 2]    [1 3]    [2 3]        [1 1] [2 2]    [3 3] };
    pcfg.edg_col = { 'orange' 'purple' 'aquamarine' 'red' 'yellow' 'blue'};
    
    pcfg.lne_mlt = 5;
    
    pcfg.sbj_num = sbj_num;
    NetworkPlotting2(pcfg)
    
    % Timekeeping
    sbj_tme(sbj_num) = toc;
    
end