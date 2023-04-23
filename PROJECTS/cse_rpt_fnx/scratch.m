clear; clc;

ovr_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/iEEG/cse_rpt_fnx';

%% Load data & Collate
% Fornix NP data
fcfg = [];
fcfg.dta_loc = [ ovr_dir '/' 'data' '/' 'Scores_2022_11_10_updated.csv' ];
fcfg.dta_col = 2;
[ fnx_dta, fnx_dta_sbj, fnx_dta_col] = ejk_dta_frm( fcfg );
    fnx_dta = [ fnx_dta_sbj fnx_dta ];
    fnx_dta_col = [ 'test' fnx_dta_col ];

fnx_dta(cellfun(@isempty,fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1'))),strcmpi(fnx_dta_col,'Subtest_T1')) = fnx_dta(cellfun(@isempty,fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1'))),strcmpi(fnx_dta_col,'Subtest_T2'));

ovr_col = strcat(fnx_dta( :, strcmpi(fnx_dta_col,'Test')),'_',fnx_dta( :, strcmpi(fnx_dta_col,'Subtest_T1')));

% Redcap
dta_dir = '/home/ekaestne/PROJECTS/DATA';
red_cap_fle = [ dta_dir '/' 'csv' '/' 'redcap' '/' 'Redcap_2022_07_22.csv'];

fcfg = [];
fcfg.red_fle = red_cap_fle;
fcfg.sep     = '|';
[~ , ~ , ~ , sbj_cog, ~] = mmil_load_redcap(fcfg);

%% Cipher
% Save out fornix half of cipher
cell2csv([ ovr_dir '/' 'data' '/' 'cipher_fnx.csv' ],ovr_col)

% Load modified cipher


%% Normative scores plotting
nor_plt = { 'Age_Sex' 'Age_Sex_Education' };
nor_nme = { { 'T1 (Age,Sex)' 'T2 (Age,Sex)' } {} };
nor_col = { { 'T1_tscore' 'T2_tscore' }       {} };

typ_nme = { 'List_Encoding' 'List_Retrieval' };
typ_scr     = { { 'CVLT-II_List A Trial 1' 'CVLT-II_List A Trial 2' 'CVLT-II_List A Trial 3' 'CVLT-II_List A Trial 4' 'CVLT-II_List A Trial 5' 'CVLT-II_List A 1-5 Total' 'CVLT-II_List B' } ...
                { 'CVLT-II_Short Delay Free Recall' 'CVLT-II_Short Delay Cued Recall' 'CVLT-II_Long Delay Free Recal' 'CVLT-II_Short Delay Free Recall'  } };
typ_scr_nme = { { 'Trial 1'                'Trial 2'                'Trial 3'                'Trial 4'                'Trial 5'                'Trial Total'              'List B'         } ...
                { 'CVLT_SDFR'                       'CVLT_SDCR'                       'CVLT_LDFR'                     'CVLT_LDCR'                        } };

% Norm Plotting %%%%%%%%%%%%%%%%%%%%%%%%            
iNO = 1;
iTY = 2;

xdt = num2cell([[1:numel(typ_scr_nme{iTY})] [1:numel(typ_scr_nme{iTY})]]);


[~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
ydt = [[(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{1}))-50)./10]' [(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{2}))-50)./10]'];

fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ydt;

fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY}))];
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = typ_scr_nme{iTY};
fcfg.ylb = { 'z-score' };
fcfg.ttl = typ_nme{iTY};

fcfg.hln     = -1;
fcfg.hln_col = rgb('red');

fcfg.jtr = 1;
fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+0.5 ];
fcfg.ylm = [ -4 4 ];

fcfg.out_dir = [ ovr_dir '/' 'plots' '/' nor_plt{iNO} '/' ];
fcfg.out_nme = typ_nme{iTY};

ejk_scatter(fcfg)

%% Normative Data Group Plotting
typ_nme = { 'CVLT' };
typ_scr     = { { 'CVLT-II_List A 1-5 Total' 'CVLT-II_Long Delay Free Recal'   } };
typ_scr_nme = { { 'CVLT_Total'               'CVLT_LDFR'           } };

raw_col = { 'T1_raw' 'T2_raw' };
nor_plt = { 'Age_Sex' 'Age_Sex_Education' };
nor_col = { { 'T1_tscore' 'T2_tscore' }       {} };

cph_typ = mmil_readtext([ ovr_dir '/' 'data' '/' 'cipher_bth.csv' ]);

% Norm Group Plotting %%%%%%%%%%%%%%%%%%%%%%%%
iNO = 1;
iTY = 1;

xdt = num2cell([[1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]+.1 ]);

[~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
ydt = [[(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{1}))-50)./10]' [(fnx_dta(ydt_row,strcmpi(fnx_dta_col,nor_col{iNO}{2}))-50)./10]'];
xlb = [];
for iRA = 1:numel(typ_scr{iTY})
    ydt = [ ydt (sbj_cog.(cph_typ{strcmpi(cph_typ(:,1),typ_scr{iTY}{iRA}),3})(string_find(sbj_cog.sbj_nme,'fc'))-50)./10 ];
    xlb = [ xlb typ_scr_nme{iTY}(iRA) {'HC'} ];
end

fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ydt;

fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY})) repmat({rgb('grey')},1,numel(typ_scr_nme{iTY}))];
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = xlb;
fcfg.ylb = { 'Z-Score' };
fcfg.ttl = typ_nme{iTY};

fcfg.hln     = -1;
fcfg.hln_col = rgb('red');

fcfg.jtr = 1;
fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+0.5 ];
fcfg.ylm = [ -4 4 ];

fcfg.out_dir = [ ovr_dir '/' 'plots' '/' nor_plt{iNO} '/' ];
fcfg.out_nme = [ typ_nme{iTY} '_' 'group' ];

ejk_scatter(fcfg)

% Raw Group Plotting %%%%%%%%%%%%%%%%%%%%%%%%
raw_col = { 'T1_raw' 'T2_raw' };

cph_typ = mmil_readtext([ ovr_dir '/' 'data' '/' 'cipher_bth.csv' ]);


iTY = 1;

xdt = num2cell([[1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]-.1 [1:numel(typ_scr_nme{iTY})]+.1 ]);

[~, ydt_row ] = ismember(typ_scr{iTY},ovr_col);
ydt = [[fnx_dta(ydt_row,strcmpi(fnx_dta_col,raw_col{1}))]' [fnx_dta(ydt_row,strcmpi(fnx_dta_col,raw_col{2}))]'];
xlb = [];
for iRA = 1:numel(typ_scr{iTY})
    ydt = [ ydt sbj_cog.(cph_typ{strcmpi(cph_typ(:,1),typ_scr{iTY}{iRA}),2})(string_find(sbj_cog.sbj_nme,'fc')) ];
    xlb = [ xlb typ_scr_nme{iTY}(iRA) {'HC'} ];
end

fcfg = [];

fcfg.xdt     = xdt;
fcfg.ydt     = ydt;

fcfg.fce_col = [ repmat({rgb('light orange')},1,numel(typ_scr_nme{iTY})) repmat({rgb('purple')},1,numel(typ_scr_nme{iTY})) repmat({rgb('grey')},1,numel(typ_scr_nme{iTY}))];
fcfg.edg_col = repmat({rgb('black')},1,numel(xdt));

fcfg.xlb = xlb;
fcfg.ylb = { 'Raw Score' };
fcfg.ttl = typ_nme{iTY};

fcfg.jtr = 1;
fcfg.xlm = [ 0.5 numel(typ_scr_nme{iTY})+1 ];

fcfg.out_dir = [ ovr_dir '/' 'plots' '/' 'Raw' '/' ];
fcfg.out_nme = typ_nme{iTY};

ejk_scatter(fcfg)








