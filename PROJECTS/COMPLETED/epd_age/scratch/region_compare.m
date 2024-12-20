clear; clc;

%% Load Data
ttt = load('/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/Graph_Scratch/Thickness_Desikan_Modified/GraphData_ClusteringCoefficient.mat');

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'Epilepsy_and_Aging';
dta_hld = 'MRI_thickness_aparc_xsplit_Epilepsy_and_Aging_updated.csv';
thk_dta = mmil_readtext( [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'ROI' '/' dta_hld ] );
    thk_roi_nme = thk_dta(1,2:end);
    clear thk_dta

hld_dta.plt_dta = ttt.plt_dta;

clear ttt
    
%%
[~,roi_lst,~]=fs_read_annotation([ '/home/ekaestne/PROJECTS/DATA/fsaverage/ROIs' '/' 'lh.aparc.split.annot' ]);

mos_ant_tmp_roi = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'rostral-ITG' ...
            'rostral-MTG' ...
            'rostral-STG' ...
            'temporalpole' ...
            'rostral-fusiform' ...
            'entorhinal' } , 'uni' , 0 );
[ ~ , ind ] = intersect( roi_lst , { 'rostral-ITG' ...
            'rostral-MTG' ...
            'rostral-STG' ...
            'temporalpole' ...
            'rostral-fusiform' ...
            'entorhinal' } )
        
ant_tmp_roi = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'middle-ITG'         'rostral-ITG' ...
            'middle-MTG'         'rostral-MTG' ...
            'middle-STG'         'rostral-STG' ...
            'temporalpole' ...
            'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' } , 'uni' , 0 );

tmp_roi = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' 'insula' } , 'uni' , 0 );

ant_frn_roi = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'parsorbitalis'       'parstriangularis'  'parsopercularis' ...
            'caudalmiddlefrontal' 'middle-middlefrontal' 'rostral-middlefrontal' ...
            'caudal-superiorfrontal' 'middle-superiorfrontal' 'rostral-superiorfrontal' ...
            'frontalpole' 'medialorbitofrontal' 'rostralanteriorcingulate' 'caudalanteriorcingulate' 'lateralorbitofrontal' }  , 'uni' , 0 );
         
frn_roi = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'inferior-precentral' 'middle-precentral' 'superior-precentral' ...
            'parsorbitalis'       'parstriangularis'  'parsopercularis' ...
            'caudalmiddlefrontal' 'middle-middlefrontal' 'rostral-middlefrontal' ...
            'caudal-superiorfrontal' 'middle-superiorfrontal' 'rostral-superiorfrontal' ...
            'paracentral' 'frontalpole' 'medialorbitofrontal' 'rostralanteriorcingulate' 'caudalanteriorcingulate' 'lateralorbitofrontal' }  , 'uni' , 0 );
        
occ_tmp = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'lateraloccipital' 'cuneus' 'lingual' 'pericalcarine' }  , 'uni' , 0 );
    [ ~ , occ_tmp_loc , ~ ] = intersect( roi_lst , occ_tmp );

par_tmp = cellfun( @(x) mmil_spec_char(x,{'-'}) , ...
          { 'inferior-postcentral' 'middle-postcentral' 'superior-postcentral' ...
            'supramarginal' 'inferiorparietal' 'superiorparietal' ...
            'precuneus' 'posteriorcingulate' 'isthmuscingulate' } , 'uni' , 0 );


[ numel( intersect( thk_roi_nme , strcat('lhs','_',tmp_roi) ) ) numel(tmp_roi) ]
[ numel( intersect( thk_roi_nme , strcat('lhs','_',frn_roi) ) ) numel(frn_roi) ]
[ numel( intersect( thk_roi_nme , strcat('lhs','_',occ_tmp) ) ) numel(occ_tmp) ]
[ numel( intersect( thk_roi_nme , strcat('lhs','_',par_tmp) ) ) numel(par_tmp) ]

setxor( thk_roi_nme , [ strcat('lhs','_',[ tmp_roi frn_roi occ_tmp par_tmp ]) strcat('rhs','_',[ tmp_roi frn_roi occ_tmp par_tmp ])] )

%%
% 1 = HC, 2 = MCI, 3 = TLE

col_num = 5;

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',tmp_roi) strcat('rhs','_',tmp_roi)] );
tmp_dta{1}{1} = hld_dta.plt_dta.dta{1}{1}( row_ind , col_num );
tmp_dta{1}{2} = hld_dta.plt_dta.dta{1}{2}( row_ind , col_num );
tmp_dta{1}{3} = hld_dta.plt_dta.dta{1}{3}( row_ind , col_num );

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_tmp_roi) strcat('rhs','_',ant_tmp_roi)] );
tmp_dta{2}{1} = hld_dta.plt_dta.dta{1}{1}( row_ind , col_num );
tmp_dta{2}{2} = hld_dta.plt_dta.dta{1}{2}( row_ind , col_num );
tmp_dta{2}{3} = hld_dta.plt_dta.dta{1}{3}( row_ind , col_num );

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_frn_roi) strcat('rhs','_',ant_frn_roi)] );
tmp_dta{3}{1} = hld_dta.plt_dta.dta{1}{1}( row_ind , col_num );
tmp_dta{3}{2} = hld_dta.plt_dta.dta{1}{2}( row_ind , col_num );
tmp_dta{3}{3} = hld_dta.plt_dta.dta{1}{3}( row_ind , col_num );


%
fcfg = [];

fcfg.ydt     = { tmp_dta{1}{1}    tmp_dta{1}{2} tmp_dta{1}{3}  tmp_dta{1}{4} };
fcfg.xdt     = { 1                2             3              4 };

fcfg.fce_col = { rgb('dark grey') rgb('purple') rgb('teal')    rgb('teal') };
fcfg.edg_col = { rgb('black')     rgb('red')    rgb('black')   rgb('black') };

fcfg.xlb = {     'HC'             'aMCI'        'TLE'          'TLE' };
fcfg.ylb = 'Plot Efficiency';

fcfg.ylm = [-.1 1.1];

ejk_scatter(fcfg)

%
mci_sub = tmp_dta{2}{2} - tmp_dta{2}{1};
tle_sub = tmp_dta{2}{3} - tmp_dta{2}{1};

fcfg = [];

fcfg.ydt     = { mci_sub          tle_sub };
fcfg.xdt     = { 1                2             };

fcfg.fce_col = { rgb('purple') rgb('teal') };
fcfg.edg_col = { rgb('red')    rgb('black') };

fcfg.xlb = {    'aMCI'        'TLE' };
fcfg.ylb = 'Plot Efficiency';

fcfg.ylm = [-1.1 1.1];

ejk_scatter(fcfg)

%
[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{3} )
[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{2} )
[ ~ , pvl ] = ttest2( tmp_dta{2}{2} , tmp_dta{2}{3} )
[ ~ , pvl ] = ttest2( mci_sub , tle_sub )

[ ~ , pvl ] = ttest2( tmp_dta{3}{1} , tmp_dta{3}{3} )
[ ~ , pvl ] = ttest2( tmp_dta{3}{1} , tmp_dta{3}{2} )
[ ~ , pvl ] = ttest2( tmp_dta{3}{2} , tmp_dta{3}{3} )

%%
% 1 = HC, 2 = MCI, 3 = TLE
clear tmp_dta

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_tmp_roi) strcat('rhs','_',ant_tmp_roi)] );
tmp_dta{1}{1} = hld_dta.plt_dta.dta{1}{1}( row_ind , col_num );
tmp_dta{1}{2} = hld_dta.plt_dta.dta{1}{2}( row_ind , col_num );
tmp_dta{1}{3} = hld_dta.plt_dta.dta{1}{3}( row_ind , col_num );

tle_tmp_sub = tmp_dta{1}{3} - tmp_dta{1}{1};
mci_tmp_sub = tmp_dta{1}{2} - tmp_dta{1}{1};

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_frn_roi) strcat('rhs','_',ant_frn_roi)] );
tmp_dta{2}{1} = hld_dta.plt_dta.dta{1}{1}( row_ind , col_num );
tmp_dta{2}{2} = hld_dta.plt_dta.dta{1}{2}( row_ind , col_num );
tmp_dta{2}{3} = hld_dta.plt_dta.dta{1}{3}( row_ind , col_num );

tle_frn_sub = tmp_dta{2}{3} - tmp_dta{2}{1};
mci_frn_sub = tmp_dta{2}{2} - tmp_dta{2}{1};

[ ~ , pvl ] = ttest2( tmp_dta{1}{1} , tmp_dta{1}{2} ) % HC v MCI
[ ~ , pvl ] = ttest2( tmp_dta{1}{1} , tmp_dta{1}{3} ) % HC v TLE
[ ~ , pvl ] = ttest2( tmp_dta{1}{2} , tmp_dta{1}{3} ) % MCI v TLE

[ ~ , pvl ] = ttest2( tle_tmp_sub , mci_tmp_sub )

[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{2} ) % HC v MCI
[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{3} ) % HC v TLE
[ ~ , pvl ] = ttest2( tmp_dta{2}{2} , tmp_dta{2}{3} ) % MCI v TLE

%%
% 1 = HC, 2 = MCI, 3 = EO-TLE, 4 = LO-TLE
clear tmp_dta

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_tmp_roi) strcat('rhs','_',ant_tmp_roi)] );
tmp_dta{1}{1} = hld_dta.plt_dta.dta{2}{1}( row_ind , col_num );
tmp_dta{1}{2} = hld_dta.plt_dta.dta{2}{2}( row_ind , col_num );
tmp_dta{1}{3} = hld_dta.plt_dta.dta{2}{3}( row_ind , col_num );
tmp_dta{1}{4} = hld_dta.plt_dta.dta{2}{4}( row_ind , col_num );

[ ~ , row_ind , ~ ] = intersect( thk_roi_nme , [ strcat('lhs','_',ant_frn_roi) strcat('rhs','_',ant_frn_roi)] );
tmp_dta{2}{1} = hld_dta.plt_dta.dta{2}{1}( row_ind , col_num );
tmp_dta{2}{2} = hld_dta.plt_dta.dta{2}{2}( row_ind , col_num );
tmp_dta{2}{3} = hld_dta.plt_dta.dta{2}{3}( row_ind , col_num );
tmp_dta{2}{4} = hld_dta.plt_dta.dta{2}{4}( row_ind , col_num );

[ ~ , pvl ] = ttest2( tmp_dta{1}{1} , tmp_dta{1}{2} ) % HC v MCI
[ ~ , pvl ] = ttest2( tmp_dta{1}{1} , tmp_dta{1}{3} ) % HC v EO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{1}{1} , tmp_dta{1}{4} ) % HC v LO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{1}{2} , tmp_dta{1}{3} ) % MCI v EO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{1}{2} , tmp_dta{1}{4} ) % MCI v LO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{1}{3} , tmp_dta{1}{4} ) % EO-TLE v LO-TLE

[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{2} ) % HC v MCI
[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{3} ) % HC v EO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{2}{1} , tmp_dta{2}{4} ) % HC v LO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{2}{2} , tmp_dta{2}{3} ) % MCI v EO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{2}{2} , tmp_dta{2}{4} ) % MCI v LO-TLE
[ ~ , pvl ] = ttest2( tmp_dta{2}{3} , tmp_dta{2}{4} ) % EO-TLE v LO-TLE


%%
 cfg.reg_inc_plt = [strcat('lhs','_',mos_ant_tmp_roi) strcat('rhs','_',mos_ant_tmp_roi)];

pct_hld(:,2) = {0};
pct_hld( [7 42  45 59 39 34] , 2) = num2cell([7 42  45 59 39 34]);


