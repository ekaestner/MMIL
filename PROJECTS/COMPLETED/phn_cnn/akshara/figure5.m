%% Make Memory Connectomes
% DEMEANED METHOD
mem_dem_cnn = 'memory_demeaned_connections.csv';
mem_cnn = mmil_readtext([prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'data' '/' mem_dem_cnn ]);

fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;
fcfg.out_nme = 'Phentotype_Memory_noBad';

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { mem_cnn(2:end,1)}; % { mem_cnn(2:end,1)    mem_cnn(2:end,3) };
fcfg.edg_wgh = { sqrt(cell2mat(mem_cnn(2:end,2)))*5 }; %  sqrt(abs(cell2mat(mem_cnn(2:end,4))))*5 };
fcfg.edg_col = { rgb('red') }; % { rgb('red')          rgb('blue') };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)