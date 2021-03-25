clear; clc;

%% Clincial
cln_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/clinical_named.csv');

[ ~ , cln_ind ] = sort( cell2mat(cln_dta(:,3)) );
cln_ind = cln_ind(end-4:end);

% Plot
fcfg = [];

for iTR = 1:5
    fcfg.xdt{iTR}     = cln_dta{cln_ind(iTR),3};
    fcfg.xdt_err{iTR} = cln_dta{cln_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = cln_dta{cln_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'ClinicalImportance';

ejk_bar_vert(fcfg)

%% HCV
hcv_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/hcv_named.csv');

[ ~ , hcv_ind ] = sort( cell2mat(hcv_dta(:,3)) );
hcv_ind = hcv_ind(end-1:end);

% Plot
fcfg = [];

for iTR = 1:2
    fcfg.xdt{iTR}     = hcv_dta{hcv_ind(iTR),3};
    fcfg.xdt_err{iTR} = hcv_dta{hcv_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = hcv_dta{hcv_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.60];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'HCVImportance';

ejk_bar_vert(fcfg)

%% Tract
trc_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/tracts_named.csv');

[ ~ , trc_ind ] = sort( cell2mat(trc_dta(:,3)) );
trc_ind = trc_ind(end-4:end);

% Plot
fcfg = [];

for iTR = 1:5
    fcfg.xdt{iTR}     = trc_dta{trc_ind(iTR),3};
    fcfg.xdt_err{iTR} = trc_dta{trc_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = trc_dta{trc_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'TractImportance';

ejk_bar_vert(fcfg)

%% Connectome
cnn_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/connectome_named.csv');

[ ~ , cnn_ind ] = sort( cell2mat(cnn_dta(:,3)) );
cnn_ind = cnn_ind(end-4:end);

% Plot
fcfg = [];

for iTR = 1:5
    fcfg.xdt{iTR}     = cnn_dta{cnn_ind(iTR),3};
    fcfg.xdt_err{iTR} = cnn_dta{cnn_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = cnn_dta{cnn_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'ConnectomeImportance';

ejk_bar_vert(fcfg)

%% Tract+HCV
trc_hcv_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/tracts_hcv_named.csv');

[ ~ , trc_hcv_ind ] = sort( cell2mat(trc_hcv_dta(:,3)) );
trc_hcv_ind = trc_hcv_ind(end-4:end);

% Plot
fcfg = [];

for iTR = 1:5
    fcfg.xdt{iTR}     = trc_hcv_dta{trc_hcv_ind(iTR),3};
    fcfg.xdt_err{iTR} = trc_hcv_dta{trc_hcv_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = trc_hcv_dta{trc_hcv_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'TractHCVImportance';

ejk_bar_vert(fcfg)

%% Connectome+HCV
cnn_hcv_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/connectome_hcv_named.csv');

[ ~ , cnn_hcv_ind ] = sort( cell2mat(cnn_hcv_dta(:,3)) );
cnn_hcv_ind = cnn_hcv_ind(end-4:end);

% Plot
fcfg = [];

for iTR = 1:5
    fcfg.xdt{iTR}     = cnn_hcv_dta{cnn_hcv_ind(iTR),3};
    fcfg.xdt_err{iTR} = cnn_hcv_dta{cnn_hcv_ind(iTR),4};
    fcfg.ydt{iTR}     = iTR;
    fcfg.fce_col{iTR} = rgb('black');
    fcfg.ylb{iTR}     = cnn_hcv_dta{cnn_hcv_ind(iTR),1};
end

fcfg.xlb = 'Feature Importance';
fcfg.xlm = [0 0.35];
  
fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/shuffle/figure';
fcfg.out_nme = 'ConnectomeHCVImportance';

ejk_bar_vert(fcfg)

