thk_dta = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/ENIGMA_EPILEPSY/Dalin/T1_3Tepi_voladj_age18_70_EM_miss_fill_No_SA.csv');

dem_dta = thk_dta(:,1:10);
    dem_dta_col = dem_dta(1,:);
    dem_dta     = dem_dta(2:end,:);
    
thk_dta = thk_dta(:,11:end);
    thk_dta_col = thk_dta(1,:);
    thk_dta     = cell2mat(thk_dta(2:end,:));
    
%%
fcfg = [];

fcfg.sbj_nme = strcat( cellfun(@num2str,dem_dta(:,1),'uni',0), dem_dta(:,2));

fcfg.dta     = thk_dta;
fcfg.dta_nme = thk_dta_col;

fcfg.cov     = cell2mat(dem_dta(:,5));
fcfg.cov_nme = dem_dta_col(:,5);

fcfg.out_dir = '/home/ekaestner/Downloads/'; % split % desikan

rsd_mtx_mlt = ejk_residual_matrix( fcfg );

cell2csv( '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/CoAuthor/ENIGMA_EPILEPSY/Dalin/ResidualMatrix.csv',[ dem_dta_col thk_dta_col ; dem_dta num2cell(rsd_mtx_mlt) ] );
