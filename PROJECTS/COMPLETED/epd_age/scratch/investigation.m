

%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcd epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epd020_232' '.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/syn02/1/data/MMILDB/carrierm/fsurf/FSURF_022_S_4173_20110913.132811.781000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('ADNI4173 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' '022S4173_238' '.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf/FSURF_EEP006_1_20090127_20090127.123234.218000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('Austin eep006 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'eep006_232' '.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epd020_238' '.png' ] , '-dpng' , '-r200' )
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('Austin epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epd020_232' 'thickness.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/syn02/1/data/MMILDB/carrierm/fsurf/FSURF_022_S_4173_20110913.132811.781000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('ADNI4173 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' '022S4173_238' 'thickness.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcd epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epd020_232' 'thickness.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epd020_238' 'thickness.png' ] , '-dpng' , '-r200' )
close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_EPK104_1_20120127_20120127.131916.531000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epk_238' '.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_EPK104_1_20120127_20120127.131916.531000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epk_238' 'thickness.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf/FSURF_EPK104_1_20120127_20120127.131916.531000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epk_232' '.png' ] , '-dpng' , '-r200' )
close all

%%
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf/FSURF_EPK104_1_20120127_20120127.131916.531000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); ylim([0 14000]); title('mmilmcdRSI epd020 thickness');
print( [ '/home/ekaestne/PROJECTS/OUTPUT/Epilepsy_and_Aging/surface/corticalthickness/Diagnosis/New Folder' '/' 'epk_232' 'thickness.png' ] , '-dpng' , '-r200' )
close all
