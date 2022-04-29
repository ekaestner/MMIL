clear; clc;

raw_fle = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/data/JohnnyData/raw_data/mwp1epd072.nii';

slc_ind = [ 84 87 90 93 96 ];

plt_out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_lat/manuscript/figures/png/Figure5';

%% Load data
mri_dta = niftiread(raw_fle);
    mri_dta = mri_dta / max(mri_dta(:));
    mri_dta = 1 - mri_dta;

imshow( rot90(squeeze(mri_dta(10:159,84,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" )
    print(gcf,[ plt_out_dir '/' 'slice_84.png'],'-dpng'); close all
imshow( rot90(squeeze(mri_dta(10:159,87,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" )
    print(gcf,[ plt_out_dir '/' 'slice_87.png'],'-dpng'); close all
imshow( rot90(squeeze(mri_dta(10:159,90,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" )
    print(gcf,[ plt_out_dir '/' 'slice_90.png'],'-dpng'); close all
imshow( rot90(squeeze(mri_dta(10:159,93,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" )
    print(gcf,[ plt_out_dir '/' 'slice_93.png'],'-dpng'); close all
imshow( rot90(squeeze(mri_dta(10:159,96,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" )
    print(gcf,[ plt_out_dir '/' 'slice_96.png'],'-dpng'); close all
    
figure('Visible','off'); hold on;
imshow( rot90(squeeze(mri_dta(48,20:195,10:159))), 'InitialMagnification',200,'Interpolation',"bilinear" );
    line([169-84 169-84],[1 149],'Color',rgb('black'),'LineWidth',2)
    line([169-87 169-87],[1 149],'Color',rgb('black'),'LineWidth',2)
    line([169-90 169-90],[1 149],'Color',rgb('black'),'LineWidth',2)
    line([169-93 169-93],[1 149],'Color',rgb('black'),'LineWidth',2)
    line([169-96 169-96],[1 149],'Color',rgb('black'),'LineWidth',2)
    print(gcf,[ plt_out_dir '/' 'slice_selection_reverse.png'],'-dpng'); close all
    



