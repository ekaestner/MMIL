clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm/';

dta_dir = [ prj_dir '/' 'Data' '/' '2023_03_27' '/'];

atl_dir = [ prj_dir '/' 'Data/atlas/aal_bonilha'];
atl_nme = 'aal';

plt_dir = [ prj_dir '/' 'Output' '/'  'hippocampal slice' '/'];
    
%%
load([ dta_dir '/' 'performance.mat' ])
load([ dta_dir '/' 'saliency.mat' ])

%% Matching slice indices
figure()
subplot(2,2,1)
iR = 20;    
atl_img_plt = squeeze(nft_atl(:,iR,:));
atl_img_plt(atl_img_plt==0) = 255;
atl_img_plt(atl_img_plt<254) = 175;
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,2)
iR = 45;    
atl_img_plt = squeeze(nft_atl(:,iR,:));
atl_img_plt(atl_img_plt==0) = 255;
atl_img_plt(atl_img_plt<254) = 175;
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,3)
iR = 72;    
atl_img_plt = squeeze(nft_atl(:,iR,:));
atl_img_plt(atl_img_plt==0) = 255;
atl_img_plt(atl_img_plt<254) = 175;
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,4)
iR = 120;    
atl_img_plt = squeeze(nft_atl(:,iR,:));
atl_img_plt(atl_img_plt==0) = 255;
atl_img_plt(atl_img_plt<254) = 175;
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

%% Saliency
figure()
subplot(2,2,1)
iR = 20;    
atl_img_plt = mean(squeeze(sal_dta.x3d_model_original_data.cn_vs_ep(:,iR,:,:)),3);
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,2)
iR = 45;    
atl_img_plt = mean(squeeze(sal_dta.x3d_model_original_data.cn_vs_ep(:,iR,:,:)),3);
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,3)
iR = 72;    
atl_img_plt = mean(squeeze(sal_dta.x3d_model_original_data.cn_vs_ep(:,iR,:,:)),3);
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

subplot(2,2,4)
iR = 120;    
atl_img_plt = mean(squeeze(sal_dta.x3d_model_original_data.cn_vs_ep(:,iR,:,:)),3);
img_hld = imagesc(atl_img_plt);
title([ 'Slice: ' num2str(iR)] );

%%
iR = 38;
[~, cor_ind ] = max(sum(squeeze(sum(nft_atl==iR,1)),2));

plt_ind = cor_ind-15:5:cor_ind+15;

for iP = 1:numel(plt_ind)
    
    figure()
    
    % make atlas image
    subplot(1,2,1)
    
    atl_img_plt = squeeze(nft_atl(:,plt_ind(iP),:));
    atl_img_plt(atl_img_plt==0) = 255;
    atl_img_plt(atl_img_plt==37 | atl_img_plt==38) = 254;
    atl_img_plt(atl_img_plt<254) = 175;
    atl_img_plt(atl_img_plt==254) = 100;
    img_hld = imagesc(rot90(atl_img_plt));
    title([ 'Slice: ' num2str(plt_ind(iP))] );
    
    % Show Saliency
    subplot(1,2,2)
    atl_img_plt = mean(squeeze(sal_dta.x3d_model_original_data.cn_vs_ep(:,plt_ind(iP),:,:)),3);
    imagesc(rot90(atl_img_plt),[0 0.40]);
    
    tightfig();
    set(gcf,'Position',[0 0 1920 1080])
    print(gcf, [ plt_dir '/' 'Slice_ ' num2str(plt_ind(iP))],'-dpng')
    close all
    
end