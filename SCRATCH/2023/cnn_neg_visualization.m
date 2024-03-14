clear; clc;

dta_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/05_Developement/cnn_neg';
nii_dir = [ dta_fld '/' 'Data' '/' 'Saliency' '/' ];
plt_fld = [ dta_fld '/' 'Out' '/' 'plots' '/'];

sal_grp = { 'MRI_neg_SF' 'MRI_pos_SF' };

%% Load
% Get dummy size
nft_fle = dir([ nii_dir '/' sal_grp{iSA} '/' '*.nii']); nft_fle = {nft_fle(1).name};
dum_dta = niftiread([ nii_dir '/' sal_grp{1} '/' nft_fle{1}]);
clear nft_fle

% Load
for iSA = 1:numel(sal_grp)
    
    nft_fle.(sal_grp{iSA}) = dir([ nii_dir '/' sal_grp{iSA} '/' '*.nii']); nft_fle.(sal_grp{iSA}) = {nft_fle.(sal_grp{iSA})(:).name};    
    fle_nme.(['x' sal_grp{iSA}]) = nft_fle.(sal_grp{iSA});
    
    sal_dta.(['x' sal_grp{iSA}]) = nan( [ size(dum_dta) numel( fle_nme.(['x' sal_grp{iSA}])) ] );
    
    for iF = 1:numel(fle_nme.(['x' sal_grp{iSA}]))
        sal_dta.(['x' sal_grp{iSA}])(:,:,:,iF) = niftiread([ nii_dir '/' sal_grp{iSA} '/' fle_nme.(['x' sal_grp{iSA}]){iF} ]);
        
        % Remove 0's
        sal_dta.(['x' sal_grp{iSA}])(sal_dta.(['x' sal_grp{iSA}])==0) = NaN;
        
    end
end

%% Normalize
low_pct = 75;

% Put together data structure
for iZ = 1:numel(sal_grp)
    
    sal_dta_nor.all.keep.(['x' sal_grp{iZ}]) = sal_dta.(['x' sal_grp{iZ}]);
    
    for iF = 1:size(sal_dta_nor.all.keep.(['x' sal_grp{iZ}]),4)
        
        % Values to play with
        sal_val_hld = squeeze(sal_dta_nor.all.keep.(['x' sal_grp{iZ}])(:,:,:,iF));
        
        % Normalize
        sal_val_men = nanmean(sal_val_hld(:));
        sal_val_std = nanstd(sal_val_hld(:));
        sal_val_hld = (sal_val_hld-sal_val_men) ./ sal_val_std;
        
        % Remove values
        sal_val_hld_rmv     = sal_val_hld;
        sal_val_hld_rmv_mid = sal_val_hld;
        
        % Remove bottom
        pct_hld = prctile(sal_val_hld(:),low_pct);
        sal_val_hld_rmv(sal_val_hld<pct_hld) = NaN;
        
        % Re-assign
        sal_dta_nor.zsc.keep.(['x' sal_grp{iZ}])(:,:,:,iF)       = sal_val_hld;
        sal_dta_nor.zsc.remove.(['x' sal_grp{iZ}])(:,:,:,iF)     = sal_val_hld_rmv;
        
        clear sal_val_hld sal_val_hld_rmv
        
    end
end

%% Average Axial plots
ejk_chk_dir([ plt_fld '/'  'axial_average' '/' ]);

cmp_nme = { 'xMRI_neg_SF'       'xMRI_pos_SF' };
nor_nme = { { 'all' 'zsc' } { 'all' 'zsc' } };
kep_nme = { { 'keep' 'keep' }   { 'keep' 'keep' } }; 
col_nme = { 'hot'   'hot'  };
grp_lim = { {[-.002 .002] [0 1.5] } {[-.002 .002] [0 1.5] } };
hot_col = hot;

for iP = 1:numel(cmp_nme)
    for iPP = 1:numel(nor_nme{iP})
        
        slc = 18:5:size(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),3)-18;
        dim = ceil(sqrt(numel(slc)));
        
        figure('Visible','off');
        
        for iSB = 1:numel(slc)
            
            one_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP})(11:end-9,11:end-8,slc(iSB),:),4)));
            
            imAlpha_one=ones(size(one_hld));
            imAlpha_one(isnan(one_hld))=0;
            
            ax(iSB) = subplot(dim,dim,iSB);
            img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ grp_lim{iP}{iPP}(1) grp_lim{iP}{iPP}(2) ]);
            colormap(ax(iSB),hot_col(1:end-45,:)); %col_nme{iP}
            set(gca,'color',rgb('white'));
            axis off;
            set(gcf, 'InvertHardcopy', 'off')
            
        end
        
        tightfig();
        set(gcf','Position',[0 0 1080 1080]);
        print( [ plt_fld '/'  'axial_average' '/' cmp_nme{iP} '_' nor_nme{iP}{iPP} '_' kep_nme{iP}{iPP}  '_' col_nme{iP} '.png'] ,'-dpng')
        close all
        
    end
end

%% Average Coronal plots
ejk_chk_dir([ plt_fld '/'  'coronal_average' '/' ]);

cmp_nme = { 'xMRI_neg_SF'       'xMRI_pos_SF' };
nor_nme = { { 'all' 'zsc' } { 'all' 'zsc' } };
kep_nme = { { 'keep' 'keep' }   { 'keep' 'keep' } }; 
col_nme = { 'hot'   'hot'  };
grp_lim = { {[-.002 .002] [0 1.5] } {[-.002 .002] [0 1.5] } };
hot_col = hot;

for iP = 1:numel(cmp_nme)
    for iPP = 1:numel(nor_nme{iP})
        
        slc = 22:6:size(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),2)-22;
        dim = ceil(sqrt(numel(slc)));
        
        figure('Visible','off');
        
        for iSB = 1:numel(slc)
            
            one_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP})(11:end-9,slc(iSB),11:end-8,:),4)));
            
            imAlpha_one=ones(size(one_hld));
            imAlpha_one(isnan(one_hld))=0;
            
            ax(iSB) = subplot(dim,dim,iSB);
            img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ grp_lim{iP}{iPP}(1) grp_lim{iP}{iPP}(2) ]);
            colormap(ax(iSB),hot_col(1:end-45,:)); %col_nme{iP}
            set(gca,'color',rgb('white'));
            axis off;
            set(gcf, 'InvertHardcopy', 'off')
            
        end
        
        tightfig();
        set(gcf','Position',[0 0 1080 1080]);
        print( [ plt_fld '/'  'coronal_average' '/' cmp_nme{iP} '_' nor_nme{iP}{iPP} '_' kep_nme{iP}{iPP}  '_' col_nme{iP} '.png'] ,'-dpng')
        close all
        
    end
end

%% Single subject visualization
ejk_chk_dir([ plt_fld '/'  'single_subjects' '/' cmp_nme{iP} '/']);

cor_slc = 76;

cmp_nme = { 'xMRI_neg_SF'       'xMRI_pos_SF' };
nor_nme = { { 'zsc' }    { 'zsc' } };
kep_nme = { { 'keep' }   { 'keep' } }; 
col_nme = { 'hot'   'hot'  };
grp_lim = { { [0 1.5] } { [0 1.5] } };
hot_col = hot;

for iP = 1:numel(cmp_nme)
    ejk_chk_dir([ plt_fld '/'  'single_subjects' '/' cmp_nme{iP} '/']);
    for iPP = 1:numel(nor_nme{iP})`
        
        axl_slc = 18:5:size(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),3)-18;
            axl_slc = axl_slc(1:8);
        
        cor_slc = 22:6:size(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),2)-22;
            cor_slc = cor_slc(5:12);
            
            for iS = 1:size(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),4)
            
                figure('Visible','off');
                
                % Axial
                for iSB = 1:numel(axl_slc)
                    
                    one_hld = rot90(squeeze(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP})(11:end-9,11:end-8,axl_slc(iSB),iS)));
                    
                    imAlpha_one=ones(size(one_hld));
                    imAlpha_one(isnan(one_hld))=0;
                    
                    ax(iSB) = subplot(dim,dim,iSB);
                    img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ grp_lim{iP}{iPP}(1) grp_lim{iP}{iPP}(2) ]);
                    colormap(ax(iSB),hot_col(1:end-45,:)); %col_nme{iP}
                    set(gca,'color',rgb('white'));
                    axis off;
                    set(gcf, 'InvertHardcopy', 'off')
                    
                end
                
                % Coronal
                for iSB = 1:numel(cor_slc)
                    
                    one_hld = rot90(squeeze(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP})(11:end-9,cor_slc(iSB),11:end-8,iS)));
                    
                    imAlpha_one=ones(size(one_hld));
                    imAlpha_one(isnan(one_hld))=0;
                    
                    ax(iSB) = subplot(dim,dim,iSB+8);
                    img_one = imagesc(one_hld,'AlphaData',imAlpha_one,[ grp_lim{iP}{iPP}(1) grp_lim{iP}{iPP}(2) ]);
                    colormap(ax(iSB),hot_col(1:end-45,:)); %col_nme{iP}
                    set(gca,'color',rgb('white'));
                    axis off;
                    set(gcf, 'InvertHardcopy', 'off')
                    
                end
                
                % Save out
                tightfig();
                set(gcf','Position',[0 0 1080 1080]);
                print( [ plt_fld '/'  'single_subjects' '/' cmp_nme{iP} '/' nft_fle.(cmp_nme{iP}(2:end)){iS}(1:end-4) '_' cmp_nme{iP} '_' nor_nme{iP}{iPP} '_' kep_nme{iP}{iPP}  '_' col_nme{iP} '.png'] ,'-dpng')
                close all
                
            end
        
    end
end





