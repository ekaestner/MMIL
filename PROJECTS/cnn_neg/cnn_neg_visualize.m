%% Save average nii
cmp_nme = { 'xMRI_neg_SF'       'xMRI_pos_SF' };
nor_nme = { { 'all' 'zsc' } { 'all' 'zsc' } };
kep_nme = { { 'keep' 'keep' }   { 'keep' 'keep' } }; 

for iP = 1:numel(cmp_nme)
    for iPP = 1:numel(nor_nme{iP})
        
        one_hld = squeeze(nanmean(sal_dta_nor.(nor_nme{iP}{iPP}).(kep_nme{iP}{iPP}).(cmp_nme{iP}),4));
        
         niftiwrite( one_hld, [ nii_dir '/'  'average' '_' 'saliency' '_' 'nii' '_' cmp_nme{iP} '_' nor_nme{iP}{iPP} '.nii' ] );
        
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

%% Difference Average Axial plots
ejk_chk_dir([ plt_fld '/'  'axial_average' '/' ]);

cmp_nme = { 'xMRI_neg_SF'       'xMRI_pos_SF' };
nor_nme = { { 'all' 'zsc' }    };
kep_nme = { { 'keep' 'keep' }  };
col_nme = { 'hot'   };
grp_lim = {[-.002 .002] [-0.05 0.05] };
hot_col = jet; hot_col(70:186,:) = repmat(rgb('grey'),size(hot_col(70:186,:),1),1);

%for iP = 1:numel(cmp_nme)
for iPP = 1:numel(nor_nme{1})
    
    slc = 18:5:size(sal_dta_nor.(nor_nme{1}{iPP}).(kep_nme{1}{iPP}).(cmp_nme{1}),3)-18;
    dim = ceil(sqrt(numel(slc)));
    
    figure('Visible','off');
    
    for iSB = 1:numel(slc)
        
        one_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{1}{iPP}).(kep_nme{1}{iPP}).(cmp_nme{1})(11:end-9,11:end-8,slc(iSB),:),4)));
        two_hld = rot90(squeeze(nanmean(sal_dta_nor.(nor_nme{1}{iPP}).(kep_nme{1}{iPP}).(cmp_nme{2})(11:end-9,11:end-8,slc(iSB),:),4)));
        dff_hld = two_hld - one_hld;
        
        imAlpha_one=ones(size(dff_hld));
        imAlpha_one(isnan(dff_hld))=0;
        
        ax(iSB) = subplot(dim,dim,iSB);
        img_one = imagesc(dff_hld,'AlphaData',imAlpha_one,[ grp_lim{1}(1) grp_lim{1}(2) ]);
        colormap(ax(iSB),hot_col); %col_nme{iP}
        set(gca,'color',rgb('white'));
        axis off;
        set(gcf, 'InvertHardcopy', 'off')
        
    end
    
    tightfig();
    set(gcf','Position',[0 0 1080 1080]);
    print( [ plt_fld '/'  'axial_diff_average' '/' 'diff' '_' cmp_nme{2} '_' 'minus' '_' cmp_nme{1} '_' nor_nme{1}{iPP} '_' kep_nme{1}{iPP}  '_' col_nme{1} '.png'] ,'-dpng')
    close all
    
end
% end

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
    for iPP = 1:numel(nor_nme{iP})
        
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
                set(gcf,'Position',[0 0 1080 1080]);
                print( [ plt_fld '/'  'single_subjects' '/' cmp_nme{iP} '/' nft_fle.(cmp_nme{iP}(2:end)){iS}(1:end-4) '_' cmp_nme{iP} '_' nor_nme{iP}{iPP} '_' kep_nme{iP}{iPP}  '_' col_nme{iP} '.png'] ,'-dpng')
                close all
                
            end
        
    end
end
