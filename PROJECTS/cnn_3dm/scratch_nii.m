clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/cnn_3dm/';

nii_dir = [ prj_dir '/' 'Data/2023_02_24/saliency_maps/' ]; % 2022_12_12

cat_nii     = { '^cn-ep_'     '^cn_'            '^ep_' };
cat_nii_nme = { 'cn_vs_ep'    'cn'              'ep' };
cat_col     = { rgb('hazel')  rgb('light blue') rgb('light red') };
cat_nii_cmp = {[2 3] [1]};

atl_dir = [ prj_dir '/' 'Data/atlas/aal_bonilha'];
atl_nme = 'aal';

plt_dir = [ prj_dir '/' 'Data/2023_02_24/saliency_plots/' ];
    ejk_chk_dir( [ plt_dir '/' 'ROIs' '/' ]);
    ejk_chk_dir( [ plt_dir '/' 'CoronalSlices' '/' ]);
        
%% Load Atlas
nft_atl = niftiread([ atl_dir '/' atl_nme '.nii' ]);
% imagesc(rot90(squeeze(nft_atl(:,90,:))));

nft_atl_lbl = mmil_readtext([ atl_dir '/' atl_nme '.txt' ]);
    nft_atl_lbl = regexpi(nft_atl_lbl,'\|','split');
    nft_atl_lbl = cat(1,nft_atl_lbl{:});
       
%% Load Data
% Get filenames of the the individual categories
nft_fle = dir(nii_dir); nft_fle = {nft_fle(:).name}; nft_fle = nft_fle(3:end);
for iN = 1:numel(cat_nii_nme)
    fle_nme.(cat_nii_nme{iN}) = nft_fle(string_find(nft_fle,cat_nii{iN}));
end

% reslice
flg = [];
flg.prefix = 'aal_reslice_';
flg.interp = 0;
for iN = 1:numel(cat_nii_nme)
    spm_reslice( [[ atl_dir '/' atl_nme '.nii' ] ; strcat(nii_dir,fle_nme.(cat_nii_nme{iN}))'],flg);
end

flg = [];
flg.prefix = 'map_reslice_';
flg.interp = 0;
spm_reslice( [ strcat(nii_dir,fle_nme.(cat_nii_nme{1})(1)); [ atl_dir '/' atl_nme '.nii' ]],flg);

% Create ROI data structure
for iR = 1:size(nft_atl_lbl,1)
    dta_hld.(nft_atl_lbl{iR,3}) = nan( numel(fle_nme.(cat_nii_nme{iN})), numel(cat_nii) );
end

% Load through: 1) Load data,
for iN = 1:numel(cat_nii_nme)
    cor_slc_hld = nan([size(nft_atl) numel(fle_nme.(cat_nii_nme{iN}))]);
    
    for iF = 1:numel(fle_nme.(cat_nii_nme{iN}))
        
        % Load Data
        nft_dta = niftiread([ nii_dir '/' 'aal_reslice_' fle_nme.(cat_nii_nme{iN}){iF} ]);
        nft_dta(nft_dta==0) = NaN;
        cor_slc_hld(:,:,:,iF) = nft_dta;
        
        % Find average ROI saliency
        for iR = 1:size(nft_atl_lbl,1)
            %
            dta_hld.(nft_atl_lbl{iR,3})(iF,iN) = nanmean( nft_dta(nft_atl==iR));            
        end
    end
    
    cor_slc.(cat_nii_nme{iN}) = nanmean(cor_slc_hld,4);
    
end

%% Statistical Tests



%% Plot data
for iR = 1:size(nft_atl_lbl,1)
   
    figure('Visible','off')
    
    for iP = 1:numel(cat_nii_cmp)
        sbp = subplot(1,numel(cat_nii_cmp),iP);
        
        fcfg = [];        

        for iD = 1:numel(cat_nii_cmp{iP})
            fcfg.ydt{iD} = dta_hld.(nft_atl_lbl{iR,3})(:,cat_nii_cmp{iP}(iD));
            fcfg.xdt{iD} = iD;
        end
        
        fcfg.fce_col     = cat_col(cat_nii_cmp{iP});
        fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(cat_nii_cmp{iP})) ];
        fcfg.box_plt_col = cat_col(cat_nii_cmp{iP});
        
        fcfg.box_plt = ones(1,numel(fcfg.xdt));
        fcfg.xlb = cat_nii_nme(cat_nii_cmp{iP});
        fcfg.xlm = [ 0.5 numel(cat_nii_cmp{iP})+2.5 ];
        fcfg.ylb = cat_nii_nme(cat_nii_cmp{iP});
        
        fcfg.ttl = nft_atl_lbl{iR,3};
        
        fcfg.mkr_sze = repmat(40,1,numel(fcfg.xdt));
        fcfg.aph_val = 0.80;
        
        fcfg.sbp = sbp;
        
        ejk_scatter(fcfg)        
    end
    
    print( [ plt_dir '/' 'ROIs' '/' nft_atl_lbl{iR,3} '.png'] ,'-dpng')
    close all
    
end

%% Plot saliency maps
for iR = 1:size(nft_atl_lbl,1)
    
    figure('Visible','off')
    
    % make atlas image
    subplot(2,2,1)
    [~, cor_ind ] = max(sum(squeeze(sum(nft_atl==iR,1)),2));
    atl_img_plt = squeeze(nft_atl(:,cor_ind,:));
    atl_img_plt(atl_img_plt==0) = 255;
    atl_img_plt(atl_img_plt==iR) = 254;
    atl_img_plt(atl_img_plt<254) = 175;
    atl_img_plt(atl_img_plt==254) = 100;
    img_hld = imagesc(rot90(atl_img_plt));
    title(mmil_spec_char(nft_atl_lbl{iR,3},{'_'},{' '}))
    
    % make rest of images
    for iN = 1:numel(cat_nii_nme)
        subplot(2,2,iN+1)
        
        atl_img_plt = squeeze(cor_slc.(cat_nii_nme{iN})(:,cor_ind,:));
        imagesc(rot90(atl_img_plt),[0 0.60]);
        title(mmil_spec_char(cat_nii_nme{iN},{'_'},{' '}))
    end
   
    tightfig();
    print( [ plt_dir '/' 'CoronalSlices' '/' nft_atl_lbl{iR,3} '.png'] ,'-dpng')
    close all
    
end

%% Plot L/R
roi_run = [1 3];

for iRR = 1:numel(roi_run)
    
    out_dir_roi = [ plt_dir '/' 'ROIs' '_' cat_nii_nme{roi_run(iRR)} '/' ];
    ejk_chk_dir(out_dir_roi)
    
    % ROI Performance %%%%%%%%%%%%%%%%
    prf_hld = nan(size(nft_atl_lbl,1),1);
    for iR = 1:size(nft_atl_lbl,1)
        prf_hld(iR,1) = nanmean(dta_hld.(nft_atl_lbl{iR,3})(:,roi_run(iRR)));
    end
    
    % Plot %%%%%%%%%%%%%%%%
    cut_off = 0:0.05:1; cut_off_hld = cut_off>max(prf_hld); cut_off_hld(find(cut_off_hld,1)) = 0; cut_off(cut_off_hld) = [];
    
    for iC = 1:numel(cut_off)-1
        
        roi_use = nft_atl_lbl(prf_hld>=cut_off(iC) & prf_hld<cut_off(iC+1),3);
        prf_use = prf_hld(prf_hld>=cut_off(iC) & prf_hld<cut_off(iC+1));
        [~, srt_ind] = sort(prf_use);
        roi_use = roi_use(srt_ind);
        
        if ~isempty(roi_use)
            fcfg = [];
            
            for iR = 1:numel(roi_use)
                fcfg.ydt{iR} = dta_hld.(roi_use{iR})(:,roi_run(iRR));
                fcfg.xdt{iR} = iR;
            end
            
            fcfg.fce_col     = [ repmat({[0.2 0.2 0.2]},1,numel(fcfg.xdt)) ];
            fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
            fcfg.box_plt_col = [ repmat({[0.2 0.2 0.2]},1,numel(fcfg.xdt)) ];
            
            fcfg.box_plt = ones(1,numel(fcfg.xdt));
            fcfg.xlb = roi_use;
            fcfg.xlm = [ 0.5 numel(fcfg.xdt)+0.5 ];
            fcfg.ylb = {'Saliency'};
            
            fcfg.ylm = [.0 .30];
            
            fcfg.mkr_sze = repmat(15,1,numel(fcfg.xdt));
            fcfg.aph_val = 0.80;
            
            fcfg.out_dir = out_dir_roi;
            fcfg.out_nme = [ '00_combine_plot' '_' num2str(cut_off(iC)) '_to_' num2str(cut_off(iC+1))];
            
            ejk_scatter(fcfg)
        end
    end
    
    cell2csv([ out_dir_roi '/' 'Saliency_Values' '_' cat_nii_nme{roi_run(iRR)} '.csv'],[ nft_atl_lbl(:,3) num2cell(prf_hld)]);
    
end

%% Hold
% % Left
% lbl_hld.nft_atl_lbl_lft = nft_atl_lbl(string_find(nft_atl_lbl(:,3),'_L$'),3);
%     lbl_hld.nft_atl_lbl_lft(string_find(lbl_hld.nft_atl_lbl_lft,'Cerebelum')) = [];
%     
% % Right
% lbl_hld.nft_atl_lbl_rgh = nft_atl_lbl(string_find(nft_atl_lbl(:,3),'_R$'),3);
%     lbl_hld.nft_atl_lbl_rgh(string_find(lbl_hld.nft_atl_lbl_rgh,'Cerebelum')) = [];    
% 

%% reslice investigate
% nft_atl_org = niftiread([ atl_dir '/' 'aal.nii' ]);
% nft_atl_avg = niftiread([ atl_dir '/' 'meanaal.nii' ]);
% nft_atl_rsl = niftiread([ atl_dir '/' 'aal_reslice_aal.nii' ]);
% nft_atl_map_rsl = niftiread([ atl_dir '/' 'map_reslice_aal.nii' ]);
% nft_dta_org = niftiread([ nii_dir '/' 'cn-ep_vanilla_backprop_r0.nii' ]);
% nft_dta_rsl = niftiread([ nii_dir '/' 'aal_reslice_cn-ep_vanilla_backprop_r0.nii' ]);
% 
% figure()
% subplot(2,2,1); imagesc(rot90(squeeze(nft_atl_map_rsl(:,69,:)))); title('AAL: Map reslice');
% subplot(2,2,2); imagesc(rot90(squeeze(nft_dta_org(:,69,:)))); title('Map: Original');
% subplot(2,2,3); imagesc(rot90(squeeze(nft_atl_org(:,95,:)))); title('AAL: original');
% subplot(2,2,4); imagesc(rot90(squeeze(nft_dta_rsl(:,95,:)))); title('Map: AAL reslice');
% colormap('jet')
% set(gcf,'Position',[0 0 1280 1080]);
% 
% 
% figure()
% iR = 37;
% subplot(1,2,1)
% [~, cor_ind ] = max(sum(squeeze(sum(nft_atl_org==iR,1)),2));
% atl_img_plt = squeeze(nft_atl_org(:,cor_ind,:));
% atl_img_plt(atl_img_plt==0) = 255;
% atl_img_plt(atl_img_plt==iR) = 254;
% atl_img_plt(atl_img_plt<254) = 175;
% atl_img_plt(atl_img_plt==254) = 100;
% img_hld = imagesc(rot90(atl_img_plt));
% title([mmil_spec_char(nft_atl_lbl{iR,3},{'_'},{' '}) ' : ' 'original' ])
% 
% subplot(1,2,2)
% [~, cor_ind ] = max(sum(squeeze(sum(nft_atl_map_rsl==iR,1)),2));
% atl_img_plt = squeeze(nft_atl_map_rsl(:,cor_ind,:));
% atl_img_plt(atl_img_plt==0) = 255;
% atl_img_plt(atl_img_plt==iR) = 254;
% atl_img_plt(atl_img_plt<254) = 175;
% atl_img_plt(atl_img_plt==254) = 100;
% img_hld = imagesc(rot90(atl_img_plt));
% title([mmil_spec_char(nft_atl_lbl{iR,3},{'_'},{' '}) ' : ' 'map reslice' ])

