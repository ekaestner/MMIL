clear; clc;

%%
prj_dir = '/home/ekaestne/PROJECTS/';
sbj_grp = {'fc063'};

fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = sbj_grp(:,1);
fcfg.dta_typ = 'dti_cnn'; % Alicia fMRI Number of Voxels
dti_cnn_dta = ejk_load_mcd_data(fcfg);

low_tri = logical(tril(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),-1));
lbl_hld = dti_cnn_dta.cll_lbl(logical(low_tri));
lbl_hld = cellfun(@(x) strsplit(x,'=='),lbl_hld,'uni',0);

% %%%%%%%%%%%%%%%%%%%%%%%%
tmp_roi_nme = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' };

bad_roi = { 'unknown' 'middletemporal' 'postcentral' 'precentral' 'superiorfrontal' 'superiortemporal' 'corpuscallosum' 'fusiform' 'inferiortemporal' 'rostralmiddlefrontal' };

[~,lhs_lbl,~]=fs_read_annotation(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' 'lh' '.aparc.split.annot']);
    lhs_lbl = lhs_lbl(1:59);
    [~,rmv_roi] = intersect(lhs_lbl,bad_roi);
    lhs_lbl(rmv_roi) = [];
    [~,tmp_roi] = intersect(lhs_lbl,tmp_roi_nme);
    lhs_tmp_lbl = lhs_lbl(tmp_roi);
    lhs_rst_lbl = setxor(lhs_lbl,lhs_tmp_lbl);

roi_lbl = [ strcat('lh.',lhs_tmp_lbl) ;  strcat('lh.',lhs_rst_lbl) ;  strcat('rh.',lhs_tmp_lbl) ; strcat('rh.',lhs_rst_lbl) ];
    
%% Load PCA
pca_wgh = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/pca/pca_weights.csv');
pca_nme = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/pca/column_names.csv');

%% Connectome Box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iPC = 1:10
    
    jet_map = jet(10000);
    jet_ind = linspace( -0.2 , 0.2 , 10000);
    
    for iR = 1:98
        for iC = 1:98
            
            if iC < iR
                patch( [ iC-1 iC-1 iC iC ] , [ iR-1 iR iR iR-1 ] , rgb('grey') )
            elseif iC == iR
                patch( [ iC-1 iC-1 iC iC ] , [ iR-1 iR iR iR-1 ] , rgb('white') )
            else
                row = string_find( tmp_roi_nme , {roi_lbl{iR}(4:end)} );
                ind = intersect( string_find( pca_nme , roi_lbl(iR) ) , string_find( pca_nme , roi_lbl(iC) ) );
                if ~isempty(row) && ~isempty(ind)
                    [ ~ , col_ind ] = min( abs(jet_ind - pca_wgh{iPC+1,ind}) );
                    patch( [ iC-1 iC-1 iC iC ] , [ iR-1 iR iR iR-1 ] , jet_map(col_ind,:) )
                else
                    patch( [ iC-1 iC-1 iC iC ] , [ iR-1 iR iR iR-1 ] , rgb('grey') )
                end
            end
            
        end
    end
    
    axis('off')
    
    print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/memory/pca/' '/' 'ConnectomePCABox' '_' num2str(iPC) '.png'],'-dpng','-r200')
    close all
    
end

    
    
    
    
    
    
    
    
    
    
    
    
    
%% Connectome Box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_row = 98:-1:98-numel(tmp_roi)+1;
    tmp_row_ind = 1;
rst_row = 49-numel(tmp_roi):-1:1;
    rst_row_ind = 1;
ovr_ind = 1;   
    
row_num = fliplr(1:49);
col_num = fliplr(1:59);
for iR = 1:59
    
    if ~any(strcmpi( bad_roi , pct_hld{iR,1} ))
        for iC = 1:49
            if any( strcmpi( tmp_roi , pct_hld{iR,1} ))
                if tmp_row(tmp_row_ind) == row_num(iC)
                    patch( [ iC-1 iC-1 iC iC ] , [ tmp_row(tmp_row_ind)-1 tmp_row(tmp_row_ind) tmp_row(tmp_row_ind) tmp_row(tmp_row_ind)-1 ]  , rgb('white'))
                elseif tmp_row(tmp_row_ind) < row_num(iC)
                    patch( [ iC-1 iC-1 iC iC ] , [ tmp_row(tmp_row_ind)-1 tmp_row(tmp_row_ind) tmp_row(tmp_row_ind) tmp_row(tmp_row_ind)-1 ] , rgb('grey') )
                else
                    patch( [ iC-1 iC-1 iC iC ] , [ tmp_row(tmp_row_ind)-1 tmp_row(tmp_row_ind) tmp_row(tmp_row_ind) tmp_row(tmp_row_ind)-1 ] , col{iR+1} )
                end
            elseif ~any(strcmpi( tmp_roi , pct_hld{iR,1} )) && ~any(strcmpi( bad_roi , pct_hld{iR,1} ))
                if rst_row(rst_row_ind) == row_num(iC)
                    patch( [ iC-1 iC-1 iC iC ] , [ rst_row(rst_row_ind)-1 rst_row(rst_row_ind) rst_row(rst_row_ind) rst_row(rst_row_ind)-1 ]  , rgb('white'))
                elseif rst_row(rst_row_ind) < row_num(iC)
                    patch( [ iC-1 iC-1 iC iC ] , [ rst_row(rst_row_ind)-1 rst_row(rst_row_ind) rst_row(rst_row_ind) rst_row(rst_row_ind)-1 ]  , rgb('grey'))
                else
                    patch( [ iC-1 iC-1 iC iC ] , [ rst_row(rst_row_ind)-1 rst_row(rst_row_ind) rst_row(rst_row_ind) rst_row(rst_row_ind)-1 ] , col{iR+1} )
                end
            end
        end        
        if any( strcmpi( tmp_roi , pct_hld{iR,1} )); tmp_row_ind = tmp_row_ind + 1; end
        if ~any(strcmpi( tmp_roi , pct_hld{iR,1} )) && ~any(strcmpi( bad_roi , pct_hld{iR,1} )); rst_row_ind = rst_row_ind + 1; end
        ovr_ind = ovr_ind + 1;
    end
    
end

axis('off')

print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/ROI/' '/' 'ConnectomeBox' '.png'],'-dpng','-r200')
close all

