clc; clear;

% WHITE MATTER TRACTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];

prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'PhenotypeConnectome';

fcfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'tracts' ];
fcfg.out_nme = 'test';

fcfg.trc_nme = { 'L_IFO' 'R_IFO' 'L_ILF' 'R_ILF' 'L_UNC' 'R_UNC' 'L_tSLF' 'R_tSLF' };

ejk_tract_plot(fcfg)

% ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_roi = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' };

sph = { 'lh' 'rh' };
sph_vew = { 'lat' 'ven' 'med' };

[~,albl,~]=fs_read_annotation(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' sph{1} '.aparc.split.annot']);
pct_hld = [albl num2cell(zeros(size(albl,1),1))];

ind = 1:59;
pct_hld = pct_hld(ind,:);
for iR = ind
    pct_hld(ind(iR),2)     = { iR };
end

srf_brn{1} = fs_read_surf(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'surf' '/' 'lh.pial']);
srf_brn{1}.surf_brain.coords = srf_brn{1}.vertices;
srf_brn{1}.surf_brain.faces = srf_brn{1}.faces;
srf_brn{2} = fs_read_surf(['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'surf' '/' 'rh.pial']);
srf_brn{2}.surf_brain.coords = srf_brn{2}.vertices;
srf_brn{2}.surf_brain.faces = srf_brn{2}.faces;

num_col = 125;
dst_col = distinguishable_colors(num_col);
for iDC = 1:num_col
    dst_mtx(iDC,1) = pdist2( rgb( 'blue' ) , dst_col(iDC,:) );
end
[ ~ , dst_mtx_srt ] = sort(dst_mtx);
tmp_col = dst_mtx_srt(1:numel(tmp_roi)*1);
rst_col = dst_mtx_srt( end - (59*1) :end );

tmp_ind = 1;
rst_ind = 1;
col{1} = rgb('grey');
for iR = 1:size(pct_hld,1)
    if any(strcmpi( tmp_roi , pct_hld{iR,1} ))
        col{iR + 1} = dst_col(tmp_col(tmp_ind),:);
        col{iR + 1} = col{iR + 1} + [ 0.00 0.05 0.125 ];
        if col{iR + 1}(1) > 1; col{iR + 1}(1) = 1; end
        if col{iR + 1}(2) > 1; col{iR + 1}(2) = 1; end
        if col{iR + 1}(3) > 1; col{iR + 1}(3) = 1; end
        tmp_ind = tmp_ind+1;
    else
        col{iR + 1} = dst_col(rst_col(rst_ind),:);
        col{iR + 1} = col{iR + 1} - [ 0.20 0.20 0.20 ];
        if col{iR + 1}(1) < 0; col{iR + 1}(1) = 0; end
        if col{iR + 1}(2) < 0; col{iR + 1}(2) = 0; end
        if col{iR + 1}(3) < 0; col{iR + 1}(3) = 0; end
        rst_ind = rst_ind+1;
    end
end

col_map = [];
top_pct = 1;
for iC = 1:numel(col)-1
    col_map = [col_map ; [linspace(col{iC}(1),col{iC+1}(1),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(2),col{iC+1}(2),ceil(1000*top_pct/(numel(col)-1)))' linspace(col{iC}(3),col{iC+1}(3),ceil(1000*top_pct/(numel(col)-1)))']; ];
end
top_pct = numel(col)-1;

for iH = 1:2
    for iSP = 1:numel(sph_vew)
        
        fig_hld(1) = figure('Visible','off','Position',[0 0 1920 1080]);
        
        pcfg = [];
        
        pcfg.surf_brain  = srf_brn{iH};
        pcfg.aparc       = ['/space/syn09/1/data/MMILDB/ABALA/Connectome/preproc' '/' 'fsaverage' '/' 'label' '/' sph{iH} '.aparc.split.annot'];
        
        pcfg.sph         = sph{iH};
        pcfg.sph_vew     = sph_vew{iSP};
        
        pcfg.label       = 0;
        pcfg.radius      = [];
        pcfg.alpha       = 1;
        
        pcfg.non_ele     = [];
        pcfg.sve_img     = 0; % ###
        
        if ~strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 1],'visible','off','Parent',fig_hld(1));
        elseif strcmpi(sph_vew{iSP},'ven')
            pcfg.axe_hnd = axes('OuterPosition',[0 0 1 0.6],'visible','off','Parent',fig_hld(1));
        end
        
        pcfg.fig_hdl = fig_hld(1);
        
        pcfg.col_map = col_map;
        
        pcfg.tbl_pct = cell2mat(pct_hld(:,2)) ./ top_pct;
        pcfg.tbl_pct(isnan(pcfg.tbl_pct)) = 0;
        
        pcfg.tbl_loc = strcat('lhs_',pct_hld(:,1));
        
        pcfg.top_pct = 1;
        
        nyu_plot2(pcfg);
        
        print(['/home/ekaestne/PROJECTS/OUTPUT/PhenotypeConnectome/ROI/' '/' 'ROI_guide' '_' sph{iH} '_' sph_vew{iSP} '.png'],'-dpng','-r200')
        close all
        
    end
end

% Connectome Box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad_roi = { 'unknown' 'middletemporal' 'postcentral' 'precentral' 'superiorfrontal' 'superiortemporal' 'corpuscallosum' 'fusiform' 'inferiortemporal' 'rostralmiddlefrontal' };

tmp_row = 49:-1:49-numel(tmp_roi)+1;
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

 

% Connectome %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.sbj_nme = sbj_grp(:,1);
fcfg.dta_typ = 'dti_cnn'; % Alicia fMRI Number of Voxels
dti_cnn_dta = ejk_load_mcd_data(fcfg);

low_tri = logical(tril(ones(size(dti_cnn_dta.dti_cnn_dta,2),size(dti_cnn_dta.dti_cnn_dta,3)),-1));
lbl_hld = dti_cnn_dta.cll_lbl(logical(low_tri));

tmp_roi = { 'caudal-ITG'   'middle-ITG'         'rostral-ITG' ...
            'caudal-MTG'   'middle-MTG'         'rostral-MTG' ...
            'caudal-STG'   'middle-STG'         'rostral-STG' ...
            'temporalpole' 'transversetemporal' 'bankssts' ...
            'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' ...
            'parahippocampal' 'entorhinal' };

lbl_spl_hld = regexp(lbl_hld,'==','split');
hld_nme = [];
for iTR = 1:numel(tmp_roi)
    hld_nme = [hld_nme find(~cellfun(@isempty,cellfun(@(x) string_find(x,tmp_roi{iTR}),lbl_spl_hld,'uni',0)))];
end
hld_nme = unique(hld_nme);
        
fcfg = [];

fcfg.prj_dir = prj_dir;
fcfg.prj_nme = prj_nme;
fcfg.out_nme = [ 'ROI/ROI_connectome.png' ];

fcfg.brn_nme = 'fsaverage';
fcfg.nde_nme = 'FSAverageDesikanMod_master';

fcfg.edg_nme = { lbl_hld(hld_nme)       };
fcfg.edg_wgh = { ones(numel(hld_nme),1)*0.5 };
fcfg.edg_col = { rgb('red')                       };
fcfg.edg_cut_off = 0;

connectome_plot(fcfg)