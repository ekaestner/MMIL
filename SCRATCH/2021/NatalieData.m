clear; clc;

ntl_dta = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/NatalieData.csv');
    ntl_col = ntl_dta(1,:);
    ntl_dta = ntl_dta(2:end,:);
sbj_dta = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/emory_slah_subjects.csv');

%%
for iS = 1:size(ntl_dta,1)
    
    sbj_nme_hld = ntl_dta{iS,2};
    sbj_nme_dsh = strfind(sbj_nme_hld,'-');
    sbj_nme_und = strfind(sbj_nme_hld,'_');
    
    if ~isempty(sbj_nme_dsh)
        sbj_cph{iS,1} = sbj_nme_hld( [ 1:sbj_nme_dsh(1)-1 sbj_nme_dsh(1)+1:sbj_nme_dsh(2)-1 ] );
    else
        sbj_cph{iS,1} = sbj_nme_hld( [ 1:sbj_nme_und(1)-1 sbj_nme_und(1)+1:end ] );
    end
    
    sbj_cph{iS,2} = ntl_dta{iS,2};
    sbj_cph{iS,3} = iS;
    
    sbj_dta_row = string_find( sbj_dta(:,1), sbj_cph{iS,1} );
    
    if isempty(sbj_dta_row)
        sbj_cph{iS,4} = '';
        sbj_cph{iS,5} = nan;
    else
        sbj_cph{iS,4} = sbj_dta{ sbj_dta_row, 1 };
        sbj_cph{iS,5} = sbj_dta_row;
    end
    
end

cell2csv( '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/subject_cipher.csv', sbj_cph)

%%
sbj_cph = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/subject_cipher_add.csv');

ntl_cat_col = find( strcmpi(ntl_col, 'Natalie Category') );
dti_col     = find( strcmpi(ntl_col, 'DTI sequence') );
srg_col     = find( strcmpi(ntl_col, 'type of surgery') );
lmm_pre_col = find( strcmpi(ntl_col, 'preop LM (ss)') );
lmm_pst_col = find( strcmpi(ntl_col, 'postop LM - 6 mth') );
ray_pre_col = find( strcmpi(ntl_col, 'preop RAVLT - immediate recall zscore') );
ray_pst_col = find( strcmpi(ntl_col, 'postop RAVLT imm recall (z-score)') );

for iS = 1:size(sbj_cph,1)
    
    sbj_dta_out{iS,1} = ntl_dta{sbj_cph{iS,3},2};
    
    if ~isempty(sbj_cph{iS,5})
        sbj_dta_out{iS,2} = sbj_dta{sbj_cph{iS,5},1};
    else
        sbj_dta_out{iS,2} = '';
    end
    
    sbj_dta_out{iS,3} = ntl_dta{sbj_cph{iS,3},ntl_cat_col};
    
    sbj_dta_out{iS,4} = ntl_dta{sbj_cph{iS,3},dti_col};
    sbj_dta_out{iS,5} = ntl_dta{sbj_cph{iS,3},lmm_pre_col};
    sbj_dta_out{iS,6} = ntl_dta{sbj_cph{iS,3},lmm_pst_col};
    sbj_dta_out{iS,7} = ntl_dta{sbj_cph{iS,3},ray_pre_col};
    sbj_dta_out{iS,8} = ntl_dta{sbj_cph{iS,3},ray_pst_col};
    sbj_dta_out{iS,9} = ntl_dta{sbj_cph{iS,3},srg_col};
    
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/combo_data.csv',sbj_dta_out)
clear sbj_dta_out

%%
sbj_cph = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/subject_cipher_add.csv');

ntl_dta = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/NatalieData.csv');
    ntl_col = ntl_dta(1,:);
    ntl_dta = ntl_dta(2:end,:);

usd_out_cme = mmil_readtext('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/Emory_Slah_Explore/ucsd_preprocessing_summary.csv');

size(sbj_cph);
size(ntl_dta);
size(usd_out_cme);

out_dta = cell( size(sbj_cph,1), 3);
chk_sbj = cell( size(sbj_cph,1), 5);
for iS = 1:size(sbj_cph,1)
    
    ntl_ind = find(strcmpi( ntl_dta(:,2), sbj_cph{iS,2} ));
    usd_ind = find(strcmpi( usd_out_cme(:,1), sbj_cph{iS,4} ));
    
    out_dta{iS,1} = sbj_cph{iS,1};
    out_dta{iS,2} = ntl_dta{ntl_ind,3};
    if ~isempty(usd_ind)
        out_dta{iS,3} = usd_out_cme{usd_ind,3};
    else
        out_dta{iS,3} = '';
    end
    
    chk_sbj{iS,1} = sbj_cph{iS,1};
    chk_sbj{iS,2} = ntl_dta{ntl_ind,2};
    if ~isempty(usd_ind)
        chk_sbj{iS,3} = usd_out_cme{usd_ind,1};
    else
        chk_sbj{iS,3} = '';
    end
    
end

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/Emory_Slah_Explore/dan_spreadsheet.csv',out_dta);
cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Emory_Slah_Explore/Emory_Slah_Explore/dan_spreadsheet_subject_check.csv',chk_sbj);


























