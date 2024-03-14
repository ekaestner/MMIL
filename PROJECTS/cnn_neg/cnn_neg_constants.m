dta_fld = '/home/ekaestner/Dropbox/McDonald Lab/Erik/02_Projects/02_Core/cnn_neg';

%%
out_fld = [ dta_fld '/' 'Out' ];

nii_dir = [ dta_fld '/' 'Data' '/' 'Saliency' '/' ];

plt_fld = [ out_fld '/' 'plots' '/'];

atl_dir = [ nii_dir '/' 'atlas'];

%%
sal_grp = { 'MRI_neg_SF' 'MRI_pos_SF' };
atl_nme = 'aal';