clear; clc;

%% Setup
img_dir     = '/home/mmilmcd/data/MST_stimuli/';
fld_dir     = { 'Set C' 'Set D' 'Set E' 'Set F'};
img_out_dir = '/home/ekaestne/Downloads/image_resized';
fld_out_dir = { 'SetC_resized' 'SetD_resized' 'SetE_resized' 'SetF_resized' };

trg_height = 400;
trg_width  = 400;
    
%% Resize
trg_height = trg_height-20:trg_height+20;
trg_width = trg_width-20:trg_width+20;

for iF  = 1:length(fld_dir)
    
    mkdir([ img_out_dir '/' fld_out_dir{iF} ])
    
    stm_fld = [ img_dir '/' fld_dir{iF} ];
    fle_nme = [ stm_fld '/' '*.jpg' ];
    all_fle = dir(fle_nme);
        all_fle = {all_fle(:).name};
    
    for iI = 1:length(all_fle)
               
        img_fle_nme = [ stm_fld '/' all_fle{iI} ];
        
        img_hld = imread(img_fle_nme);
%         if ~any( size(img_hld,1)==trg_width ) && ~any( size(img_hld,2)==trg_height ); error(); ; end
        
        if any( size(img_hld,1)==trg_width ) || any( size(img_hld,2)==trg_height )
            imwrite( img_hld, [ img_out_dir '/' fld_out_dir{iF} '/' all_fle{iI}(1:end-4)  '_400x400_' fld_out_dir{iF}  all_fle{iI}(end-3:end)] )
        else
            
            act_sze = size(img_hld);
            [ max_sze, max_ind] = max(act_sze);
            
            scl_fct = max_sze / mean(trg_height);
            img_rsz = imresize( img_hld, [ round(act_sze(1)/scl_fct) round(act_sze(2)/scl_fct) ]);
            
            imwrite( img_hld, [ img_out_dir '/' fld_out_dir{iF} '/' all_fle{iI}(1:end-4)  '_400x400_' fld_out_dir{iF}  all_fle{iI}(end-3:end)] )
        end
        
    end
end