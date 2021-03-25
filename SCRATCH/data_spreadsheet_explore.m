clear; clc;

dta_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/RSI/';

%% Get subjects
rcn_hld = mmil_readtext('/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer/mmilmcdRSI_freesurfer_recons.csv');

sbj_hld = mmil_readtext([ '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD' '/' 'Redcap_2021_02_03.csv' ]);
    sbj_hld_nme = sbj_hld(1,[1 16 48 70 80 82 100 102]);
    sbj_hld_dta = sbj_hld(2:end,[1 16 48 70 80 82 100 102]);

mcd_rsi_dti_prj = mmil_readtext( [ dta_dir '/' 'MCD_RSI' '/' 'DTI_all.csv' ] );
    mcd_rsi_dti_col = [];

mcd_rsi_flx_prj = mmil_readtext( [ dta_dir '/' 'MCD_RSI' '/' 'DTI_flex.csv' ] );
    mcd_rsi_flx_col = [];

cnn_sbj_nme = mmil_readtext([ '/home/ekaestner/gitrep/MMIL/EXTERNAL/McD' '/' 'connectomes_norm.txt' ]);
    
dta_chk = cell( size(sbj_hld_dta,1), size(sbj_hld_nme,2)+2 );
dta_chk_nme = { 'sbj_nme' 'epd_sde' 'dti_red_cap' 'dti_spd_sht' 'rsi_red_cap' 'rsi_spd_sht' 'pre_lm2' 'pst_lm2' };
for iS = 1:size(sbj_hld_dta,1) 
    
    dta_chk(iS,1:2) = sbj_hld_dta(iS,1:2);
    if isempty(dta_chk{iS,2}); dta_chk{iS,2} = 'control'; end
    
    rcn_nme_ind = find(strcmpi( rcn_hld(:,1), dta_chk{iS,1}));
    if ~isempty(rcn_nme_ind) && ~isempty(rcn_hld{ rcn_nme_ind, 3})
        rcn_nme_sbj = rcn_hld{ rcn_nme_ind, 3};
        rcn_nme_sbj_cut = strfind( rcn_nme_sbj, '_');
        rcn_nme_sbj = rcn_nme_sbj(1:rcn_nme_sbj_cut(end)-1);
        
        if strcmpi(sbj_hld_dta{iS,3},'yes'); dta_chk{iS,3} = 1; else dta_chk{iS,3} = 0;  end
        dti_ind = find( strcmpi( dta_chk{iS,1}, mcd_rsi_dti_prj(:,1) ));
        if numel(dti_ind) > 1;  use_ind = string_find( mcd_rsi_dti_prj(dti_ind,2), rcn_nme_sbj );
            if ~isnan(mcd_rsi_dti_prj{dti_ind(use_ind),38}); dta_chk{iS,4} = 1; else dta_chk{iS,4} = 0; end;
        elseif isempty(dti_ind); dta_chk{iS,4} = 0;
        else; if ~isnan(mcd_rsi_dti_prj{dti_ind,38}); dta_chk{iS,4} = 1; else dta_chk{iS,4} = 0; end; end
        
        
        if strcmpi(sbj_hld_dta{iS,4},'yes'); dta_chk{iS,5} = 1; else dta_chk{iS,5} = 0;  end
        rsi_ind = find( strcmpi( dta_chk{iS,1}, mcd_rsi_flx_prj(:,1) ));
        if numel(rsi_ind) > 1;  use_ind = string_find( mcd_rsi_flx_prj(rsi_ind,2), rcn_nme_sbj );
            if ~isnan(mcd_rsi_flx_prj{rsi_ind(use_ind),38}); dta_chk{iS,6} = 1; else dta_chk{iS,6} = 0; end;
        elseif isempty(rsi_ind); dta_chk{iS,6} = 0;
        else; if ~isnan(mcd_rsi_flx_prj{rsi_ind,38}); dta_chk{iS,6} = 1; else dta_chk{iS,6} = 0; end; end
        
        
        if ~isempty(sbj_hld_dta{iS,5}); dta_chk{iS,7} = 1; else dta_chk{iS,7} = 0;  end
        if ~isempty(sbj_hld_dta{iS,6}); dta_chk{iS,8} = 1; else dta_chk{iS,8} = 0;  end
        if ~isempty(sbj_hld_dta{iS,7}); dta_chk{iS,9} = 1; else dta_chk{iS,9} = 0;  end
        if ~isempty(sbj_hld_dta{iS,8}); dta_chk{iS,10} = 1; else dta_chk{iS,10} = 0;  end
        
        if ~isempty(strcmpi( cnn_sbj_nme, [dta_chk{iS,1} '_' 'norm.csv' ])); dta_chk{iS,11} = 1; else dta_chk{iS,11} = 0;  end
        
    else
        
        dta_chk{iS,3} = 0;
        dta_chk{iS,4} = 0;
        dta_chk{iS,5} = 0;
        dta_chk{iS,6} = 0;
        dta_chk{iS,7} = 0;
        dta_chk{iS,8} = 0;
        dta_chk{iS,9} = 0;
        dta_chk{iS,10} = 0;
        dta_chk{iS,11} = 0;
        
    end
    
end

%% Patient #'s
ltl_ind = find(strcmpi( dta_chk(:,2), 'left'));
rtl_ind = find(strcmpi( dta_chk(:,2), 'right'));
con_ind = find(strcmpi( dta_chk(:,2), 'control'));

% Redcap vs MMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redcap DTI
numel( intersect(find(cell2mat(dta_chk(:,3))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,3))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,3))), con_ind) )

% MMPS DTI
numel( intersect(find(cell2mat(dta_chk(:,4))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,4))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,4))), con_ind) )

% Redcap RSI
numel( intersect(find(cell2mat(dta_chk(:,5))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,5))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,5))), con_ind) )

% MMPS RSI
numel( intersect(find(cell2mat(dta_chk(:,6))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,6))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,6))), con_ind) )

% Redcap LM2
numel( intersect(find(cell2mat(dta_chk(:,7))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,7))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,7))), con_ind) )

% MMPS RSI
numel( intersect(find(cell2mat(dta_chk(:,8))), ltl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,8))), rtl_ind) )
numel( intersect(find(cell2mat(dta_chk(:,8))), con_ind) )

% LM2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMPS LM2 Pre DTI
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), ltl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), rtl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), con_ind) , find(cell2mat(dta_chk(:,7))) ) )

% MMPS LM2 Post DTI
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), ltl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), rtl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), con_ind) , find(cell2mat(dta_chk(:,8))) ) )

% MMPS LM2 Pre RSI
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), ltl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), rtl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), con_ind) , find(cell2mat(dta_chk(:,7))) ) )

% MMPS LM2 Post RSI
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), ltl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), rtl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), con_ind) , find(cell2mat(dta_chk(:,8))) ) )

% MMPS LM2 Pre CNN
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), ltl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), rtl_ind) , find(cell2mat(dta_chk(:,7))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), con_ind) , find(cell2mat(dta_chk(:,7))) ) )

% MMPS LM2 Post CNN
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), ltl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), rtl_ind) , find(cell2mat(dta_chk(:,8))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), con_ind) , find(cell2mat(dta_chk(:,8))) ) )

% BNT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMPS BNT Pre DTI
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), ltl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), rtl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), con_ind) , find(cell2mat(dta_chk(:,9))) ) )

% MMPS BNT Post DTI
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), ltl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), rtl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,4))), con_ind) , find(cell2mat(dta_chk(:,10))) ) )

% MMPS BNT Pre RSI
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), ltl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), rtl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), con_ind) , find(cell2mat(dta_chk(:,9))) ) )

% MMPS BNT Post RSI
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), ltl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), rtl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,6))), con_ind) , find(cell2mat(dta_chk(:,10))) ) )

% MMPS BNT Pre CNN
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), ltl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), rtl_ind) , find(cell2mat(dta_chk(:,9))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), con_ind) , find(cell2mat(dta_chk(:,9))) ) )

% MMPS BNT Post CNN
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), ltl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), rtl_ind) , find(cell2mat(dta_chk(:,10))) ) )
numel( intersect(intersect(find(cell2mat(dta_chk(:,11))), con_ind) , find(cell2mat(dta_chk(:,10))) ) )

%% Check Overlaps
ind_hld = {  ltl_ind   rtl_ind   con_ind };
ind_nme = { 'ltl_ind' 'rtl_ind' 'con_ind' };

dta_chk_nme = { 'sbj_nme' 'epd_sde' 'dti_red_cap' 'dti_spd_sht' 'rsi_red_cap' 'rsi_spd_sht' 'pre_lm2' 'pst_lm2' 'pre_bnt' 'pst_bnt' 'cnn_spd_sht' };

fst_ind = 5;
scd_ind = 6;

clear num_ovr dff_ovr
num_ovr{1,1} = '';
num_ovr{1,2} = dta_chk_nme{fst_ind};
num_ovr{1,3} = dta_chk_nme{scd_ind};
num_ovr{1,4} = 'overlap';
for iG = 1:numel(ind_hld)
    
    num_ovr{iG+1,1} = ind_nme{iG};   
    num_ovr{iG+1,2} =  numel( dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind))), ind_hld{iG}), 1 ) );
    num_ovr{iG+1,3} =  numel( dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind))), ind_hld{iG}), 1 ) );
    num_ovr{iG+1,4} = numel( intersect( dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind))), ind_hld{iG}), 1 ), ...
                                      dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind))), ind_hld{iG}), 1 ) ));
    
    fst_nme_hld = dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind))), ind_hld{iG}), 1 );
    scd_nme_hld = dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind))), ind_hld{iG}), 1 );
    dff_hld = setxor( dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind))), ind_hld{iG}), 1 ), ...
        dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind))), ind_hld{iG}), 1 ));
    fst_sbj = intersect(fst_nme_hld,dff_hld);
    scd_sbj = intersect(scd_nme_hld,dff_hld);
    
    dff_ovr{iG} = [ fst_sbj repmat( dta_chk_nme(fst_ind), numel(fst_sbj), 1) ; ...
                    scd_sbj repmat( dta_chk_nme(scd_ind), numel(scd_sbj), 1) ] ;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fst_ind = [ 3 5 ];
scd_ind = [ 4 6 ];

for iN = 1:numel(fst_ind)
    for iG = 1:numel(ind_hld)
        
        
        
        num_ovr.([dta_chk_nme{fst_ind(iN)} '_AND_' dta_chk_nme{scd_ind(iN)} '_FOR_' ind_nme{iG}]) = numel( intersect( dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind(iN)))), ind_hld{iG}), 1 ), ...
                                                                                                                      dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind(iN)))), ind_hld{iG}), 1 ) ));
                                                                                            
        fst_nme_hld = dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind(iN)))), ind_hld{iG}), 1 );
        scd_nme_hld = dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind(iN)))), ind_hld{iG}), 1 );
        dff_hld = setxor( dta_chk( intersect(find(cell2mat(dta_chk(:,scd_ind(iN)))), ind_hld{iG}), 1 ), ...
                          dta_chk( intersect(find(cell2mat(dta_chk(:,fst_ind(iN)))), ind_hld{iG}), 1 ));
        fst_sbj = intersect(fst_nme_hld,dff_hld);
        scd_sbj = intersect(scd_nme_hld,dff_hld);
        
        dff_ovr.([dta_chk_nme{fst_ind(iN)} '_AND_' dta_chk_nme{scd_ind(iN)} '_FOR_' ind_nme{iG}]) = [ fst_sbj repmat( dta_chk_nme(fst_ind(iN)), numel(fst_sbj), 1) ; ...
                                                                                                      scd_sbj repmat( dta_chk_nme(scd_ind(iN)), numel(scd_sbj), 1) ] ;  
        
    end
end


% DTI (Redcap/MMPS)


              
              
              dta_chk( intersect(find(cell2mat(dta_chk(:,3))), rtl_ind), 1  )
dta_chk( intersect(find(cell2mat(dta_chk(:,3))), con_ind), 1  )

% RSI (Redcap/MMPS)






%% MCD_RSI - DTI_all
mcd_rsi_dti_prj = mmil_readtext( [ dta_dir '/' 'MCD_RSI' '/' 'DTI_all.csv' ] );
    
mcd_rsi_dti_prj_nme = mcd_rsi_dti_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), mcd_rsi_dti_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = mcd_rsi_dti_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = mcd_rsi_dti_prj_nme{iC};
    end
end
mcd_rsi_dti_prj_nme_col_typ = unique(col_typ);

size(mcd_rsi_dti_prj)
numel(find(~isnan(cell2mat(mcd_rsi_dti_prj(2:end,30)))))

%% MCD_RSI - DTI_flex
mcd_rsi_flx_prj = mmil_readtext( [ dta_dir '/' 'MCD_RSI' '/' 'DTI_flex.csv' ] );
    
mcd_rsi_flx_prj_nme = mcd_rsi_flx_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), mcd_rsi_flx_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = mcd_rsi_flx_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = mcd_rsi_flx_prj_nme{iC};
    end
end
mcd_rsi_flx_prj_nme_col_typ = unique(col_typ);

%% MCD_RSI - RSI_res_hind_free
mcd_rsi_rsi_hnd_prj = mmil_readtext( [ dta_dir '/' 'MCD_RSI' '/' 'RSI_res_hind_free.csv' ] );
    
mcd_rsi_rsi_hnd_prj_nme = mcd_rsi_rsi_hnd_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), mcd_rsi_rsi_hnd_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = mcd_rsi_rsi_hnd_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = mcd_rsi_rsi_hnd_prj_nme{iC};
    end
end
mcd_rsi_rsi_hnd_prj_nme_col_typ = unique(col_typ);

cell2csv('/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr_hpp/RCI/measures',[ mcd_rsi_rsi_hnd_prj_nme_col_typ(11:36)  ; ...
                                                                                                mcd_rsi_rsi_hnd_prj_nme_col_typ(37:62)  ; ...
                                                                                                mcd_rsi_rsi_hnd_prj_nme_col_typ(115:140) ; ...
                                                                                                mcd_rsi_rsi_hnd_prj_nme_col_typ(145:170) ]);

%% How many subjects?
sbj_lm2 = sbj_hld(:,[1 5]);

[ ~, sbj_lm2_ind ]      = intersect( mcd_rsi_rsi_hnd_prj(:,1), sbj_lm2(:,1) );
mcd_rsi_rsi_hnd_prj_lm2 = mcd_rsi_rsi_hnd_prj([1 sbj_lm2_ind'],:);
mcd_rsi_rsi_hnd_prj_lm2(:,[1 115:140])

mcd_rsi_rsi_hnd_prj_lm2([1 ; find(~isnan(cell2mat(mcd_rsi_rsi_hnd_prj_lm2(2:end,120))))+1], 1)
mcd_rsi_rsi_hnd_prj_lm2([1 ; find(isnan(cell2mat(mcd_rsi_rsi_hnd_prj_lm2(2:end,120))))+1], 1)

size(mcd_rsi_rsi_hnd_prj_lm2)
numel(find(~isnan(cell2mat(mcd_rsi_rsi_hnd_prj_lm2(2:end,120)))))

size(mcd_rsi_rsi_hnd_prj)
numel(find(~isnan(cell2mat(mcd_rsi_rsi_hnd_prj(2:end,120)))))

setxor( mcd_rsi_rsi_hnd_prj_lm2(:,1), sbj_lm2(:,1) )

%% DTI
% dti_dta_all = mmil_readtext( [ dta_dir '/' 'DTI_all_mmilmcdRSI.csv' ] );    
% 
% dti_col_nme = dti_dta_all(1,:);
% 
% dsh_ind_hld = cellfun( @(x) strfind(x,'-'), dti_col_nme, 'uni', 0);
% col_typ = cell( 1, numel(dsh_ind_hld));
% for iC = 1:numel(dsh_ind_hld)
%     if ~isempty(dsh_ind_hld{iC})
%         col_typ{iC} = dti_col_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
%     else
%         col_typ{iC} = dti_col_nme{iC};
%     end
% end
% dti_dta_all_col_typ = unique(col_typ);

%% Dropbox - RSI_Project_Database.csv
drp_box_rsi_prj = mmil_readtext( [ dta_dir '/' 'Dropbox' '/' 'RSI_Project_Database.csv' ] );
    
drp_box_rsi_prj_nme = drp_box_rsi_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), drp_box_rsi_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = drp_box_rsi_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = drp_box_rsi_prj_nme{iC};
    end
end
drp_box_rsi_prj_nme_col_typ = unique(col_typ);

[ drp_box_rsi_prj_nme_col_typ(14:33) ; ...
  drp_box_rsi_prj_nme_col_typ(34:53) ; ...
  drp_box_rsi_prj_nme_col_typ(94:113) ;
  drp_box_rsi_prj_nme_col_typ(118:137) ]

%%
rsi_dta_all = mmil_readtext('/home/ekaestne/PROJECTS/DATA/ROIHOLD/RSI_res_ind_free_mmilmcdRSI.csv');

rsi_col_nme = rsi_dta_all(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), rsi_col_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = rsi_col_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = rsi_col_nme{iC};
    end
end
col_typ = unique(col_typ);

%%
mcd_rsi_dti_prj = mmil_readtext( [ '/home/ekaestne/PROJECTS/DATA/ROIHOLD' '/' 'DTI_all_mmilmcdRSI_2020_12_17.csv' ] );
    
mcd_rsi_dti_prj_nme = mcd_rsi_dti_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), mcd_rsi_dti_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        if numel(dsh_ind_hld{iC})>1
            col_typ{iC} = mcd_rsi_dti_prj_nme{iC}(1:dsh_ind_hld{iC}(2)-1);
        else
            col_typ{iC} = mcd_rsi_dti_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
        end
    else
        col_typ{iC} = mcd_rsi_dti_prj_nme{iC};
    end
end
mcd_rsi_dti_prj_nme_col_typ = unique(col_typ);

size(mcd_rsi_dti_prj)
numel(find(~isnan(cell2mat(mcd_rsi_dti_prj(2:end,30)))))

%%
mcd_rsi_dti_prj = mmil_readtext( [ '/home/ekaestne/PROJECTS/DATA/ROIHOLD' '/' 'MRI_all_mmilmcdRSI_2020_12_17.csv' ] );
    
mcd_rsi_dti_prj_nme = mcd_rsi_dti_prj(1,:);

dsh_ind_hld = cellfun( @(x) strfind(x,'-'), mcd_rsi_dti_prj_nme, 'uni', 0);
col_typ = cell( 1, numel(dsh_ind_hld));
for iC = 1:numel(dsh_ind_hld)
    if ~isempty(dsh_ind_hld{iC})
        col_typ{iC} = mcd_rsi_dti_prj_nme{iC}(1:dsh_ind_hld{iC}(1)-1);
    else
        col_typ{iC} = mcd_rsi_dti_prj_nme{iC};
    end
end
mcd_rsi_dti_prj_nme_col_typ = unique(col_typ);

size(mcd_rsi_dti_prj)
numel(find(~isnan(cell2mat(mcd_rsi_dti_prj(2:end,30)))))