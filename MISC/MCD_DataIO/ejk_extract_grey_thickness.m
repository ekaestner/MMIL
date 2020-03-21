% fcfg = [];
%
% fcfg.srf_dir = '/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/';
% fcfg.sbj_nme = {'fc096' 'epd082'};
% fcfg.sbj_dir = {'fc096_fmri_170728_20170728.171240_1' 'epd082_fmri2_160720_20160720.161406_1'};
%
% fcfg.prc_nme = '';

function [ thk_avg , col_nme ] = ejk_extract_grey_thickness(cfg)

if ~isfield(cfg,'min_val'); cfg.min_val = 1e-06; end

%%
sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.ovr_dir));
sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];

if ~isempty(sbj_dir_lst)
    
    try
        
        [ lhs_col_loc , lhs_albl , lhs_actbl ] = fs_read_annotation( [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'lh.aparc' cfg.prc_nme '.annot' ] );
        [ rhs_col_loc , rhs_albl , rhs_actbl ] = fs_read_annotation( [ cfg.prj_dir '/' 'DATA' '/' cfg.sbj_nme '/' 'ROIs' '/' 'rh.aparc' cfg.prc_nme '.annot' ] );
        
        lhs_srf_dta = mmil_rowvec(fs_load_mgh( [ cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'analysis' '/' 'thickness-lh.mgz' ]));
        rhs_srf_dta = mmil_rowvec(fs_load_mgh( [ cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'analysis' '/' 'thickness-rh.mgz' ]));
        
        for iFC = 1:size( lhs_albl , 1 )
            
            lhs_loc = find( lhs_col_loc == find(strcmpi( lhs_albl , lhs_albl{iFC} )));
            rhs_loc = find( rhs_col_loc == find(strcmpi( rhs_albl , rhs_albl{iFC} )));
            
            if isempty(lhs_loc) && isempty(rhs_loc)
                sbj_gry_dta_lhs(1,iFC) = nan;
                sbj_gry_dta_rhs(1,iFC) = nan;
            else
                sbj_gry_dta_lhs(1,iFC) = nanmean(lhs_srf_dta(lhs_loc));
                sbj_gry_dta_rhs(1,iFC) = nanmean(rhs_srf_dta(rhs_loc));          
            end
            
            lhs_col_nme(1,iFC) = lhs_albl(iFC);
            rhs_col_nme(1,iFC) = rhs_albl(iFC); 
            
        end
        
        col_nme = [ strcat('lhs','_',lhs_col_nme) strcat('rhs','_',rhs_col_nme) ];
        thk_avg = [ sbj_gry_dta_lhs        sbj_gry_dta_rhs ];
        
    catch
        
        col_nme = '';
        thk_avg = nan(1,size(col_nme,2));
        
    end
    
else
    
    col_nme = '';
    thk_avg = nan(1,size(col_nme,2));
    
end


end
