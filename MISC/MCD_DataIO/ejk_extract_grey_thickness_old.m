% fcfg = [];
%
% fcfg.srf_dir = '/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/';
% fcfg.sbj_nme = {'fc096' 'epd082'};
% fcfg.sbj_dir = {'fc096_fmri_170728_20170728.171240_1' 'epd082_fmri2_160720_20160720.161406_1'};
%
% fcfg.prc_nme = '';

function [ thk_avg , col_nme ] = ejk_extract_grey_thickness(cfg)

%%
sbj_dir_lst = dir(sprintf('%s/FSURF_*',cfg.ovr_dir));
sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];

%%
if ~isempty(sbj_dir_lst)
    
    try
        
        lhs_gry_dta = mmil_readtext([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'stats' '/' 'lh' '.' 'aparc' cfg.prc_nme '.' 'stats'],['\t']);
        rhs_gry_dta = mmil_readtext([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'stats' '/' 'rh' '.' 'aparc' cfg.prc_nme '.' 'stats'],['\t']);
        
        lhs_beg_row = string_find(lhs_gry_dta,'# ColHeaders StructName');
        lhs_gry_hed = regexp(lhs_gry_dta(lhs_beg_row),' +','split'); lhs_gry_hed = lhs_gry_hed{1}(3:end);
        
        rhs_beg_row = string_find(rhs_gry_dta,'# ColHeaders StructName');
        rhs_gry_hed = regexp(rhs_gry_dta(rhs_beg_row),' +','split'); rhs_gry_hed = rhs_gry_hed{1}(3:end);
        
        lhs_gry_dta = lhs_gry_dta(lhs_beg_row+1:end);
        lhs_gry_dta = regexp(lhs_gry_dta,' +','split');
        lhs_gry_dta = vertcat(lhs_gry_dta{:});
        
        rhs_gry_dta = rhs_gry_dta(rhs_beg_row+1:end);
        rhs_gry_dta = regexp(rhs_gry_dta,' +','split');
        rhs_gry_dta = vertcat(rhs_gry_dta{:});
        
        gry_col = find(strcmpi(lhs_gry_hed,'ThickAvg'));
        nme_col = find(strcmpi(lhs_gry_hed,'StructName'));
        
        for iFC = 1:size(lhs_gry_dta,1)
            sbj_gry_dta_lhs(1,iFC) = str2num(lhs_gry_dta{iFC,gry_col});
            gry_nme_lhs(1,iFC) = lhs_gry_dta(iFC,nme_col);
        end
        for iFC = 1:size(rhs_gry_dta,1)
            sbj_gry_dta_rhs(1,iFC) = str2num(rhs_gry_dta{iFC,gry_col});
            gry_nme_rhs(1,iFC) = rhs_gry_dta(iFC,nme_col);
        end
        
        col_nme = [ strcat('lhs','_',gry_nme_lhs) strcat('rhs','_',gry_nme_rhs) ];
        thk_avg = [ sbj_gry_dta_lhs        sbj_gry_dta_rhs ];
        
    catch
        
        col_nme = '';
        thk_avg = nan(1,size(col_nme,2));
        
    end
    
else
    
    col_nme = '';
    thk_avg = nan(1,size(col_nme,2));
    
end

