% cfg = [];
%
% cfg.srf_dir = '/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/';
% cfg.sbj_nme = {'fc096' 'epd082'};
% cfg.sbj_dir = {'fc096_fmri_170728_20170728.171240_1' 'epd082_fmri2_160720_20160720.161406_1'};

function [ sbj_vol_dta , vol_nme ] = ejk_extract_volumes(cfg)

%%
sbj_dir_lst = dir(sprintf('%s/FSURF*',cfg.ovr_dir));
sbj_dir_lst = regexp({sbj_dir_lst.name},['FSURF_' cfg.sbj_fsr_dir '.+_1$'],'match'); sbj_dir_lst = [sbj_dir_lst{:}];

if ~isempty(sbj_dir_lst)
    
    try
        
        vol_dta = mmil_readtext([cfg.ovr_dir '/' sbj_dir_lst{1} '/' 'stats' '/' 'aseg.stats'],['\t']);
        
        int_crn_vol_ind = find(~cellfun(@isempty,strfind( vol_dta , 'Estimated Total Intracranial Volume' )));
        if isempty(int_crn_vol_ind)
            int_crn_vol_ind = find(~cellfun(@isempty,strfind( vol_dta , 'ICV, Intracranial Volume, ' )));
        end
            int_crn_vol_str = strfind(vol_dta{int_crn_vol_ind},'Intracranial Volume, ');
            int_crn_vol_end = strfind(vol_dta{int_crn_vol_ind},', mm^3');
            int_crn_vol  = str2num(vol_dta{int_crn_vol_ind}(int_crn_vol_str+21:int_crn_vol_end-1));
            
        beg_row = string_find(vol_dta,'# ColHeaders  Index');
        vol_hed = regexp(vol_dta(beg_row),' +','split'); vol_hed = vol_hed{1}(3:end);
        
        vol_dta = vol_dta(beg_row+1:end);
        vol_dta = regexp(vol_dta,' +','split');
        vol_dta = vertcat(vol_dta{:});
        
        vol_col = find(strcmpi(vol_hed,'NVoxels'));
        nme_col = find(strcmpi(vol_hed,'StructName'));
        
        for iFC = 1:size(vol_dta,1)
            sbj_vol_dta(1,iFC) = str2num(vol_dta{iFC,vol_col});
            vol_nme(1,iFC) = vol_dta(iFC,nme_col);
        end
        
        sbj_vol_dta(1,end+1) = int_crn_vol;
        vol_nme = [ vol_nme 'ICV' ];
        
    catch
        
        sbj_vol_dta = [];
        vol_nme     = [];
        
    end
else
    
    sbj_vol_dta = [];
    vol_nme     = [];
    
end

end