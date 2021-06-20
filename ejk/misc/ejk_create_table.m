function tbl_out = ejk_create_table( cfg )

% load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )
% 
% fcfg = [];
% fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv'];
% fcfg.dta_col = 2;
% [ cln_dta, cln_dta_sbj, cln_dta_col] = ejk_dta_frm( fcfg );
% 
% dta_inp{1} = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Clinical.csv']);
% dta_inp{2} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/ANOVA/Pre_Omnibus_anova/TLE.Controls.pre.anova/output_table.csv');
% dta_inp{3} = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final/SpecificCor/Clinical/Fisher/TLE_Controls_pre_fishers/output_table.csv');
% 
% cfg = [];
% 
% cfg.tbl(1,:) = { 'mean/std,1,controls_pre_3T_allSurg_all,Educ' ...
%                  'mean/std,1,tle_pre_3T_allSurg_left,Educ' ...
%                  'mean/std,1,tle_pre_3T_allSurg_right,Educ' ...
%                  'copy,2,Educ,report'};
% 
% cfg.tbl(2,:) = { 'count,1,controls_pre_3T_allSurg_all,Sex,M/F' ...
%                  'count,1,tle_pre_3T_allSurg_left,Sex,M/F' ...
%                  'count,1,tle_pre_3T_allSurg_right,Sex,M/F' ...
%                  'copy,3,Sex,report'};
%              
% cfg.dta = dta_inp;
% cfg.grp = grp;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbl_out = cell( size(cfg.tbl) );

for iR = 1:size(cfg.tbl,1)
    for iC = 1:size(cfg.tbl,2)
        
        cut_hld = regexp( cfg.tbl{iR,iC},',','split');
        
        switch cut_hld{1}
            case 'mean/std'
                dta_use = cell2mat(cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})));
                tbl_out{iR,iC} = [ num2str(roundsd(nanmean(dta_use),3)) ' (' num2str(roundsd(nanstd(dta_use),2)) ')' ];
            case 'count'
                dta_use = cfg.dta{str2num(cut_hld{2})}(cfg.grp.(cut_hld{3})+1,strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4}));               
                dta_use( cellfun(@isempty,dta_use) ) = '';
                fac_cat = regexp( cut_hld{5},'/','split');
                for iF = 1:numel(fac_cat)
                    num_hld{iF} = num2str(sum(ismember(dta_use,fac_cat{iF})));
                end
                num_hld = strcat(num_hld,'/');
                num_hld = [num_hld{:}];
                tbl_out{iR,iC} = num_hld(1:end-1);
                clear num_hld
            case 'copy'
                if ~any(strcmpi(cfg.dta{str2num(cut_hld{2})}(:,1),cut_hld{3})); cut_hld{3} = mmil_spec_char(cut_hld{3},{'_'},{'.'}); end
                tbl_out{iR,iC} = cfg.dta{str2num(cut_hld{2})}{strcmpi(cfg.dta{str2num(cut_hld{2})}(:,1),cut_hld{3}),strcmpi(cfg.dta{str2num(cut_hld{2})}(1,:),cut_hld{4})};
            case 'empty'
                tbl_out{iR,iC} = '-';
        end
        
    end
end

end