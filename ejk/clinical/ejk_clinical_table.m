% cfg = [];
% 
% cfg.sbj_nme = phe_grp_fle;
% cfg.grp_col = 4; % 1; % 2; % 3; % 4;
% cfg.grp_nme = {'HC' 'EPD'};  % { 'Language Impairment' 'No Language Impairment'  }; % { 'left' 'right' } % {} % { 'HC' 'EPD' }
% 
% cfg.mes_sbj_nme = { sbj_dem.sbj_nme sbj_dem.sbj_nme sbj_dem.sbj_nme     sbj_sze.sbj_nme  sbj_sze.sbj_nme     sbj_sze.sbj_nme     sbj_sze.sbj_nme      sbj_sze.sbj_nme     sbj_cog.sbj_nme      sbj_cog.sbj_nme         sbj_cog.sbj_nme   }; % { sbj_dem.sbj_nme sbj_dem.sbj_nme };
% cfg.mes_nme     = { 'Education'     'Sex: M/F'      'Handedness: L/R/A' 'MTS: Yes/No'    'Side: L/R/B'       'Age of Onset'      'Number of AEDs'     'Seizure Frequency' 'Boston Naming Test' 'Auditory Naming Test'  'D-KEFS Category'           }; % { 'Education'     'Sex: M/F' };
% cfg.mes_lbl     = { sbj_dem.sbj_edu sbj_dem.sbj_sex sbj_dem.sbj_hnd     sbj_sze.sbj_mts  sbj_sze.sbj_sde_ons sbj_sze.sbj_age_ons sbj_sze.sbj_aed_num  sbj_sze.sbj_sze_frq sbj_cog.bnt_nor_scr  sbj_cog.ant_mem_raw_scr sbj_cog.cat_flu_nor_scr  }; % { sbj_dem.sbj_edu sbj_dem.sbj_sex };
% cfg.mes_typ     = { 'con'           'ord'           'ord'               'ord'            'ord'               'con'               'con'                'con'               'con'                'con'                   'con' }; % { 'con'           'ord' };
% cfg.ord_ord     = { {}              { 'M' 'F' }     {'L' 'R' }          {'L' 'R' 'N/A' } {'L' 'R'}           {}                  {}                   {}                  {}                   {}                      {} }; % { {}              { 'M' 'F' } };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ejk_clinical_table(cfg)

clear out_tbl

for iTB = 1:numel(cfg.mes_nme)
    
    dta_hld = [];
    grp_hld = cell(0);
    
    for iGR = 1:numel(cfg.grp_nme)
        
        sbj_nme_int = cfg.sbj_nme(strcmpi(cfg.sbj_nme(:,cfg.grp_col+1),cfg.grp_nme{iGR}),1);
        sbj_nme_ind = ismember(cfg.mes_sbj_nme{iTB},sbj_nme_int);
        
        out_tbl{1,iGR+1} = cfg.grp_nme{iGR};
        out_tbl{2,1}     = 'N';
        out_tbl{2,iGR+1} = sum(sbj_nme_ind);
        
        out_tbl{iTB+2,1} = cfg.mes_nme{iTB};
        
        if strcmpi(cfg.mes_typ{iTB},'con')
            
            men_hld = nanmean(cfg.mes_lbl{iTB}(sbj_nme_ind,1));
            std_hld = nanstd(cfg.mes_lbl{iTB}(sbj_nme_ind,1));
            out_tbl{iTB+2,iGR+1} = [num2str(men_hld) ' (' num2str(std_hld) ')'];
            
            dta_hld = [ dta_hld cfg.mes_lbl{iTB}(sbj_nme_ind,1)' ];
            grp_hld = [ grp_hld repmat({[ 'grp' num2str(iGR) ]} , 1 , numel(cfg.mes_lbl{iTB}(sbj_nme_ind,1)') ) ];
            
        elseif strcmpi(cfg.mes_typ{iTB},'ord')
            
            tab_hld = tabulate(cfg.mes_lbl{iTB}(sbj_nme_ind,1));
            for iTU = 1:numel(cfg.ord_ord{iTB}) 
                if sum(strcmpi(tab_hld(:,1),cfg.ord_ord{iTB}{iTU}))~=0
                    out_tbl{iTB+2,iGR+1}(iTU) = tab_hld{strcmpi(tab_hld(:,1),cfg.ord_ord{iTB}{iTU}),2};
                    dta_hld(iGR,iTU) = tab_hld{strcmpi(tab_hld(:,1),cfg.ord_ord{iTB}{iTU}),2};
                else out_tbl{iTB+2,iGR+1}(iTU) = 0;
                end
            end
            out_tbl{iTB+2,iGR+1} = num2str(out_tbl{iTB+2,iGR+1});
            out_tbl{iTB+2,iGR+1} = regexprep(out_tbl{iTB+2,iGR+1},' *','/');
            
        end
        
    end
    
    if strcmpi(cfg.mes_typ{iTB},'con')
        
        [ pvl_hld , ~ , stt_hld ] = anova1( dta_hld , grp_hld , 'off' );
        pvl_hld = num2str(roundsd(pvl_hld,2));
        out_tbl{iTB+2,iGR+2} = [ 'p = ' pvl_hld(2:end)];
        
    elseif strcmpi(cfg.mes_typ{iTB},'ord')
        
        pvl_hld = myfisher( dta_hld );
        pvl_hld = num2str(roundsd(pvl_hld,2));
        out_tbl{iTB+2,iGR+2} = [ 'p = ' pvl_hld(2:end)];
        
    end
    
end

if ~isdir([ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'demographics' '/' ]); mkdir([ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'demographics' '/' ]); end
cell2csv( [ cfg.prj_dir '/' 'OUTPUT' '/' cfg.prj_nme '/' 'demographics' '/' cfg.tbl_nme ] , out_tbl )

end