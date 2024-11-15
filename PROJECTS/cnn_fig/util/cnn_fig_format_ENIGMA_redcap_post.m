% [ {'sbj_nme'} {'red_cap_nme'} key_fle(1,2:end) ; ses_dta(:,1) cov_dta_out ]

chk_col = [{'red_cap_nme'} key_fle(1,2:end)];
chk_sbj = ses_dta(:,1);

%% Fix missing variables
red_col = strcmpi(chk_col,'red_cap_nme');
img_col = strcmpi(chk_col,'Imaging_Name');
epd_col = strcmpi(chk_col,'Diagnosis');
typ_col = strcmpi(chk_col,'Focal');
loc_col = strcmpi(chk_col,'Localization');
sde_col = strcmpi(chk_col,'Lateralization');
mts_col = strcmpi(chk_col,'MTS');

sdx_col = strcmpi(chk_col,'Legacy');
sdx_fil = ~cellfun(@isempty,cov_dta_out(:,sdx_col));

epd_cnt = 0;
mss_fcl = cell(numel(chk_sbj),4);
mss_loc = cell(numel(chk_sbj),4);
mss_sde = cell(numel(chk_sbj),4);
wrg_sde = cell(numel(chk_sbj),4);
mss_mts = cell(numel(chk_sbj),4);
wrg_mts = cell(numel(chk_sbj),4);

prb_sbj = cell(0); prb_ind = 1;
for iS = 1:numel(chk_sbj)
    
    % Check Diagnosis %%%   
    if isempty(cov_dta_out{iS,epd_col}) && sdx_fil(iS)
        prb_sbj{prb_ind} = [ 'Missing Expected Diagnosis' ' | ' 'subject: ' chk_sbj{iS}];
        prb_ind = prb_ind + 1;
    elseif ~isempty(cov_dta_out{iS,epd_col}) && sdx_fil(iS)
        if strcmpi(cov_dta_out{iS,epd_col},'HC') && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')
            prb_sbj{prb_ind} = ['Diagnosis: Expected HC' ' | ' 'subject: ' chk_sbj{iS}];
            prb_ind = prb_ind + 1;
        elseif ~strcmpi(cov_dta_out{iS,epd_col},'HC') && strcmpi(cov_dta_out{iS,sdx_col},'HC')
            prb_sbj{prb_ind} = ['Diagnosis: Expected EPD' ' | ' 'subject: ' chk_sbj{iS}];
            prb_ind = prb_ind + 1;
        end
    end

    % Check Focal %%%
    if isempty(cov_dta_out{iS,typ_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')
        mss_fcl{iS,1} = cov_dta_out{iS,red_col};
        mss_fcl{iS,2} = cov_dta_out{iS,img_col};
        mss_fcl{iS,3} = cov_dta_out{iS,typ_col};
        mss_fcl{iS,4} = cov_dta_out{iS,sdx_col};
    elseif ~isempty(cov_dta_out{iS,typ_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')
        if strcmpi(cov_dta_out{iS,typ_col},'GGE') && ~strcmpi(cov_dta_out{iS,sdx_col},'generalized')
            error('Focal: Expected GGE')
        end
    end
    
    % Check Localization %%% 
    if isempty(cov_dta_out{iS,loc_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC') && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE')
        mss_loc{iS,1} = cov_dta_out{iS,red_col};
        mss_loc{iS,2} = cov_dta_out{iS,img_col};
        mss_loc{iS,3} = cov_dta_out{iS,loc_col};
        mss_loc{iS,4} = cov_dta_out{iS,sdx_col};
    elseif ~isempty(cov_dta_out{iS,loc_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')  && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE')
        if strcmpi(cov_dta_out{iS,loc_col},'temporal') && ~( strcmpi(cov_dta_out{iS,sdx_col},'bilateral_TLE') || strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_R') || strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_L') || strcmpi(cov_dta_out{iS,sdx_col},'TLE_L') || strcmpi(cov_dta_out{iS,sdx_col},'TLE-R'))
            error('Localization: Expected temporal')
        elseif strcmpi(cov_dta_out{iS,loc_col},'ex-tle-frontal') && ~( strcmpi(cov_dta_out{iS,sdx_col},'extra_TLE_focal') ) 
            error('Localization: Expected ex-temporal')
        end
    end

    % Check Side %%%
    if isempty(cov_dta_out{iS,sde_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC') && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE') && ~strcmpi(cov_dta_out{iS,loc_col},'extra_TLE_focal')
        mss_sde{iS,1} = cov_dta_out{iS,red_col};
        mss_sde{iS,2} = cov_dta_out{iS,img_col};
        mss_sde{iS,3} = cov_dta_out{iS,sde_col};
        mss_sde{iS,4} = cov_dta_out{iS,sdx_col};
    elseif ~isempty(cov_dta_out{iS,sde_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')  && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE') && ~strcmpi(cov_dta_out{iS,loc_col},'extra_TLE_focal')
        if strcmpi(cov_dta_out{iS,sde_col},'left') && ~( strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_L') || strcmpi(cov_dta_out{iS,sdx_col},'TLE_L') )
            wrg_sde{iS,1} = cov_dta_out{iS,red_col};
            wrg_sde{iS,2} = cov_dta_out{iS,img_col};
            wrg_sde{iS,3} = cov_dta_out{iS,sde_col};
            wrg_sde{iS,4} = cov_dta_out{iS,sdx_col};
        elseif strcmpi(cov_dta_out{iS,sde_col},'right') && ~( strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_R') || strcmpi(cov_dta_out{iS,sdx_col},'TLE-R'))
            wrg_sde{iS,1} = cov_dta_out{iS,red_col};
            wrg_sde{iS,2} = cov_dta_out{iS,img_col};
            wrg_sde{iS,3} = cov_dta_out{iS,sde_col};
            wrg_sde{iS,4} = cov_dta_out{iS,sdx_col};
        elseif strcmpi(cov_dta_out{iS,sde_col},'bilateral') && ~( strcmpi(cov_dta_out{iS,sdx_col},'bilateral_TLE') )
            wrg_sde{iS,1} = cov_dta_out{iS,red_col};
            wrg_sde{iS,2} = cov_dta_out{iS,img_col};
            wrg_sde{iS,3} = cov_dta_out{iS,sde_col};
            wrg_sde{iS,4} = cov_dta_out{iS,sdx_col};
        end
    end

    % Check MTS %%%
    if isempty(cov_dta_out{iS,mts_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC') && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE') && ~strcmpi(cov_dta_out{iS,loc_col},'extra_TLE_focal') && ~strcmpi(cov_dta_out{iS,sdx_col},'extra_TLE_focal') && ~strcmpi(cov_dta_out{iS,sde_col},'bilateral') && ~strcmpi(cov_dta_out{iS,sdx_col},'bilateral_TLE')
        mss_mts{iS,1} = cov_dta_out{iS,red_col};
        mss_mts{iS,2} = cov_dta_out{iS,img_col};
        mss_mts{iS,3} = cov_dta_out{iS,mts_col};
        mss_mts{iS,4} = cov_dta_out{iS,sdx_col};
    elseif ~isempty(cov_dta_out{iS,mts_col}) && sdx_fil(iS) && ~strcmpi(cov_dta_out{iS,sdx_col},'HC')  && ~strcmpi(cov_dta_out{iS,sdx_col},'GGE') && ~strcmpi(cov_dta_out{iS,loc_col},'extra_TLE_focal') && ~strcmpi(cov_dta_out{iS,sdx_col},'extra_TLE_focal') && ~strcmpi(cov_dta_out{iS,sde_col},'bilateral') && ~strcmpi(cov_dta_out{iS,sdx_col},'bilateral_TLE')
        if ( strcmpi(cov_dta_out{iS,mts_col},'HS_pathology') && strcmpi(cov_dta_out{iS,mts_col},'HS_mri') ) && ~( strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_L') || strcmpi(cov_dta_out{iS,sdx_col},'TLE_MTS_R') )
            wrg_mts{iS,1} = cov_dta_out{iS,red_col};
            wrg_mts{iS,2} = cov_dta_out{iS,img_col};
            wrg_mts{iS,3} = cov_dta_out{iS,mts_col};
            wrg_mts{iS,4} = cov_dta_out{iS,sdx_col};
        elseif strcmpi(cov_dta_out{iS,mts_col},'noHS') && ~( strcmpi(cov_dta_out{iS,sdx_col},'TLE_L') || strcmpi(cov_dta_out{iS,sdx_col},'TLE-R') )
            wrg_mts{iS,1} = cov_dta_out{iS,red_col};
            wrg_mts{iS,2} = cov_dta_out{iS,img_col};
            wrg_mts{iS,3} = cov_dta_out{iS,mts_col};
            wrg_mts{iS,4} = cov_dta_out{iS,sdx_col};
        end
    end

end

mss_fcl( cellfun(@isempty,mss_fcl(:,1)),:) = [];
mss_loc( cellfun(@isempty,mss_loc(:,1)),:) = [];
mss_sde( cellfun(@isempty,mss_sde(:,1)),:) = [];
wrg_sde( cellfun(@isempty,wrg_sde(:,1)),:) = [];
mss_mts( cellfun(@isempty,mss_mts(:,1)),:) = [];
wrg_mts( cellfun(@isempty,wrg_mts(:,1)),:) = [];



