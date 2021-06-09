clear; clc;

prj_dir = '/home/ekaestne/PROJECTS/OUTPUT';
prj_nme = 'PostOperative/Naming';

tot_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/TOT';

%% Load work so far
cog_dta     = mmil_readtext([ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive_QC.csv']);
    cog_dta_col = cog_dta(1,2:end);
    cog_dta_sbj = cog_dta(2:end,1);
    cog_dta     = cog_dta(2:end,2:end);
    
dti_dta = mmil_readtext( [prj_dir '/' prj_nme '/' 'Data' '/' 'fiber_FA' '_QC.csv'] );
    dti_dta_col = ejk_fix_column_names(dti_dta(1,5:end));
    dti_dta_sbj = dti_dta(2:end,1);
    dti_dta     = dti_dta(2:end,5:end);

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

sbj_inc = mmil_readtext( [ tot_dir '/' 'UCSF_subject_include.csv' ] );

sbj_usd = mmil_readtext( [ tot_dir '/' 'UCSD_TOT_Anny.csv'] );
    sbj_usd( string_find( sbj_usd(:,1), 'ucsf' ), :) = [];
    sbj_usd_col = sbj_usd(1,2:end);
    sbj_usd_sbj = sbj_usd(2:end,1);
    sbj_usd     = sbj_usd(2:end,2:end);

sbj_usf = mmil_readtext( [ tot_dir '/' 'UCSF_TOT_Sanam_v2.csv' ] );
    sbj_usf( cellfun( @isempty, sbj_usf(:,2) ), : ) = [];
    sbj_usf_col = sbj_usf(1,2:end);
    sbj_usf_sbj = sbj_usf(2:end,2);
    sbj_usf     = sbj_usf(2:end,2:end);

%% Collate

sbj_hld_col = { 'SubjID' 'DTI' 'ANT_Pre' 'BNT_Pre' 'ANT_Post' 'BNT_Post'  ...
                'Pre_ANT_TOT'  'Pre>2'  'Pre_pc_correct'  ...
                'Post_ANT_TOT' 'Post>2' 'Post_pc_correct' ...
                'Pre_BNT_cue' ...
                'Post_BNT_cue' };

sbj_hld = cell( numel( grp.tle_controls_pre_3T_allSurg_all ), numel(sbj_hld_col)+1 );


for iS = 1:numel( grp.tle_controls_pre_3T_allSurg_all )

    % Understand which subjects %%%%%%%%%%%%%%%%%
    sbj_hld{iS, 1} = dti_dta_sbj{ grp.tle_controls_pre_3T_allSurg_all(iS), 1 }; % sbj_exs
    
    sbj_hld{iS, 2} = ~isnan( dti_dta{ grp.tle_controls_pre_3T_allSurg_all(iS), 1 } ); % dti_exs
    
    sbj_hld{iS, 3} = ~isnan( cog_dta{ grp.tle_controls_pre_3T_allSurg_all(iS), 2 } ); %ant_exs
    sbj_hld{iS, 4} = ~isnan( cog_dta{ grp.tle_controls_pre_3T_allSurg_all(iS), 1 } ); %bnt_exs
 
    sbj_hld{iS, 5} = ~isnan( cog_dta{ grp.tle_controls_pre_3T_allSurg_all(iS), 5 } ); %ant_exs
    sbj_hld{iS, 6} = ~isnan( cog_dta{ grp.tle_controls_pre_3T_allSurg_all(iS), 4 } ); %bnt_exs
    
    % Check for pre-scores %%%%%%%%%%%%%%%%%
    if isempty(strfind(sbj_hld{iS, 1}, 'ucsf'))
        
        sbj_ind = find( strcmpi( sbj_usd_sbj , sbj_hld{iS, 1}) );
        
        if ~isempty(sbj_ind)
           sbj_hld{iS,7}  = sbj_usd{ sbj_ind, 4 };
           sbj_hld{iS,8}  = sbj_usd{ sbj_ind, 2 };
           sbj_hld{iS,9}  = sbj_usd{ sbj_ind, 3 };
           sbj_hld{iS,10} = sbj_usd{ sbj_ind, 8 };
           sbj_hld{iS,11} = sbj_usd{ sbj_ind, 6 };
           sbj_hld{iS,12} = sbj_usd{ sbj_ind, 7 };
           sbj_hld{iS,13} = nan;
           sbj_hld{iS,14} = nan;
        elseif isempty(sbj_ind)
            sbj_hld(iS, 7:numel(sbj_hld_col)) = {nan};
        end
        
    elseif ~isempty(strfind(sbj_hld{iS, 1}, 'ucsf'))
        
        sbj_ind = find( strcmpi( sbj_usf_sbj , sbj_hld{iS, 1}) );
        
        if ~isempty(sbj_ind)
           sbj_hld{iS,7}  = sbj_usf{ sbj_ind, 2 };
               if strcmpi(sbj_hld{iS,7},'N'); sbj_hld{iS,7}=nan; end
               str_hld = strfind(sbj_hld{iS,7},';');
               if ~isempty(str_hld); sbj_hld{iS,7} = str2num(sbj_hld{iS,7}(1:str_hld(1)-1)); end
           sbj_hld{iS,8}  = nan;
           sbj_hld{iS,9}  = nan;
           sbj_hld{iS,10} = sbj_usf{ sbj_ind, 7 };
               if strcmpi(sbj_hld{iS,10},'no post op data'); sbj_hld{iS,10} = nan; end
               str_hld = strfind(sbj_hld{iS,10},';');
               if ~isempty(str_hld); sbj_hld{iS,10} = str2num(sbj_hld{iS,10}(1:str_hld(1)-1)); end
           sbj_hld{iS,11} = nan;
           sbj_hld{iS,12} = nan;
           sbj_hld{iS,13} = sbj_usf{ sbj_ind, 6 };
                if isempty(sbj_hld{iS,13}); sbj_hld{iS,13}=nan; end
                if strcmpi(sbj_hld{iS,13},'N'); sbj_hld{iS,13}=nan; end
                str_one = strfind(sbj_hld{iS,13},'+');
                str_two = strfind(sbj_hld{iS,13},'ph');
                if ~isempty(str_one) && ~isempty(str_two); sbj_hld{iS,13} = str2num(sbj_hld{iS,13}(str_one(1)+1:str_two(1)-1)); end
                if ~isempty(str_one) && isempty(str_two);  sbj_hld{iS,13} = str2num(sbj_hld{iS,13}(str_one(1)+1:end)); end
           sbj_hld{iS,14} = sbj_usf{ sbj_ind, 8 };
                if strcmpi(sbj_hld{iS,14},'no post op data'); sbj_hld{iS,14} = nan; end
                if isempty(sbj_hld{iS,14}); sbj_hld{iS,14}=nan; end
                if strcmpi(sbj_hld{iS,14},'N'); sbj_hld{iS,14}=nan; end
                str_one = strfind(sbj_hld{iS,14},'+');
                str_two = strfind(sbj_hld{iS,14},'ph');
                if ~isempty(str_one) && ~isempty(str_two); sbj_hld{iS,14} = str2num(sbj_hld{iS,14}(str_one(1)+1:str_two(1)-1)); end
                if ~isempty(str_one) && isempty(str_two);  sbj_hld{iS,14} = str2num(sbj_hld{iS,14}(str_one(1)+1:end)); end
        elseif isempty(sbj_ind)
            sbj_hld(iS, 7:numel(sbj_hld_col)) = {nan};
        end
        
    end
    
    
end

sbj_hld(cellfun(@isempty,sbj_hld)) = {nan};

%% Check what we have
for iS = 1:numel( grp.tle_controls_pre_3T_allSurg_all )
    
    mss_ind = [];
    
    if isnan(sbj_hld{iS,7})  && sbj_hld{iS,2} && sbj_hld{iS,3}; mss_ind = [ mss_ind 7 ]; end
    if isnan(sbj_hld{iS,8})  && sbj_hld{iS,2} && sbj_hld{iS,3}; mss_ind = [ mss_ind 8 ]; end
    if isnan(sbj_hld{iS,9})  && sbj_hld{iS,2} && sbj_hld{iS,3}; mss_ind = [ mss_ind 9 ]; end
    if isnan(sbj_hld{iS,10}) && sbj_hld{iS,2} && sbj_hld{iS,5}; mss_ind = [ mss_ind 10 ]; end
    if isnan(sbj_hld{iS,11}) && sbj_hld{iS,2} && sbj_hld{iS,5}; mss_ind = [ mss_ind 11 ]; end
    if isnan(sbj_hld{iS,12}) && sbj_hld{iS,2} && sbj_hld{iS,5}; mss_ind = [ mss_ind 12 ]; end
    if isnan(sbj_hld{iS,13}) && sbj_hld{iS,2} && sbj_hld{iS,4}; mss_ind = [ mss_ind 13 ]; end
    if isnan(sbj_hld{iS,14}) && sbj_hld{iS,2} && sbj_hld{iS,6}; mss_ind = [ mss_ind 14 ]; end
        
    mss_col = strcat( sbj_hld_col(mss_ind), '; ');
    if isempty(mss_col)
        mss_col = '';
    else
        mss_col = cat(2,mss_col{:});
    end
    
    sbj_hld{iS,15} = mss_col;  
    
end

cell2csv([ tot_dir '/' 'Collatd_Data.csv'], sbj_hld)

%%
mss_sbj_pre_tot = sbj_hld( string_find( sbj_hld(:,15), {'Pre_ANT_TOT'} ), 1);
mss_sbj_pre_cue = sbj_hld( string_find( sbj_hld(:,15), {'Pre_BNT_cue'} ), 1);
mss_sbj_pst_tot = sbj_hld( string_find( sbj_hld(:,15), {'Post_ANT_TOT'} ), 1);
mss_sbj_pst_cue = sbj_hld( string_find( sbj_hld(:,15), {'Post_BNT_cue'} ), 1);

cell2csv([ tot_dir '/' 'Missing_PreTOT.csv'], mss_sbj_pre_tot)
cell2csv([ tot_dir '/' 'Missing_PreBNTCue.csv'], mss_sbj_pre_cue)

cell2csv([ tot_dir '/' 'Missing_PostTOT.csv'], mss_sbj_pst_tot)
cell2csv([ tot_dir '/' 'Missing_PostBNTCue.csv'], mss_sbj_pst_cue)






