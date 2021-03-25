clear; clc;

%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'scratch';

%
bld_dta = '/home/mmilmcd/data/MCD_BOLD/subjects';

%
red_fle = 'sbj000_total_2019_03_27.csv';

%% Redcap Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);

sbj_dem.sbj_nme = regexprep(sbj_dem.sbj_nme , '_ucsf' , 'SF');

%% Freesurfer Names
bld_nme = dir( [ bld_dta '/' ] );
bld_nme = { bld_nme(:).name };

bld_nme( string_find(bld_nme,'RM_') ) = [];

%%
bld_out = cell( size( sbj_dem.sbj_nme , 1) , 5);
for iS = 1:size( sbj_dem.sbj_nme , 1) 
    
    pot_fld = bld_nme( string_find( bld_nme , sbj_dem.sbj_nme{iS} ) );
    
    bld_out{iS,1} = iS; 
    bld_out{iS,2} = sbj_dem.sbj_nme{iS};
    
    if numel(pot_fld) > 1
        
        bld_out{iS,4} = pot_fld(string_find( pot_fld , {'post'} ));
        pot_fld(string_find( pot_fld , {'post'} )) = [];
        
        if numel(pot_fld) == 1
            bld_out{iS,3} = pot_fld{1};
            bld_out{iS,5} = 2;
        elseif numel(pot_fld) > 1
            bld_out{iS,3} = {pot_fld{:}};
            bld_out{iS,5} = 99;
        end
    elseif numel(pot_fld) == 1
        bld_out{iS,3} = pot_fld{1};
        bld_out{iS,4} = [];
        bld_out{iS,5} = 1;
    else
        error('')
    end
    
end

bld_out( cellfun( @isempty , bld_out(:,3) ) , 3 ) = {''};
bld_out( cellfun( @isempty , bld_out(:,4) ) , 4 ) = {''};
bld_out( cellfun( @isempty , bld_out(:,5) ) , 5 ) = {nan};

%%
num = 5; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 42; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 51; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 65; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 69; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 91; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 92; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 93; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 94; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 95; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 130; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 131; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 136; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 227; num_chs = 1; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 234; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};
num = 235; num_chs = 2; tot_num = 1:numel(bld_out{num,3}); tot_num(num_chs) = []; bld_out{num,4} = {bld_out{num,3}{tot_num}}; bld_out{num,3} = bld_out{num,3}{num_chs};

% bld_out( cellfun( @iscell , bld_out(:,3) ) , : )
% 
% num = 136;
% bld_out_hld = bld_out{num,3}

%% Finalize
for iS = 1:size( sbj_dem.sbj_nme , 1)
    
    if iscell( bld_out{iS,4} )
        bld_out{iS,4} = ejk_cell_2_str(bld_out{iS,4});
    end
    
end

%% Save
cell2csv( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_bold_names.csv'] , bld_out);





