
clear; clc;

%
prj_dir = '/home/ekaestne/PROJECTS/';
prj_nme = 'scratch';

%
frs_dta = '/home/mmilmcd/data/FSRECONS/';

%
red_fle = 'sbj000_total_2019_03_27.csv';

%% Redcap Data
fcfg = [];
fcfg.prj_dir = prj_dir;
fcfg.red_fle = red_fle;
[sbj_dem , sbj_sze , sbj_scn , sbj_cog] = mmil_load_redcap(fcfg);

%% Freesurfer Names
frs_nme = dir( [ frs_dta '/' ] );
frs_nme = { frs_nme(:).name };
frs_nme = frs_nme(string_find( frs_nme , { 'FSURF_' } ));

%%
frs_out = cell( size( sbj_dem.sbj_nme , 1) , 5);
for iS = 1:size( sbj_dem.sbj_nme , 1) 
    
    pot_fld = frs_nme( string_find( frs_nme , sbj_dem.sbj_nme{iS} ) );
    
    frs_out{iS,1} = iS; 
    frs_out{iS,2} = sbj_dem.sbj_nme{iS};
    
    if numel(pot_fld) > 1
        
        frs_out{iS,4} = pot_fld(string_find( pot_fld , {'post'} ));
        pot_fld(string_find( pot_fld , {'post'} )) = [];
        
        if numel(pot_fld) == 1
            frs_out{iS,3} = pot_fld{1};
            frs_out{iS,5} = 2;
        elseif numel(pot_fld) > 1
            frs_out{iS,3} = {pot_fld{:}};
            frs_out{iS,5} = 99;
        end
    elseif numel(pot_fld) == 1
        frs_out{iS,3} = pot_fld{1};
        frs_out{iS,4} = [];
        frs_out{iS,5} = 1;
    else
        error('')
    end    
end

frs_out( cellfun( @isempty , frs_out(:,3) ) , 3 ) = {''};
frs_out( cellfun( @isempty , frs_out(:,4) ) , 4 ) = {''};
frs_out( cellfun( @isempty , frs_out(:,5) ) , 5 ) = {nan};

%%
num = 5; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 6; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 18; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 19; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 26; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 34; num_chs = 4; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 40; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 41; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 42; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 45; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 49; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 50; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 51; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 54; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 57; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 58; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 63; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 69; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 76; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 81; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 87; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 89; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 95; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 124; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 127; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 130; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 145; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 147; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 148; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 150; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 160; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 161; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 163; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 166; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 168; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 169; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 170; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 171; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 174; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 175; num_chs = 5; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 177; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 179; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 180; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 181; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 182; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 183; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 189; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 190; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 192; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 195; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 197; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 198; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 199; num_chs = 3; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 222; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 225; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
num = 234; num_chs = 1; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};

% num = 40; num_chs = 2; tot_num = 1:numel(frs_out{num,3}); tot_num(num_chs) = []; frs_out{num,4} = { frs_out{num,4}{:} frs_out{num,3}{tot_num} }; frs_out{num,3} = frs_out{num,3}{num_chs};
% 
% frs_out
% 
% num = 199;
% frs_out_hld = frs_out{num,3}
% fsr_fld{ string_find( fsr_fld(:,1) , frs_out{num,2} ) , 2 }{ fsr_fld{ string_find( fsr_fld(:,1) , frs_out{num,2}) , 3 } }

%% Finalize
frs_out = [ frs_out(:,1:2) cell(size(frs_out,1),1) frs_out(:,3:5) ];

for iS = 1:size( sbj_dem.sbj_nme , 1)
    
    %
    if ~strcmpi(frs_out{iS,4},'')
        und_ind = strfind( frs_out{iS,4} , '_' );
        fst_und_ind = und_ind(1);
        per_ind = strfind( frs_out{iS,4} , '.' );
        per_ind = per_ind(1);
        frs_out{iS,3} = frs_out{iS,4}( fst_und_ind+1 : per_ind-1 );
        und_ind = strfind( frs_out{iS,3} , '_' );
        und_ind = und_ind(end);
        frs_out{iS,3} = frs_out{iS,3}(1:und_ind-1);
    else
        frs_out{iS,3} = '';
    end
    
    %
    if iscell( frs_out{iS,5} )
        frs_out{iS,5} = ejk_cell_2_str(frs_out{iS,5});
    end
    
end

%% Save
cell2csv( [ prj_dir '/' 'SUBJECTS' '/' 'mmilmcdRSI_freesurfer_names.csv'] , frs_out);





