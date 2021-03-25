function [ dta_out , dta_lbl ] = ejk_create_laterality_index( dta_inp , lbl_inp)

%%
rhs_col = string_find( lbl_inp , {'rh_'} ); und_typ = 1;
if isempty(rhs_col); rhs_col = string_find( lbl_inp , {'Right_'} ); und_typ = 5; end
if isempty(rhs_col); rhs_col = string_find( lbl_inp , {'rhs_'} ); und_typ = 4; end
if isempty(rhs_col); rhs_col = string_find( lbl_inp , {'rh-'} ); und_typ = 2; end
if isempty(rhs_col); rhs_col = string_find( lbl_inp , {'R_'} ); und_typ = 3; end
if isempty(rhs_col); error(':p'); end
lhs_col = string_find( lbl_inp , {'lh_'} );
if isempty(lhs_col); lhs_col = string_find( lbl_inp , {'Left_'} ); end
if isempty(lhs_col); lhs_col = string_find( lbl_inp , {'lhs_'} ); end
if isempty(lhs_col); lhs_col = string_find( lbl_inp , {'lh-'} ); end
if isempty(lhs_col); lhs_col = string_find( lbl_inp , {'L_'} ); end
if isempty(lhs_col); error(':p'); end

%%
rhs_nme = sort( lbl_inp(rhs_col) );
lhs_nme = sort( lbl_inp(lhs_col) );

dta_out = nan( size( dta_inp , 1 ) , size( rhs_nme , 2 ) );
dta_lbl = cell( 1 , size(dta_out,2) );   

for iN = 1:size(dta_out,2)
    
    rhs_num = dta_inp( : , strcmpi( lbl_inp , rhs_nme{iN} ) );
        if iscell(rhs_num); rhs_num=cell2mat(rhs_num); end
    lhs_num = dta_inp( : , strcmpi( lbl_inp , lhs_nme{iN} ) );   
        if iscell(lhs_num); lhs_num=cell2mat(lhs_num); end
    
    dta_out( :,iN ) = ( lhs_num - rhs_num ) ./ ( lhs_num + rhs_num );
    
    if     und_typ==1
        dta_lbl{iN} = regexprep( lbl_inp{ strcmpi( lbl_inp , rhs_nme{iN} ) } , 'rh_' , '');
    elseif und_typ==2
        dta_lbl{iN} = regexprep( lbl_inp{ strcmpi( lbl_inp , rhs_nme{iN} ) } , 'rh-' , '');
    elseif und_typ==3
        dta_lbl{iN} = regexprep( lbl_inp{ strcmpi( lbl_inp , rhs_nme{iN} ) } , 'R_' , '');
    elseif und_typ==4
        dta_lbl{iN} = regexprep( lbl_inp{ strcmpi( lbl_inp , rhs_nme{iN} ) } , 'rhs_' , '');
    elseif und_typ==5
        dta_lbl{iN} = regexprep( lbl_inp{ strcmpi( lbl_inp , rhs_nme{iN} ) } , 'Right_' , '');
    end
    
end


end