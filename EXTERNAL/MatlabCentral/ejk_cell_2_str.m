function cll_out = ejk_cell_2_str(cll_inp)

    cll_inp = strcat( cll_inp , ';' );
    cll_inp = [cll_inp{:}];
    cll_out = cll_inp(1:end-1);

end