%% Clinical
ejk_cln = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/Clinical.csv');
    ejk_cln_col = ejk_cln(1,:);
    ejk_cln     = ejk_cln(2:end,:);

aln_cln = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/AlenaFindings/leftTLE_3T.csv');
    aln_cln_col = aln_cln(1,:);
    aln_cln     = aln_cln(2:end,:);

[ ~, ejk_cln_cph] = intersect(ejk_cln(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));



ejk_col_num = 15;  ejk_cln_col(ejk_col_num)
aln_col_num = 47;

[ ejk_cln_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_cln(ejk_cln_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% DTI-FA
ejk_fib = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/fiber_FA.csv');
    ejk_fib_col = ejk_fib(1,:);
    ejk_fib     = ejk_fib(2:end,:);

[ ~, ejk_fib_cph] = intersect(ejk_fib(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 19;  ejk_fib_col(ejk_col_num)
aln_col_num = 73;

[ ejk_fib_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_fib(ejk_fib_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% DTI-MD LI
ejk_fib = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/fiber_MD_LI.csv');
    ejk_fib_col = ejk_fib(1,:);
    ejk_fib     = ejk_fib(2:end,:);

[ ~, ejk_fib_cph] = intersect(ejk_fib(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 19;  ejk_fib_col(ejk_col_num)
aln_col_num = 83;

[ ejk_fib_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_fib(ejk_fib_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% DTI-FA wmparc
ejk_wmp = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/wmparc_FA_wm_aparc_annot_LI.csv');
    ejk_wmp_col = ejk_wmp(1,:);
    ejk_wmp     = ejk_wmp(2:end,:);

[ ~, ejk_wmp_cph] = intersect(ejk_wmp(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 9;  ejk_wmp_col(ejk_col_num)
aln_col_num = 104;

[ ejk_wmp_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_wmp(ejk_wmp_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% fMRI
ejk_fmr = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/alicia.csv');
    ejk_fmr_col = ejk_fmr(1,:);
    ejk_fmr     = ejk_fmr(2:end,:);

[ ~, ejk_fmr_cph] = intersect(ejk_fmr(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 9;  ejk_fmr_col(ejk_col_num)
aln_col_num = 58;

[ ejk_fmr_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_fmr(ejk_fmr_cph,1) ejk_fmr(ejk_fmr_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% fMRI-LI
ejk_fmr = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/alicia_LI.csv');
    ejk_fmr_col = ejk_fmr(1,:);
    ejk_fmr     = ejk_fmr(2:end,:);

[ ~, ejk_fmr_cph] = intersect(ejk_fmr(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 5;  ejk_fmr_col(ejk_col_num)
aln_col_num = 59;

[ ejk_fmr_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_fmr(ejk_fmr_cph,1) ejk_fmr(ejk_fmr_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]

%% Cognitive
ejk_cog = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Language/Data/Cognitive.csv');
    ejk_cog_col = ejk_cog(1,:);
    ejk_cog     = ejk_cog(2:end,:);

[ ~, ejk_cog_cph] = intersect(ejk_cog(:,1) ,aln_cln(:,2));
[ ~, aln_cln_cph] = sort(aln_cln(:,2));    
    
ejk_col_num = 6;  ejk_cog_col(ejk_col_num)
aln_col_num = 17;

[ ejk_cog_col(ejk_col_num)         aln_cln_col(aln_col_num) ]
[ ejk_cog(ejk_cog_cph,1) ejk_cog(ejk_cog_cph,ejk_col_num) aln_cln(aln_cln_cph,aln_col_num)   ]






