[~, pvl_hld, ~, stt_hld] = ttest2( scr_num(grp.side.L), scr_num(grp.side.R));
stt_tbl(1,:) = { 'Number of Times: L vs R' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };

[ rvl_hld, pvl_hld ]  = corrcoef(scr_num,scr_hld,'Rows','complete');
    deg_fre      = numel(intersect( find(~isnan(scr_num)), find(~isnan(scr_hld)) ));
    stt_tbl(5,:) = { 'Correct -BY- Number: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

[ pvl_hld, tst_hld, stt_hld ] = anova1(scr_hld, dem_dta(:,strcmpi(dem_dta_col,'Site')),'off');
    stt_tbl(25,:) = { 'Correct: Site, all' ['F(' num2str(tst_hld{2,3}) ',' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
    
%%
stt_tbl(1,:) = { 'Number of Times: L vs R' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };

stt_tbl(3,:) = { 'Correct: L vs R' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };

stt_tbl(5,:) = { 'Correct -BY- Number: LTLE' ['r(' num2str(deg_fre) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(6,:) = { 'Correct -BY- Number: LTLE' ['r(' num2str(deg_fre) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(7,:) = { 'Correct -BY- Number: RTLE' ['r(' num2str(deg_fre) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(9,:)  = { 'Correct -BY- QC: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(10,:) = { 'Correct -BY- QC: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(11,:) = { 'Correct -BY- QC: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(13,:) = { 'Correct -BY- L Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(14,:) = { 'Correct -BY- L Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(15,:) = { 'Correct -BY- L Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(17,:) = { 'Correct -BY- R Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(18,:) = { 'Correct -BY- R Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(19,:) = { 'Correct -BY- R Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(21,:) = { 'Correct -BY- LI Hip Vol: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(22,:) = { 'Correct -BY- LI Hip Vol: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(23,:) = { 'Correct -BY- LI Hip Vol: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(25,:) = { 'Correct: Site, all' ['F(' num2str(tst_hld{2,3}) ',' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(26,:) = { 'Correct: Site, LTLE' ['F(' num2str(tst_hld{2,3}) ',' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(27,:) = { 'Correct: Site, RTLE' ['F(' num2str(tst_hld{2,3}) ',' num2str(tst_hld{3,3}) ')' '=' num2str(roundsd(tst_hld{2,5},3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };

stt_tbl(29,:) = { 'Correct -BY- Age: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(30,:) = { 'Correct -BY- Age: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(31,:) = { 'Correct -BY- Age: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
 
stt_tbl(33,:) = { 'Score: TLE, M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(34,:) = { 'Score: LTLE, M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(35,:) = { 'Score: RTLE, M vs F' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };

stt_tbl(37,:) = { 'Score: TLE, Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(38,:) = { 'Score: LTLE, Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(39,:) = { 'Score: RTLE, Handedness' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
 
stt_tbl(41,:) = { 'Correct -BY- Education: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };   
stt_tbl(42,:) = { 'Correct -BY- Education: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };    
stt_tbl(43,:) = { 'Correct -BY- Education: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };    

stt_tbl(45,:) = { 'Correct -BY- AgeOfOnset: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(46,:) = { 'Correct -BY- AgeOfOnset: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(47,:) = { 'Correct -BY- AgeOfOnset: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };

stt_tbl(49,:) = { 'Score: TLE, MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(50,:) = { 'Score: LTLE, MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };   
stt_tbl(51,:) = { 'Score: RTLE, MTS' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };  

stt_tbl(53,:) = { 'Correct -BY- ASMs: TLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(54,:) = { 'Correct -BY- ASMs: LTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };
stt_tbl(55,:) = { 'Correct -BY- ASMs: RTLE' ['r(' num2str(deg_fre-2) ')' 'r=' num2str(roundsd(rvl_hld(1,2),2)) ',' 'p=' num2str(roundsd(pvl_hld(1,2),2))] };   
    
stt_tbl(57,:) = { 'Score: TLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(58,:) = { 'Score: LTLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(59,:) = { 'Score: RTLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };  
    
stt_tbl(61,:) = { 'Score: TLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(62,:) = { 'Score: LTLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };
stt_tbl(63,:) = { 'Score: RTLE, Surgery' ['t(' num2str(stt_hld.df) ')' 't=' num2str(roundsd(stt_hld.tstat,3)) ',' 'p=' num2str(roundsd(pvl_hld,2))] };     