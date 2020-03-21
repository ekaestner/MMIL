%% Masking
tts_men = plt_men_dat{2}{dta_ind};
anc_men = men_dff_rhs;

tts_msk = find(tts_men==0); % 13887
anc_msk = find(anc_men==0); % 13916
sum([ tts_msk'== anc_msk]) / numel(anc_msk)

%% p-values
scatter( pvl_rhs , plt_pvl_dat{2}{1} )
xlabel('ANCOVA'); ylabel('TTEST');
line([0 1] , [0 1])

size(plt_pvl_dat{1}{dta_ind})
size(pvl_lhs)

size(plt_pvl_dat{1}{dta_ind})
size(men_dff_lhs)
size(men_dff_rhs)

% Lowest pvlaue
[ ~ , low_pvl ] = min(plt_men_dat{2}{1});

plt_men_dat{2}{1}(low_pvl)
plt_pvl_dat{2}{1}(low_pvl)

men_dff_rhs(low_pvl)
pvl_rhs(low_pvl)

FDR(pvl_rhs,.05)
FDR(plt_pvl_dat{2}{1},.05)

%
[~,lhs_ind] = min(plt_pvl_dat{1}{1});
[~,rhs_ind] = min(plt_pvl_dat{2}{1});

plt_pvl_dat{1}{1}(lhs_ind)
plt_pvl_dat{2}{1}(rhs_ind)

pvl_lhs(lhs_ind)
pvl_rhs(rhs_ind)

%
[~,lhs_ind] = min(plt_men_dat{1}{1});
   [~,ind] = sort(plt_men_dat{1}{1});
   lhs_ind = ind(3);
[~,rhs_ind] = min(plt_men_dat{2}{1});

plt_pvl_dat{1}{1}(lhs_ind)
plt_pvl_dat{2}{1}(rhs_ind)

pvl_lhs(lhs_ind)
pvl_rhs(rhs_ind)

%
figure()
subplot(2,1,1)
hist(grp_dta.lhs.Diagnosis.EPD_Old.dta(:,lhs_ind),100); mean(grp_dta.lhs.Diagnosis.EPD_Old.dta(:,lhs_ind)); xlim([0 5]);
subplot(2,1,2)
hist(grp_dta.lhs.Diagnosis.HC.dta(:,lhs_ind),100); mean(grp_dta.lhs.Diagnosis.HC.dta(:,lhs_ind)); xlim([0 5]);

figure()
subplot(2,1,1)
hist(grp_dta.rhs.Diagnosis.EPD_Old.dta(:,rhs_ind),100); mean(grp_dta.rhs.Diagnosis.EPD_Old.dta(:,rhs_ind)); xlim([0 5]);
subplot(2,1,2)
hist(grp_dta.rhs.Diagnosis.HC.dta(:,rhs_ind),100); mean(grp_dta.rhs.Diagnosis.HC.dta(:,rhs_ind)); xlim([0 5]);

%
bdd_sbj_ind = find(grp_dta.lhs.Diagnosis.EPD_Old.dta(:,lhs_ind)==0);
grp_dta.lhs.Diagnosis.EPD_Old.sbj_nme(bdd_sbj_ind)

find(grp_dta.lhs.Diagnosis.EPD_Old.dta(bdd_sbj_ind(1),:)==0)
find(grp_dta.lhs.Diagnosis.EPD_Old.dta(bdd_sbj_ind(1)+1,:)==0)

numel(find(grp_dta.lhs.Diagnosis.EPD_Old.dta(bdd_sbj_ind(3),:)==0))
numel(find(grp_dta.lhs.Diagnosis.EPD_Old.dta(bdd_sbj_ind(3)+1,:)==0))

%% Mean values
scatter( men_dff_rhs, plt_men_dat{2}{dta_ind} )
xlabel('ANCOVA'); ylabel('TTEST');
line([-.2 0.2] , [-.2 0.2])


figure(); hist( men_dff_rhs(1:10000) - plt_men_dat{2}{dta_ind}(1:10000) , 1000 )
