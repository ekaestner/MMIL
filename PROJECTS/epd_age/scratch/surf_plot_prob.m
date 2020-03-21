%%
ttt = fs_load_mgh( [ '/home/mmilmcdRSI/MetaData/ADNI/SurfGroupAvgs/DysnomicVNC/Cortthickness_diff_tval_Dysnomic_Normal Control-sm2819-lh.mgh' ] );
figure(); hist(ttt , 1000); ylim([ 0 900])
size(ttt)
sum(isnan(ttt))
sum(ttt==0)

figure(); hist(plt_dat{1}{4} , 1000); ylim([ 0 900])

hist([sort(ttt)' - sort(plt_dat{1}{4})'],1000)
hist([ttt' - plt_dat{1}{4}'],1000)

%%
ttt = fs_load_mgh( [ '/home/mmilmcdRSI/MetaData/ADNI/SurfGroupAvgs/DysexecutiveVNC/Cortthickness_diff_tval_Dysexecutive_Normal Control-sm2819-lh.mgh' ] );
figure(); hist(ttt , 1000); ylim([ 0 900]); xlim([-9 9]);
size(ttt)
sum(isnan(ttt))
sum(ttt==0)

figure(); hist(plt_dat{2}{5} , 1000); ylim([ 0 900]); xlim([-9 9]);
size(plt_dat{2}{5})
sum(isnan(plt_dat{2}{5}))
sum(plt_dat{2}{5}==0)
plt_dat{2}{5}(isnan(plt_dat{2}{5})) = 0;

sum(sort(ttt)==sort(plt_dat{2}{5}))
hist([sort(ttt)' - sort(plt_dat{2}{5})'],1000)
hist([ttt' - plt_dat{1}{5}'],1000)

%%
ttt = fs_load_mgh( [ '/home/mmilmcdRSI/MetaData/ADNI/SurfGroupAvgs/AmnesticVNC/Cortthickness_diff_tval_Amnestic_Normal Control-sm2819-lh.mgh' ] );
figure(); hist(ttt , 1000); ylim([ 0 900])
size(ttt)
sum(isnan(ttt))
sum(ttt==0)

figure(); hist(plt_dat{1}{3},1000)

%%
figure(); hist(plt_dat{1}{3},1000)
sum(isnan(plt_dat{1}{1}))

figure(); hist(plt_dat{2}{1},1000)



%%
figure(); hist(vecdiff,1000)
figure(); hist(vecstderr,1000)

figure(); hist(plt_dat{1}{1},1000)
figure(); hist(plt_dat{1}{2},1000)
figure(); hist(plt_dat{1}{3},1000)
figure(); hist(plt_dat{1}{4},1000)
figure(); hist(plt_dat{1}{5},1000)

%%
load('/home/ekaestne/PROJECTS/OUTPUT/EmTest/data_hold.mat')

% Pseudo-Don Amnestic
info{1}.info{1}.data

% EJK Amnestic
grp_dta.lhs.Diagnosis.Amnestic

[ info{1}.info{1}.data(1:20,1) grp_dta.lhs.Diagnosis.Amnestic.dta(1:20,1) ]


% Pseudo-Don Control
info{2}.info{1}

% EJK Control
grp_dta.lhs.Diagnosis.Amnestic

[ sort(info{2}.info{1}.data(:,1)) sort(grp_dta.lhs.Diagnosis.NormalControl.dta(:,1)) ]

sum(sort(info{2}.info{1}.data(:,1))==sort(grp_dta.lhs.Diagnosis.NormalControl.dta(:,1))) / numel(info{2}.info{1}.data(:,1))

sum(isnan(grp_dta.lhs.Diagnosis.NormalControl.dta(:,1)))
sum(isnan(grp_dta.lhs.Diagnosis.EPD.dta(:,1)))
sum(isnan(grp_dta.lhs.Diagnosis.ClusterDerivedNormal.dta(:,1)))
sum(isnan(grp_dta.lhs.Diagnosis.Amnestic.dta(:,1)))
sum(isnan(grp_dta.lhs.Diagnosis.Dysnomic.dta(:,1)))
sum(isnan(grp_dta.lhs.Diagnosis.Dysexecutive.dta(:,1)))

%%
subplot(3,1,1)
hist(grp_dta.lhs.Diagnosis.Dysexecutive.dta(:,100000),30)
subplot(3,1,2)
hist(grp_dta.lhs.Diagnosis.NormalControl.dta(:,100000),30)
subplot(3,1,3)
hist(grp_dta.lhs.Diagnosis.EPD.dta(:,100000),30)

subplot(3,1,1)
hist(grp_dta.lhs.Diagnosis.Dysexecutive.dta(:,5000),30)
subplot(3,1,2)
hist(grp_dta.lhs.Diagnosis.NormalControl.dta(:,5000),30)
subplot(3,1,3)
hist(grp_dta.lhs.Diagnosis.EPD.dta(:,5000),30)

%
subplot(3,1,1)
hist(grp_dta.lhs.Diagnosis.Dysexecutive.dta(20,:),30)
subplot(3,1,2)
hist(grp_dta.lhs.Diagnosis.NormalControl.dta(20,:),30)
subplot(3,1,3)
hist(grp_dta.lhs.Diagnosis.EPD.dta(20,:),30)

subplot(3,1,1)
hist(grp_dta.lhs.Diagnosis.Dysexecutive.dta(7,:),30)
subplot(3,1,2)
hist(grp_dta.lhs.Diagnosis.NormalControl.dta(7,:),30)
subplot(3,1,3)
hist(grp_dta.lhs.Diagnosis.EPD.dta(7,:),30)

subplot(3,1,1)
hist(grp_dta.lhs.Diagnosis.Dysexecutive.dta(51,:),30)
subplot(3,1,2)
hist(grp_dta.lhs.Diagnosis.NormalControl.dta(51,:),30)
subplot(3,1,3)
hist(grp_dta.lhs.Diagnosis.EPD.dta(51,:),30)

%%
figure()
subplot(5,1,1)
hist(plt_dat{1}{1},1000)
subplot(5,1,2)
hist(plt_dat{1}{2},1000)
subplot(5,1,3)
hist(plt_dat{1}{3},1000)
subplot(5,1,4)
hist(plt_dat{1}{4},1000)
subplot(5,1,5)
hist(plt_dat{1}{5},1000)

%%
figure()
s
hist(plt_dat{1}{1},1000)

%%
string_find(lhs_dta.srf_dta_sbj,{'UCSD'})
string_find(lhs_dta.srf_dta_sbj,{'epd'})
string_find(lhs_dta.srf_dta_sbj,{'_S_'})

figure(); 
subplot(3,2,1)
hist(lhs_dta.srf_dta(675,:),1000)
subplot(3,2,2)
hist(lhs_dta.srf_dta(676,:),1000)
subplot(3,2,3)
hist(lhs_dta.srf_dta(812,:),1000)
subplot(3,2,4)
hist(lhs_dta.srf_dta(728,:),1000)
subplot(3,2,5)
hist(lhs_dta.srf_dta(555,:),1000)
subplot(3,2,6)
hist(lhs_dta.srf_dta(525,:),1000)

%%
figure()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
subplot(5,5,1)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI fc027 thickness');

subplot(5,5,2)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI fc027 sphere-sm2819');

subplot(5,5,4)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd fc027 thickness');

subplot(5,5,5)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd fc027 thickness-sphere-sm2819-lh');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,5,6)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI ucsf013 thickness');

subplot(5,5,7)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI ucsf013 sphere-sm2819');
%
subplot(5,5,9)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd ucsf013 thickness');

subplot(5,5,10)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd_ucsf013_fmri_141119_20141119.101553_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd ucsf013 sphere-sm2819');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,5,11)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd095_fmri_180625_20180625.171911_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI epd095 thickness');

subplot(5,5,12)
ttt = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_epd095_fmri_180625_20180625.171911_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI epd095 sphere-sm2819');
%
subplot(5,5,14)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd095_fmri_180625_20180625.171911_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd epd095 thickness');

subplot(5,5,15)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd095_fmri_180625_20180625.171911_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd epd095 sphere-sm2819');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,5,16)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf/FSURF_EEP006_1_20090127_20090127.123234.218000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('Austin EEP006 thickness');

subplot(5,5,17)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/MRI_aging/fsurf/FSURF_EEP006_1_20090127_20090127.123234.218000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('Austin EEP006 sphere-sm2819');
%
subplot(5,5,19)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_EEP006-1_20090127_20090127.123234.218000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd EEP006 thickness');

subplot(5,5,20)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_EEP006-1_20090127_20090127.123234.218000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd EEP006 sphere-sm2819');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,5,21)
ttt = fs_load_mgh('/space/syn02/1/data/MMILDB/carrierm/fsurf/FSURF_022_S_4173_20110913.132811.781000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI ADNI4173 thickness');

subplot(5,5,22)
ttt = fs_load_mgh('/space/syn02/1/data/MMILDB/carrierm/fsurf/FSURF_022_S_4173_20110913.132811.781000_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI ADNI4173 sphere-sm2819');
%
subplot(5,5,24)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd epd020 thickness');

subplot(5,5,25)
ttt = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_epd020_070330_20070331.150414_1/analysis/thickness-sphere-sm2819-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcd epd020 sphere-sm2819');

tightfig();

%%
tt1 = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-lh.mgz');
tt2 = fs_load_mgh('/space/syn09/1/data/MMILDB/MCD_RSI/fsurf/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-sphere-sm2819-lh.mgz');

tt3 = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-lh.mgz');
tt4 = fs_load_mgh('/space/md18/3/data/MMILDB/EPIPROJ/FSRECONS/FSURF_fc027_meg_071105_20071105.132339_1/analysis/thickness-sphere-sm2819-lh.mgz');

%
size(tt1)
size(tt2)
size(tt3)
size(tt4)

%%
ttt = fs_load_mgh('/space/syn02/1/data/MMILDB/carrierm/fsurf/FSURF_022_S_4173_20110913.132811.781000_1/analysis/thickness-lh.mgz');
hist(ttt,1000); xlim([-0.25 5]); title('mmilmcdRSI ADNI4173 thickness');
print()
close





