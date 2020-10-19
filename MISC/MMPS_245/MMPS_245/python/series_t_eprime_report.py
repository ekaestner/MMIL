print('Plot histograms of time differences between Eprime and Imaging-Series data')
print('by redading data from a file like ~/MetaData/DAL_ABCD_QC/DAL_ABCD_QC_eprime.csv')
print('')
print('                                          Octavio Ruiz, 2017aug-29, sep28-oct05')
print('Some time_threshold interesting values: 730, 17, 12.5 minutes.', '\n')

#-----------------------------------------------------------------------
# Import pandas, avoiding warning
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
pd.set_option('display.width', 512)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Get data
fname = input("File name (full): ")
dat = pd.read_csv(fname, header=0 )
# Keep variables of interest
dat = dat[['pGUID', 'VisitID', 'StudyDate', 'SeriesTime', 'SiteName', 'Completed', 'AdditionalInfo', 'taskname',
            'eprime_file_name',
            'eprime_datetime_found', 'eprime_file_res', 'eprime_date', 'eprime_time', 'eprime_t_series_mode', 'eprime_t_rel_to_series', 'eprime_t_series_match']]
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
if input('Consider (a)ll runs or only (f)irst run of each task? ') == 'a':
    optn_sel = 'All task runs'
    pass

else:
    optn_sel = 'First runs of each task'
    # Find and keep first rows of a task, per subject and session

    # Sort data
    dat.sort_values( ['pGUID', 'taskname', 'StudyDate', 'SeriesTime'],  inplace=True)

    # Remove second (last) task instance, if more than one was run
    dat_tsk1 = dat.drop_duplicates( ['pGUID', 'taskname'] )
    dat = dat_tsk1
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Condition data

import matplotlib.pyplot as plt

tdif_str = dat['eprime_t_rel_to_series']
tdif = tdif_str.convert_objects(convert_numeric=True)
dat['tdif'] = tdif

thrs_s = input('Threshold (minutes, []=>end): ')
thrs = float(thrs_s)
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Plot global histograms of large and small time values
fig1 = plt.figure(1)
tdif_small = tdif[ tdif.abs() <= thrs]
tdif_large = tdif[ tdif.abs()  > thrs]
plt.subplot(211);  tdif_small.hist( bins=481 )
plt.subplot(212);  tdif_large.hist( bins=481 )
plt.xlabel('eprime_t_rel_to_series (minutes)')
fig1.suptitle( optn_sel + '.  Time differences less (top) or larger (below) than' + ' %.1f min.'%(thrs) )

# Calculate fraction of valid time differences within and beyond trheshold
tshn = len(tdif_small)
tlgn = len(tdif_large)
tsl_rat = tshn / (tshn + tlgn)
print('Num.of cases per t_diff: large:', tlgn, ', short:', tshn,
      '.  proportion_within_t_max = %5.2f' % (tsl_rat*100), '%' )
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Plot histograms of small time values, broken by diverse variables

fig2 = plt.figure(2);   ax2  = fig2.add_subplot(1,1,1)
fig3 = plt.figure(3);   ax3  = fig3.add_subplot(1,1,1)
fig4 = plt.figure(4);   ax4  = fig4.add_subplot(1,1,1)

while len(thrs_s) > 0:
    thrs = float(thrs_s)

    dat_tsml  =  dat[ abs(dat['tdif']) <= thrs ]
    
    dat_tsml.hist(column='tdif', by='SiteName', ax=ax2, sharex=True, bins=481 )
    fig2.suptitle(optn_sel + '.  By Site.')
    
    dat_tsml.hist(column='tdif', by='taskname', ax=ax3, sharex=True, bins=481 )
    fig3.suptitle(optn_sel + '.  By Task.')

    dat_tsml.hist(column='tdif', by='eprime_file_res', ax=ax4, sharex=True, bins=481 )
    fig4.suptitle(optn_sel + '.  By eprime_file_res.')

    plt.show( block=False)

    thrs_s = input('\nThreshold (minutes, []=>end): ')
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
    # dat_lrge  =  dat[ abs(dat['tdif']) >  thrs ]
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# I cannot control Figure size:
# These doesn't work, when I set it before plotting
#   fig1.set_figwidth(16.0)
#   fig1.set_figheight(8.5)
# If I set it after plotting, the result is worse: Window doesn't change size, plot gets cropped
# 
# This doesn't work either (produces no size changes):
#   dat_tsml.hist(column='tdif', by='eprime_file_res',  ax=ax2, figsize=(10,20), sharex=True, bins=481 )

# I cannot make a figures to show a title.  Tried:
#   ...
#   ax1.set_title('per Site')
#   ...
#   plt.title('per Task')

# fig1 = plt.figure(1)   It is easier not to use .hist for fig 1, but then the plot does not update when I change thresh.
# Why ?
