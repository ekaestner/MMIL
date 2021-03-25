import sys

#-----------------------------------------------------------------------
cmnd_syntx_ok = True
optn = ''
fname_out = ''

if 3 <= len(sys.argv) and len(sys.argv) <= 4:
    fname = sys.argv[1]
    optn = sys.argv[2]
    
    if optn in ['Encoding', 'Summary', 'ExportFile', 'DateTime', 'DateTimeDiagnos', 'Diagnos']:
        if optn == 'ExportFile':
            if len(sys.argv) == 4:
                fname_out = sys.argv[3]
            else:
                cmnd_syntx_ok = False
    else:
        print('Wrong option')
        cmnd_syntx_ok = False

else:
    cmnd_syntx_ok = False

if not cmnd_syntx_ok:
    print('Reads an ASCII file containing an E-Prime spreadsheet. Extract contents while fixing encoding and format issues.')
    print('Returns diagnostic code, extracts particular information like experiment date and time, or writes a file with clean spreadsheet.')
    print('                                                                                     Octavio Ruiz. 2017jun05-sep26, oct17-nov01')
    print('Usage:')
    print('  ./eprime_spreadsheet_read.py')
    print('  ./eprime_spreadsheet_read.py  file  option')
    print('  ./eprime_spreadsheet_read.py  file  ExportFile  out_file (without .ext)\n')
    print('Options:')
    print('                    Prints this help')
    print('   Encoding         Returns encoding, file format, and row organization of input file')
    print('   Summary          As Encoding, plus a summary of column names and last row')
    print("   ExportFile       Save a tab-separated output file with clean contents of input file; we'll make extension = '.txt'")
    print('   DateTime         Returns ISO starting date and time of experiment reported in input file')
    print('   DateTimeDiagnos  As DateTime, plus a file diagnostic code (see below):')
    print('   Diagnos          Returns a diagnostic number, formed by summing:')
    print("      File found?                           0  -->  [Unable to find file]")
    print("      Encoding:                     1,      2  -->  ['utf-8','utf-16']")
    print("      Separator:                   10,     20  -->  [Tab, Comma]")
    print("      Quoted rows?                  0,    100  -->  [No, Yes]")
    print("      Start_time_info:              0,   1000  -->  [Found, Unable to extract]")
    print("      File read?                    0,  10000  -->  [Yes,   Unable to read]")
    print("      Filename contains '()':          100000  -->  [Yes]  (Reserved, calculated elsewhere)")
    print("      Experiment matches filename?    1000000  -->  [No]")
    print()

    sys.exit()
#-----------------------------------------------------------------------


# Acceptable extreme values
rows_n_min = 2
cols_n_min = 2


#-----------------------------------------------------------------------
# Determine file encoding (UTF-16 or -8)

def file_check_encoding( fname ):

    # Return variables
    f_encoding  = ''
    f_num_lines = None
    ienc = 0

    encd_optns_list = ['utf-8','utf-16']
    for i, encd in enumerate(encd_optns_list):
        try:
            ff = open(fname, encoding=encd)
            txt = ff.read()
            ff.close()
            num_rets = txt.count('\n')
            if num_rets >= rows_n_min:
                f_encoding = encd
                f_num_lines = num_rets
                ienc = i+1
        except:
            # File reading failed => this is not the right encoding option
            num_rets = -1
        # print('encd =', encd, ',  num_rets =', num_rets, ',  f_encoding =', f_encoding, ',  f_num_lines =', f_num_lines)

    return f_encoding, f_num_lines, ienc
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Find file encoding, and if it contains a spreadsheet as expected

def file_check( fname ):

    # Return variables
    f_encoding  = ''
    sep_tab     = None
    hoffs       = None   # header offset (rows)
    f_num_lines = None
    f_txt_lngth = None
    quoted_rows = None
    diagns      = float('nan')

    # Local variables
    ienc  = 0
    itab  = 0
    iquot = 0

    # Find file encoding (UTF-16 or -8)
    f_encoding, f_num_lines, ienc = file_check_encoding( fname )

    # If able to determine file encoding, read file and proceed with more tests
    if len(f_encoding) > 0:
        ff = open(fname, encoding=f_encoding)
        txt = ff.read()
        ff.close()
        f_txt_lngth = len(txt)

        # Tab- or comma-separated values?
        if txt.count('\t') > cols_n_min:
            sep_tab = True
            itab = 1
        else:
            sep_tab = False
            itab = 2

        # Rows encased by quotation marks (as we found in a few files)?

        ff = open(fname, encoding=f_encoding)
        # get a "typical" line, after the min. number of initial rows
        i_typ = int( (rows_n_min + f_num_lines) / 2 )
        # for i in range(0,rows_n_min):
        for i in range(0,i_typ):
            line = ff.readline()
            # seps_in_line_num = line.count('\t')
        seps_in_line_num = line.count('\t')
        ff.close()
        aa = line.rstrip('\n')    # If line has a carriage return at the end, as it should, remove it
        if aa.startswith('"') and aa.endswith('"'):
            quoted_rows = True
            iquot = 1
        else:
            quoted_rows = False
            iquot = 0

        # Find the first line in file that seems to contain the column names
        # (some files have comments in the first line or lines, instead of the column names).
        # use identified separator
        seps_in_line_num_typic = seps_in_line_num

        ff = open(fname, encoding=f_encoding)
        # for hoffs in range(0,rows_n_min):
        for hoffs in range(0, f_num_lines-1):
            line = ff.readline()
            if 0 < line.count('\t') and line.count('\t') >= seps_in_line_num_typic:
                break
        ff.close()

        # Calculate a "diagnostic" number summarizing file's econding and organization
        diagns = ienc + itab*10 + iquot*100
    else:
        diagns = 0

    return f_encoding, sep_tab, hoffs, f_num_lines, f_txt_lngth, quoted_rows, diagns
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Import pandas, avoiding warning
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

import pandas as pd
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
def sprdhst_read( fname, encd, offs, sep_tab ):
    if sep_tab:
        dat = pd.read_csv(fname, encoding=encd, header=offs, sep='\t')
    else:
        # Separator = comma
        dat = pd.read_csv(fname, encoding=encd, header=offs,)
    return dat
#-----------------------------------------------------------------------


# ============================================================================================================================
#                                                          Main
import io
import datetime
import dateutil.parser

pd.set_option('display.width', 512)

# ------------------------------------------------------------------------
# Find file, find and report file encoding and format

encoding, sep_tab, hoffs, lines_num, txt_lngth, quoted_rows, diagns  =  file_check( fname )

if diagns <= 0:
    if optn in ['Diagnos']:
        print( diagns )
    else:
        print('File not found:', fname)
    sys.exit()

if len(encoding) <= 0:
    diagns = diagns + 10000
    print('eprime_sprdsht_get.py:  Unable to read file ', fname)

if optn in ['Encoding', 'Summary', 'ExportFile']:
    print('encoding =', encoding, ',  sep_tab =', sep_tab, ',  hoffs =', hoffs,
          ',  lines_num =', lines_num, ',  txt_lngth =', txt_lngth, ',  quoted_rows =', quoted_rows, ',  diagns =', diagns)

if diagns >= 10000:
    print( diagns )
    sys.exit()
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Read spreadsheet into a pandas DataFrame
if quoted_rows:
    ff  = open(fname)
    txt = ''
    for line in ff:
        aa = line.rstrip('\n')
        bb = aa.strip('"')
        aa = bb + '\n'
        txt = txt + aa
    ff.close()
    # Store "de-quoted" text as a readable memory object, that can be read into a pd.DataFrame
    buf = io.StringIO( txt )
    dat = sprdhst_read( buf, encoding, hoffs, sep_tab )
else:
    dat = sprdhst_read( fname, encoding, hoffs, sep_tab )
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Check that file contains a known task and that this task corresponds to its identifier in the file name

# Possible task keys, and corresponding file-name identifiers
tsk_fid_dict = {'MID':    'MID',
                'Nback': ['NBACK','WM'],
                'SST':    'SST' }
exper = ''
exp_cont_name_match = False

# if 'ExperimentName' in dat.columns[0]:   # It can be 'ExperimentName' or '0-ExperimentName' or ...?...
#     # Find experiment reported in spreadsheet
#     ExpNameCol = dat.columns[0]
#     for tsk in tsk_fid_dict.keys():
#         if tsk in dat[ExpNameCol].unique()[0]:
#             exper = tsk
#             exp_fnameID_list = tsk_fid_dict[tsk]
#             break
#     # Check if experiment matches file name identifier
#     for exp_fn in exp_fnameID_list:
#         if fname.rfind(exp_fn) > (len(fname)-10):   # (look for task id near the end of file name)
#             exp_cont_name_match = True
#             break

if 'ExperimentName' in dat.columns[0]:   # It can be 'ExperimentName' or '0-ExperimentName' or ...?...
    # Find experiment reported in spreadsheet
    ExpNameCol = dat.columns[0]

    exp_in_file = dat[ExpNameCol].unique()[0]
    if isinstance(exp_in_file, str):
        for tsk in tsk_fid_dict.keys():
            if tsk in exp_in_file:
                exper = tsk
                exp_fnameID_list = tsk_fid_dict[tsk]
                break
        # Check if experiment matches file name identifier
        for exp_fn in exp_fnameID_list:
            if fname.rfind(exp_fn) > (len(fname)-10):   # (look for task id near the end of file name)
                exp_cont_name_match = True
                break
    else:
        exper = ''
        exp_cont_name_match = False

if exp_cont_name_match:
    pass
else:
    diagns = diagns + 1000000


if diagns >= 1000000  and  not optn in ['DateTime', 'DateTimeDiagnos', 'Diagnos']:
    if len(exper) <= 0:
        print('\nUnable to find experiment type in file contents\n')
    else:
        if not exp_cont_name_match:
            print('\nExperiment type in file does not match file name')

    print('encoding =', encoding, ',  sep_tab =', sep_tab, ',  hoffs =', hoffs,
        ',  lines_num =', lines_num, ',  txt_lngth =', txt_lngth, ',  quoted_rows =', quoted_rows, ',  diagns =', diagns)

# if optn in ['Encoding', 'Summary', 'ExportFile']:
#     if len(exper) <= 0:
#         print('\nUnable to find experiment type in file contents\n')
#     else:
#         if not exp_cont_name_match:
#             print('\nExperiment type in file does not match file name')

#     print('encoding =', encoding, ',  sep_tab =', sep_tab, ',  hoffs =', hoffs,
#           ',  lines_num =', lines_num, ',  txt_lngth =', txt_lngth, ',  quoted_rows =', quoted_rows, ',  diagns =', diagns)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
if optn == 'ExportFile':
    fxpos = fname_out.rfind('.')
    if fxpos > 0:
        if fxpos > (len(fname_out) - 5):   # The rightmost '.' is near the end of filename, and thus consistent with an extension
            fname_bas = fname_out[0:fxpos]
        else:
            fname_bas = fname_out
    else:
        if fxpos < 0:
            fname_bas = fname_out
        else:
            print('Error: invalid output file name')
            exit()

    if len(fname_bas) > 0:
        fname_out = fname_bas + '.txt'
    print( '\nWriting tab-separated file:', fname_out, '...' )

    dat.to_csv( fname_out, index=False, sep='\t' )

    print('done')
    sys.exit()
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
if optn == 'Summary':
    print('dat:  ndim =', dat.ndim,  ',  shape(rows,cols) =', dat.shape)
    print('\nTask reported in file =', exper, ',  matches file-name identifier:', exp_cont_name_match)
    print('\nCOLUMN NAMES:\n', dat.columns )
    print('\nLAST ROW:\n',  dat.iloc[-1:])
    print('\n')
    if 'GetReady.RTTime' in dat.columns:
        print( dat.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'GetReady.RTTime']] )
    elif 'Wait4Scanner.RTTime' in dat.columns:
        print( dat.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'Wait4Scanner.RTTime']] )
    else:
        print( dat.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime']] )
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Extract time information and task starting time
exp_date_str = dat.iloc[0]['SessionDate']
exp_time_str = dat.iloc[0]['SessionTime']
exp_datime = dateutil.parser.parse( exp_date_str + ' ' + exp_time_str )

delay = float('nan')
if 'GetReady.RTTime' in dat.columns:
    delay = dat.iloc[0]['GetReady.RTTime']
elif 'Wait4Scanner.RTTime' in dat.columns:
    delay = dat.iloc[0]['Wait4Scanner.RTTime']

if delay == delay:
    delay = int(delay / 1000)   # seconds
    exp_t0 = exp_datime + datetime.timedelta(seconds=delay)
else:
    # delay is nan, => unable to extract task starting time delay
    exp_t0 = exp_datime
    diagns = diagns + 1000

if optn == 'Summary':
    print('exp_datime =', exp_datime, ',  t_del =', delay, 's ,  exp_t0 =', exp_t0, '\n')
elif optn == 'DateTime':
    print( exp_t0 )
elif optn == 'DateTimeDiagnos':
    print( exp_t0, ',  ', diagns )
elif optn == 'Diagnos':
    print( diagns )
# ------------------------------------------------------------------------
# ============================================================================================================================
