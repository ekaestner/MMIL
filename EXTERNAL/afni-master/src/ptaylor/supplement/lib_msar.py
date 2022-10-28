
import sys, os
import lib_apqc_tcsh as lat
import afni_base as ab


#ver = '1.2'; date = 'June 14, 2019'
# + [PT] revamped: have new features like title and reflink
#
#ver = '1.3'; date = 'June 17, 2019'
# + [PT] work on textblock and images therein
#
#ver = '1.4'; date = 'June 19, 2019'
# + [PT] no more REFLINK block, is just put in on command line
# + [PT] also introducing a special case of SECTION+TEXTBLOCK called
#        INTRO (bc the shebang is always at top of script)
#
ver = '1.5'; date = 'June 19, 2019'
# + [PT] lots of tiny things updated/fixed/added.  Now in pretty good
#        working order
#
##########################################################################

# Default "mode" (or environs) is CODE; this is a list of *other*
# special modes, each of which is either a single line, or has a pair
# of start/end keywords; we go through finding the start of each mode
# based on keyword (or by being default CODE), and then know what to
# look for for its ending mode.  SPACE is now also a mode, because it
# is useful to have as filler.
OTHER_MODES = [ 'TITLE',       # title of script/demo
                'SECTION',     # can breakup demo into parts
                'TEXTBLOCK',   # basic description, just for RST
                'TEXTINTRO',   # special case of SECTION+TEXTBLOCK
                'SHEBANG',     # should always be first line of code
                'HCODE'        # script: normal code; RST: hidden+expandable
               ]

# --------------------------------------------------------------------------

def get_line_start_mode(x) :
    '''Input a string; if it is the start of a field, get the name of that
field; other, returns empty string'''

    ss = x.split()

    if not(ss):
        return 'SPACE'

    elif ss[0].__contains__("#:TITLE") :
        return 'TITLE'

    elif ss[0].__contains__("#:SECTION") :
        return 'SECTION'

    elif ss[0].__contains__('cat') and \
         x.__contains__('<<')  and \
         x.__contains__('TEXTBLOCK') :
        return 'TEXTBLOCK'
    
    elif ss[0].__contains__('cat') and \
         x.__contains__('<<')  and \
         x.__contains__('TEXTINTRO') :
        return 'TEXTINTRO'

    elif ss[0].__contains__("#!") :
        return 'SHEBANG'
    
    elif ss[0].__contains__('#:HIDE_ON') :
        return 'HCODE'

    else:
        return 'CODE'

# --------------------------------------------------------------------------

def find_mode_and_span(W, lstart, prev_mode=''):
    '''W is list of strings.

    lstart is the index of the line we start with.

    prev_mode is last/existing type of mode/block we were in (this
    affects how an empty space is interpreted, for example)

    '''

    X = W[lstart:]

    Nx    = len(X)
    lline = X[0]
    ss    = lline.split()

    start_mode = get_line_start_mode(lline)

    if start_mode == 'SPACE':
        tmode = start_mode
        tspan = 1

    elif (start_mode == 'CODE') or (start_mode == 'SHEBANG') :
        # empty lines not already in one of the OTHER_MODES should
        # just be code
        tmode = start_mode

        tspan = 0
        for jj in range(1, Nx):
            tt    = get_line_start_mode(X[jj])
            if OTHER_MODES.__contains__(tt) :
                tspan = jj
                break
        if not(tspan):
            # everything is code
            tspan = Nx

    elif start_mode == 'SECTION' :
        tmode = start_mode

        # we have to count empty spaces after this as part of the
        # mode, for accounting purposes
        tspan = 1
        for jj in range(1, Nx):
            tt = X[jj].split()
            if tt :
                tspan = jj
                break

    elif start_mode == 'TITLE' :
        tmode = start_mode

        # we have to count empty spaces after this as part of the
        # mode, for accounting purposes
        tspan = 1
        for jj in range(1, Nx):
            tt = X[jj].split()
            if tt :
                tspan = jj
                break

    elif (start_mode == 'TEXTBLOCK') or (start_mode == 'TEXTINTRO') :
        # these basically obey the same rules
        tmode = start_mode

        tspan = 0
        for jj in range(1, Nx):
            tt    = X[jj].split()
            if tt:
                if tt[0] == tmode :
                    tspan = jj+1 # bc ends with keyword
                    break
        if not(tspan) :
            print('''** Parse error:  can't find end of "cat <<{kw}"\n'''
                  '''   started at line {lnum}'''
                  '''   Need a "{kw}" somewhere.'''
                  ''.format(lnum=lstart, kw=tmode))
            sys.exit(2)

    elif start_mode == 'HCODE' :
        tmode = start_mode

        tspan = 0
        for jj in range(1, Nx):
            tt    = X[jj].split()
            if tt:
                if tt[0].__contains__("#:HIDE_OFF"):
                    tspan = jj+1 # +1 bc of keyword ending
                    break
        if not(tspan) :
            print('''** Parse error:  can't find end of "#:HIDE_ON"\n'''
                  '''   started at line {}.'''
                  '''   Need a "#:HIDE_OFF" somewhere.'''
                  ''.format(lstart))
            sys.exit(3)

    else:
        print(start_mode)

    return tmode, tspan

# -----------------------------------------------------------------------

def add_in_space( X, 
                  tstart,
                  tspan ):
    '''X is a list of strings.

    '''

    tscript = ''
    trst    = ''

    # Note the '-1', because there is a closing keyword
    Nx = tstart + tspan - 1 

    for ii in range(tstart+1, Nx):
        trst    += X[ii]
        tscript += X[ii] 

    return tscript, trst

# -----------------------------------------------------------------------

def add_in_textblock( X, 
                      tstart,
                      tspan,
                      iopts,
                      tmode = 'TEXTBLOCK' ):
    '''X is a list of strings.

    This is only for the RST file.

    tmode can also be TEXTINTRO, which is the same except for
    prepending a section heading for TEXTINTRO (and the output jumps
    to the top o' the queue in the RST page).

    '''

    tscript = ''
    trst    = ''
    
    Nx = tstart + tspan - 1

    ii = tstart+1
    while ii < Nx :

        x  = X[ii]
        ss = x.split()

        if not(ss):
            # empty case: just add line.
            trst += x
            ii   += 1
            
        elif ss[0].__contains__("#:IMAGE") :
            # IMAGE subblock must be contiguous-- empty line declares
            # end of it
            im_span = 1
            for jj in range(ii+1, Nx): # search for end (next empty line)
                y  = X[jj]
                uu = y.split()
                if not(uu) :
                    break
                else:
                    im_span+= 1
            imtxt = add_in_textblock_image( X,
                                            ii,
                                            im_span,
                                            iopts )
            #print(imtxt)
            trst+= imtxt
            ii+= im_span

        elif ss[0].__contains__("#:INCLUDE") :
            # INCLUDE subblock must be contiguous-- empty line declares
            # end of it
            im_span = 1
            for jj in range(ii+1, Nx): # search for end (next empty line)
                y  = X[jj]
                uu = y.split()
                if not(uu) :
                    break
                else:
                    im_span+= 1
            imtxt = add_in_textblock_include( X,
                                              ii,
                                              im_span,
                                              iopts )
            #print(imtxt)
            trst+= imtxt
            ii+= im_span
            
        else:
            # simple case, just add line
            trst += x
            ii   += 1

    return tscript, trst

# -----------------------------------------------------------------------

def add_in_textblock_include( X, 
                              tstart,
                              tspan,
                              iopts ):
    '''X is a list of strings.

    This is only for the RST file.  It can include any of the
    "include" directive phrases (start line, end line, etc.) from here:
    http://docutils.sourceforge.net/docs/ref/rst/directives.html#including-an-external-document-fragment

    The general format is:
    
    #:INCLUDE:  filename
        :lines:    1-10
        :language: none

    '''

    trst    = ''
    
    Nx = tstart + tspan

    fname     = ''
    oplist    = []

    # (opt) title is everything else include in this first line
    ss = X[tstart].split()
    if len(ss) == 2 :
        fname = ss[1]
    else:
        print('''** ERROR in TEXTBLOCK_INCLUDE: wrong number ({}) '''
              '''of entries'''.format(len(ss)))
        print('''   Main line must look like this, with just 1 filename:'''
              '''     #:INCLUDE:  filename''')
        sys.exit(5)

    for ii in range(tstart+1, Nx):
        ss = X[ii].strip()
        oplist.append( ss )

    # Now build table; some formatting necessary-- calc basename for
    # RST, bc name might include path
    inc_bname = os.path.basename( fname )
    trst += '''
.. literalinclude:: {}/{}
'''.format( iopts.subdir_rst, inc_bname )

    for op in oplist:
        trst+= 3*' ' + op + '\n'

    iopts.add_media( fname )

    return trst
    
# -----------------------------------------------------------------------

def add_in_textblock_image( X, 
                            tstart,
                            tspan,
                            iopts ):
    '''X is a list of strings.

    This is only for the RST file.

    The general format is:
    
    #:IMAGE:  Optional Title || and other col title!
        IMAGE1 IMAGE2
        IMAGE3 
    #:IMCAPTION: optional caption for any image(s)

    '''

    trst    = ''
    
    Nx = tstart + tspan

    imtitle   = []
    Nheader   = 0
    imcaption = []
    imlist    = []

    # (opt) title is everything else include in this first line
    ss = X[tstart].split()
    if len(ss) > 1 :
        tt =  ' '.join(ss[1:])
        imtitle = tt.split('||') # can have more than one title
        Nheader = 1
        Ntitle  = len(imtitle)

    # Now check rest for images, and/or captions; by construct, no
    # line can be empty.  We read everything as images until we get to
    # the IMCAPTION keyword, after which all is CAPTION
    lmode = 'IMAGES'
    for ii in range(tstart+1, Nx):
        ss = X[ii].split()

        if ss[0].__contains__("#:IMCAPTION") :
            lmode = 'CAPTION'
            imcaption.append(ss[1:])
        elif lmode == 'IMAGES' :
            imlist.append(ss)
        elif lmode == 'CAPTION' :
            imcaption.append(ss[1:])

    # Calc max dims of image matrix
    Ncol = 0
    Nrow = len( imlist )
    for rr in imlist:
        if len(rr) > Ncol:
            Ncol = len(rr)
    allwid = str(100 // Ncol) + ' '
    
    # Now build table;  some formatting necessary
    
    trst += '''
.. list-table:: 
   :header-rows: {}
   :widths: {}

'''.format(Nheader, Ncol * allwid)

    if imtitle :
        for cc in range(Ncol) :
            if cc :
                symb = ' '
            else:
                symb = '*'

            if cc >= Ntitle :
                trst+= '   {} - {}\n'.format(symb, ' ')
            else:
                
                trst+= '   {} - {}\n'.format(symb, imtitle[cc].strip())
            
    for rr in range(Nrow):
        for cc in range(Ncol) :
            if cc :
                symb = ' '
            else:
                symb = '*'

            if cc >= len(imlist[rr]) :
                trst+= '   {} -\n'.format(symb)
            elif imlist[rr][cc] == 'NULL':
                trst+= '   {} -\n'.format(symb)
            else:
                # have a list of these to copy later
                iopts.add_media( imlist[rr][cc] )
                # calc basename for RST, bc img name might include path
                img_bname = os.path.basename( imlist[rr][cc] )
                trst+= '''   {} - .. image:: {}/{}
          :width: 100%   
          :align: center\n'''.format( symb, iopts.subdir_rst,
                                      img_bname )
                
                
    return trst

# -----------------------------------------------------------------------

def copy_images_to_subdir(iopts):

    count_fail = 0
    list_fail = []
    
    for ii in range(iopts.nmedia):
        iname = iopts.list_media[ii]
        oname = iopts.subdir + '/' + os.path.basename(iopts.list_media[ii])
        a, b, c = ab.simple_shell_exec("\cp {} {}".format(iname, oname))
        if a :
            count_fail+= 1
            list_fail.append( iname )

    return count_fail, list_fail

# --------------------------------------------------------------------------

def add_in_code ( X, 
                  tstart,
                  tspan, 
                  tmode='CODE',
                  cmode='Tcsh',
                  add_shebang_rst='' ):
    '''X is a list of strings.

    At the moment, just these modes: cmode = {none|Tcsh}

    '''

    tscript = ''
    trst    = ''

    Nx = tstart + tspan
    N0 = tstart 

    if tmode == 'HCODE' :
        Nx-= 1 # because there is a closing keyword
        N0+= 1 # because there is an opening keyword

        trst+= '''

.. hidden-code-block:: {}
   :starthidden: True
   :label: - show code y/n -

'''.format(cmode)

    elif tmode == 'SHEBANG' :
        trst+= '' 

    else: # should just be code
        trst+= '''

.. code-block:: {}

'''.format(cmode)

    if add_shebang_rst :
        trst    += add_shebang_rst # should already be indented
                
    for ii in range(N0, Nx):

        # try to format ending \s nicely
        text = X[ii]
        if len(text) >= 2:
            if text[-2:] == '\\\n' :
                text = lineup_eol(text)
        
        tscript += text
        trst    += 3*' ' + text

    return tscript, trst

# -----------------------------------------------------------------------

def add_in_section( X, 
                    tstart,
                    tspan ):
    '''X should just be a list of a single string, but may generalize...'''

    tscript = ''
    trst    = ''

    lline = X[tstart]
    ss    = lline.split()
    text  = ' '.join(ss[1:])
    ltext = len(text)
    
    tscript+= lat.bannerize( text, fullwid=72 )
    
    trst   += text + "\n"
    trst   += (ltext  + 2 ) * '-' 

    # might have empty lines/padding
    for ii in range(1, tspan):
        jj = tstart + ii
        tscript += X[jj]
        trst    += X[jj]

    return tscript, trst 

# --------------------------------------------------------------------------

def create_text_title( X, 
                       tstart,
                       tspan ):
    '''X should just be a list of a single string, but may generalize...'''

    tscript = ''
    trst    = ''

    lline = X[tstart]
    ss    = lline.split()
    text  = ' '.join(ss[1:])
    ltext = len(text)

    tscript += '# '
    tscript += text + '\n'
    
    trst    += (ltext  + 2 ) * '*'  + "\n" 
    trst    += text + "\n"
    trst    += (ltext  + 2 ) * '*' 

    # might have empty lines/padding
    for ii in range(1, tspan):
        jj = tstart + ii
        tscript += X[jj]
        trst    += X[jj]

    return tscript, trst

# --------------------------------------------------------------------------

def interpret_MSAR_list( X, iopts ):

    N = len(X)

    USED_A_CODE_BLOCK = 0 # this lets us know where the shebang/top of
                          # orig code can be placed
    
    oscript_txt = ''
    orst_txt    = ''
    
    # header for RST
    orst_txt+= '''.. _{reflink}:

TO_BE_THE_TITLE

.. contents:: :local:

Introduction
-------------

**Download script:** :download:`{script} <{spath}/{script}>`

TO_BE_THE_INTRO


'''.format(reflink=iopts.reflink, script=iopts.oname_script,
           spath=iopts.subdir_rst)

    # header for RST
    oscript_txt+= '''TO_BE_THE_SHEBANG
'''
    
    ii      = 0
    tmode   = ''
    
    while ii < N :

        tmode, tspan = find_mode_and_span(X, ii, prev_mode=tmode)
        print("++ {:10s} range: [{:4d}, {:4d})".format(tmode,ii+1,ii+tspan+1))

        if tmode == 'SHEBANG' :
            # special case: calculate and save to prepend later to
            # the first usage of CODE or HIDDEN CODE laterz
            tscript_sheb, trst_sheb = add_in_code( X,
                                                   ii,
                                                   tspan,
                                                   tmode=tmode )
            oscript_txt = oscript_txt.replace( "TO_BE_THE_SHEBANG",
                                               tscript_sheb )
            ii         += tspan
        
        elif (tmode == 'CODE') or (tmode == 'HCODE') :

            # This silliness is because we may have to prepend the
            # shebang to a code block, and the hidden-code-block has
            # some prepending+indenting of its own.
            add_shebang_rst = ''
            if not(USED_A_CODE_BLOCK) :
                add_shebang_rst = trst_sheb
                USED_A_CODE_BLOCK = 1
                
            tscript, trst = add_in_code( X,
                                         ii,
                                         tspan,
                                         tmode=tmode,
                                         add_shebang_rst=add_shebang_rst )
            oscript_txt+= tscript
            orst_txt   += trst 
            ii         += tspan

        elif tmode == 'SECTION' :
            # SECTIONs are one line of text, but can span several
            # lines of whitespace
            tscript, trst = add_in_section( X,
                                            ii,
                                            tspan )
            oscript_txt+= tscript
            orst_txt   += trst 
            ii         += tspan

        elif tmode == 'TITLE' :
            # TITLEs are one line of text, but can span several
            # lines of whitespace
            # NOTE the special output here
            tscript, trst = create_text_title( X,
                                               ii,
                                               tspan )
            
            oscript_txt+= tscript
            orst_txt    = orst_txt.replace( "TO_BE_THE_TITLE", trst )
            ii         += tspan

        elif tmode == 'TEXTBLOCK' :
            tscript, trst = add_in_textblock( X,
                                              ii, 
                                              tspan,
                                              iopts )
            #oscript_txt+= tscript
            orst_txt   += trst 
            ii         += tspan

        elif tmode == 'TEXTINTRO' :
            tscript, trst = add_in_textblock( X,
                                              ii, 
                                              tspan,
                                              iopts,
                                              tmode = 'TEXTINTRO' )
            #oscript_txt+= tscript
            orst_txt    = orst_txt.replace( "TO_BE_THE_INTRO", trst )
            ii         += tspan
            
        elif tmode == 'SPACE' :
            tscript, trst = add_in_space( X,
                                          ii, 
                                          tspan )
            oscript_txt+= tscript
            orst_txt   += trst 
            ii         += tspan

        else:
            print("** ERROR. Unrecognized:")
            print(X[ii])
            sys.exit(2)

    # this goes here on the off chance that no Intro is provided--
    # replace the place holder with emptiness.
    orst_txt    = orst_txt.replace( "TO_BE_THE_INTRO", '' )
            
    return oscript_txt, orst_txt

# ----------------------------------------------------------------------

def lineup_eol( x, llen = 72 ) :

    Nx = len(x)
    diff = llen - Nx

    if diff > 0 :
        y = x[:-2] + diff * ' ' + x[-2:]
    else:
        y = x

    return y

# --------------------------------------------------------------------------

def read_text_to_list( fname ): 

    fff = open(fname, 'r')
    x = fff.readlines()
    fff.close()

    return x

# -------------------------------------------------------------------

def ARG_missing_arg(arg):
    print("** ERROR: missing argument after option flag: {}".format(arg))
    sys.exit(1)

# --------------------------------------------------------------------------

class init_opts_MSAR:

    script_type   = 'tcsh' # can generalize to other types later
    
    # req input
    infile        = ""
    prefix_rst    = ""    # can/should include path (a/b/FILE.rst)
    oname_script  = ""    # should just be name of file (SCRIPT.tcsh)
    reflink       = ""    # name of subdir in media/, and RST link name

    # opt input
    do_execute    = 0     # opt to run created script (e.g., to gen imgs)
    
    # made by finish_defs()
    outdir        = ""    # where everything will go
    meddir        = ""    # inside outdir, where media/ dir is
    subdir        = ""    # inside meddir, where images/script will go
    prefix_script = ""    # full name for outputting script
    subdir_rst    = ""    # local dir for RST path

    # made later during parsing, potentially
    list_media    = []    # list of images to copy into media/reflink/.
    nmedia        = 0
    
    # ----------- req -----------------

    def set_infile(self, fff):
        self.infile = fff

    def set_prefix_rst(self, prefix):
        self.prefix_rst = prefix
        
    def set_oname_script(self, oname):
        self.oname_script = oname
        
    def set_reflink(self, reflink):
        self.reflink = reflink

    # ----------- opt -----------------
    
    def set_execute(self):
        self.do_execute = 1 

    def add_media(self, X):
        self.list_media.append( X )
        self.nmedia+= 1
        
    # ------- complete some var defs/names --------

    def finish_defs(self):
        pp = os.path.dirname(self.prefix_rst)
        if not(pp) :  
            pp = "."  # bc this means "here"

        # these are all potentially full/relative paths
        self.outdir = pp
        self.meddir = '/'.join([self.outdir, 'media'])
        self.subdir = '/'.join([self.meddir, self.reflink])
        self.prefix_script = '/'.join([self.subdir, self.oname_script])
        
        # just local/partial path within RST dir
        self.subdir_rst = '/'.join(['media', self.reflink])

    # ---------- check ----------------

    def check_req(self):
        MISS = 0
        if not(self.infile) :
            print("** missing: infiles")
            MISS+=1
            
        if self.prefix_rst == "" :
            print("** missing: prefix_rst")
            MISS+=1
        elif self.prefix_rst[-4:] != '.rst' :
            print("** '-prefix_rst ..' entry MUST end with '.rst'")
            MISS+=1
            
        if self.reflink == "" :
            print("** missing: reflink")
            MISS+=1
        elif self.reflink[0] == "_" :
            print("** don't start the reflink name with an underscore '_'!")
            MISS+=1

        if self.oname_script == "" :
            print("** missing: prefix_script")
            MISS+=1
        elif self.oname_script.__contains__("/") :
            print("** don't include path in the '-oname_script ..' entry")
            MISS+=1

        return MISS


# --------------------------------------------------------------------------

help_string_MSAR = '''

PURPOSE ~1~ 

Program to take a script with some special (~simple) markup and turn
it into both an RST page and a script for the online Sphinx
documentation.


INPUTS ~1~ 

-prefix_rst AA    :(req) output filename, including any path, of the
                   RST/Sphinx file.  AA must include file extension 
                   '.rst'.  E.g.:  tutorial/fun_3dcalc.rst

-prefix_script BB :(req) output filename, *without* any path, of the
                   script file.  BB probably should include file extension, 
                   such as '.tcsh'.  E.g.:  fun_3dcalc.tcsh


-reflink       CC :(req) a string tag that will be 1) subdirectory name
                   holding images for the given demo, and 2) the RST 
                   internal reference label, as '.. _CC:'.  First character
                   of CC must be alphabetic.

-execute_script   :(req/opt) flag to not just create the RST+script, but
                   to execute the script as well.  IF the script
                   generates images that will be copied to the
                   media/CC/. directory, then this flag should be used
                   at least the first time the script is run (so the
                   files can be copied); it may not be necessary to
                   execute on later runs.
 

OUTPUTS ~1~ 

+ an RST file, which is basically a Sphinx-formatted page, that can be
  placed in a separate directory

+ an output directory to put into the Sphinx tree, called
  [rst-path]/media/CC, where [rst-path] is the location of the output
  RST file and CC is the reflink name.

+ a script file, both locally (where the script is run, so that it can
  be executed) and in [rst-path]/media/CC (which will be shown in the
  RST pages).

+ images made by the script which are flagged to be show in the RST
  pages will be copied to [rst-path]/media/CC/. 


EXAMPLES ~1~

   1) First time through, execute script to make images:
   adjunct_make_script_and_rst.py                                          \
       -input          ex_afni11_roi_cmds.tcsh                             \
       -reflink        afni11_roi_cmds                                     \
       -prefix_script  afni11_roi_cmds.tcsh                                \
       -prefix_rst ~/afni_doc/tutorials/rois_corr_vis/afni11_roi_cmds.rst  \
       -execute_script

   2) Second time through, if "only" text changes/formatting:
   adjunct_make_script_and_rst.py                                          \
       -input          ex_afni11_roi_cmds.tcsh                             \
       -reflink        afni11_roi_cmds                                     \
       -prefix_script  afni11_roi_cmds.tcsh                                \
       -prefix_rst ~/afni_doc/tutorials/rois_corr_vis/afni11_roi_cmds.rst 

'''

def parse_MSAR_args(argv):

    Narg = len(argv)

    if not(Narg):
        print(help_string_MSAR)
        sys.exit(0)

    # initialize objs
    iopts  = init_opts_MSAR()    # input-specific

    # check through inputs
    i = 0
    while i < Narg:
        if argv[i] == "-ver":
            print(ver)
            sys.exit(0)

        elif argv[i] == "-help" or argv[i] == "-h":
            print(help_string_apqc_1dplot)
            sys.exit(0)

        # ---------- req ---------------

        elif argv[i] == "-input":
            if i >= Narg:
                ARG_missing_arg(argv[i])
            i+= 1
            iopts.set_infile(argv[i])
        
        elif argv[i] == "-prefix_rst":
            if i >= Narg:
                ARG_missing_arg(argv[i])
            i+= 1
            iopts.set_prefix_rst(argv[i])
            
        elif argv[i] == "-prefix_script":
            if i >= Narg:
                ARG_missing_arg(argv[i])
            i+= 1
            iopts.set_oname_script(argv[i])

        elif argv[i] == "-reflink":
            if i >= Narg:
                ARG_missing_arg(argv[i])
            i+= 1
            iopts.set_reflink(argv[i])

        # ---------- req ---------------

        elif argv[i] == "-execute_script":
            iopts.set_execute()

        # --------- finish -------------

        else:
            print("** ERROR: unknown opt: '{}'".format(argv[i]))
            sys.exit(2)
        i+=1

               
    if iopts.check_req():
        print("** ERROR with input arguments")
        sys.exit(1)

    iopts.finish_defs()  # make some paths we need later

    return iopts

