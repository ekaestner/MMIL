The functions in this folder are for p-file reading.




The following functions are from or modified from Brian Hargreaves' matlab script package 'geX'.
        rawheadX
        rawloadX
        getoffsetsX (modified)
        rawloadX_tseries(modified from rawloadX)
        load_raw_tseries(modified from rawloadX to work with RDS)
The 'geX' package is for reading/writing GE raw-data P-files. 
When reading p-file headers, the functions in the 'geX' package only reads in a few fields.

	


The following functions are from the GE company.
	fread_vartype
	freadc
	hdr_sizes
 	read_MR_headers (commented out some 'disp' functions and specified the directory for the hdr_sizes.txt file)
	read_data_acq_tab
	read_exam_header
	read_headerinfo	(commented out some 'disp' functions)
	read_image_header
	read_rdb_hdr
	read_series_header
These functions read in much more p-file header information than the 'geX' package. 
These functions are part of the GE code for multiplexed EPI reconstruction.


Kangrong Zhu,   Stanford University     Aug 2012

