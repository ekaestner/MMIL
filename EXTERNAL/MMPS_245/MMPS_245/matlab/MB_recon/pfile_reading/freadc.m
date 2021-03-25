function C = freadc (fid, len)
% Read and deblank character strings
  C = deblank( char( fread( fid, [1, len], 'uchar')));
