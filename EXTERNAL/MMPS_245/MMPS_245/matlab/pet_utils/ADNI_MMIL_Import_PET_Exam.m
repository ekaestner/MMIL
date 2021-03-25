function N = ADNI_MMIL_Import_PET_Exam(PETdir)

try
  N = Read_ADNI_PET_dir(PETdir);
catch
  fprintf(1,'Reading PET data from directory %s failed\n',PETdir);
end

