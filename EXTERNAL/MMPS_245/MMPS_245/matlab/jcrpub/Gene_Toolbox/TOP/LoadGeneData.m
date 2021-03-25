function [SNP,Pheno,SubjInfo] = LoadGeneData(gene_name,pheno_class,indatdir,sexchkfile)

%  LoadGeneData.m
%
%  Usage: [SNP,Pheno,SubjInfo] = LoadGeneData(gene_name,pheno_class,indatdir)
%
%  LoadGeneData loads into memory 3 data structs necessary for doing per gene SNP/phenotype 
%   association tests.
%
%  Inputs: 
%     gene_name     - char string (e.g. 'MECP2')
%     pheno_class   - 'ROI' (morphometry) or 'COG' (cognitive measures)
%     indatdir      - input data directory name
%     sexchkfile    - sexcheck results file
%
%  Outputs:
%     SNP     - struct w/ the following fields:  .dat  .snp_names  .subjID
%     Pheno   - struct w/ the following fields:  .dat  .names  .subjID         
%     SubInfo - struct w/ the following fields:  .ID .age .sex .eth .diagn

snpfile = sprintf('%s/SNP_table_%s.csv',indatdir,gene);
if ~exist(snpfile,'file')
  error(sprintf('Input file %s does not exist.',snpfile));
end

pheno_file = sprintf('%s/%s_table.csv',indatdir,pheno_class);
if ~exist(pheno_file,'file')
  error(sprintf('Input file %s does not exist.',pheno_file));
end

[SNPcell,filestats] = readtext_amd(snpfile,',','#','');

[phenocell,filestats2] = readtext_amd(pheno_file,',','#','');

SNP_mat = str2double(SNPcell(2:end,2:end));
SNP_subjID = str2double(rmquotes(SNPcell(2:end,1)));
pheno_subjID = str2double(rmquotes(phenocell(2:end,1)));

if ~isequal(SNP_subjID,pheno_subjID)
   error('Mismatched subject_IDs detected.');
end
subjID = SNP_subjID;
clear SNP_subjID pheno_subjID;

% SNP rs# names
rsSNP_names = rmquotes(SNPcell(1,2:end));

% Load phenotype data into numeric matrix
phenocell_colnames = rmquotes(phenocell(1,:));
if strcmp(pheno_class,'ROI')
   % Get data for our original 4 phenotypes
   pheno_names = [{'ICV'}, {'BrainVol'}, {'CorticalThickness'}, {'CorticalArea'}];
   for ii=1:length(pheno_names),
     pheno_col(ii) = strmatch(pheno_names{ii},phenocell_colnames);
   end
   pheno = str2double(phenocell(2:end,pheno_col));
else
   % Use all phenotypes
   pheno_names = phenocell_colnames(6:end);
   pheno_col = [6:length(phenocell_colnames)];
end
pheno = str2double(phenocell(2:end,pheno_col));

% Get demographic data
agecol = strmatch('Age',phenocell_colnames,'exact');
sexcol = strmatch('Sex',phenocell_colnames,'exact');
ethcol = strmatch('Ethnicity',phenocell_colnames,'exact');
diagcol = strmatch('Diag',phenocell_colnames, 'exact');
age = str2double(phenocell(2:end,agecol));
sex = str2double(phenocell(2:end,sexcol));
eth = str2double(phenocell(2:end,ethcol));
diagn = str2double(phenocell(2:end,diagcol));

% Remove subjects w/ mismatched sex labels, non-Scandinavian subjects, and
% subjects w/ missing demographic data.
[sex_geno_mismatch] = TOP_checksex(sexchkfile);
demog_NaN = sum(isnan([age sex eth diagn]),2);
keep_ind = [~ismember(subjID, sex_geno_mismatch) & eq(eth,1) & ~demog_NaN];
SNP_mat = SNP_mat(keep_ind,:);
subjID = subjID(keep_ind);
pheno = pheno(keep_ind,:);
sex = sex(keep_ind);
age = age(keep_ind);
eth = eth(keep_ind);
diagn = diagn(keep_ind);

% Reorganize data here into 3 structures: SNP, Pheno, and SubjInfo
SNP.dat = SNP_mat;
SNP.snp_names = rsSNP_names;
SNP.subjID = subjID;

Pheno.dat = pheno;
Pheno.names = pheno_names;
Pheno.subjID = subjID;

SubjInfo.ID = subjID;
SubjInfo.age = age;
SubjInfo.sex = sex;
SubjInfo.eth = eth;
SubjInfo.diagn = diagn;
