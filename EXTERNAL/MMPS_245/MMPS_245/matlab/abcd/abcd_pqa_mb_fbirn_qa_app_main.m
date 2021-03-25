function abcd_pqa_mb_fbirn_qa_app_main(input,output)
%fBIRN_QA_APP_MAIN Summary of this function goes here
%   Detailed explanation goes here


%================Check inputs===========% 
tic

if nargin <2
    fprintf('Requires input and output directory. Exiting...\n');
    return;
end

%================Read dicom files====================%

[vol,meta,fwhm] = abcd_pqa_read_files(input,output);
  
%==============Calls fBIRN QA routine=================%


abcd_pqa_mb_fbirn(vol, meta, output,fwhm); 


%====================================


fprintf('Finished\n');
toc
end

