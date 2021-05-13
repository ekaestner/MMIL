%VisitID = 'FAMOUSFACESEPK183_1102_BITC20160825_';
%VisitID = 'FAMOUSFACESSubjectEPK206_BITC20151204_';
VisitID = 'FAMOUSFACESSUBJEPK173_BITC20141210';

RootDirs = [];
RootDirs.raw = '/space/mcdonald-syn01/1/data/MCD_MRI/raw/';
ContainerPath = MMIL_Get_Container(RootDirs,VisitID,'raw');
fname_classify = '/home/ekaestne/gitrep/MMIL/SCRATCH/2021/MCD_MRI_Series_Classify_EmoryInclude.csv';

errcode = MMIL_Classify_Dicoms_ejk_test(ContainerPath,[],fname_classify,1);

