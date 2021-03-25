function Transform_Tri_Atlas_to_Subject(input_trifile,output_trifile,regfile);
%function Transform_Tri_Atlas_to_Subject(input_trifile,output_trifile,regfile);
% Last Mod: 01/08/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,3)) return; end;

fprintf('%s: reading trifile...\n',mfilename);
surf = fs_read_trisurf(input_trifile);

fprintf('%s: loading regfile...\n',mfilename);
load(regfile);

fprintf('%s: transforming coordinates from atlas to subject...\n',mfilename);
surf.vertices = Transform_Coords_Atlas_to_Subject(surf.vertices,regStruct);

fprintf('%s: writing trifile...\n',mfilename);
fs_write_trisurf(surf,output_trifile);

fprintf('%s: finished.\n',mfilename);

