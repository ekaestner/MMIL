function [dim, inm] = utilmesh_check_source_space(so,C,E);
% so is the source space
% C,E is the linear mesh
% dim gives the vector of distances
% inm gives the vector of inside sources (logical)
%
% NOTE: this function is a modified version of that from the NFT toolbox
%       changed to not display a waitbar
%
% copied:   09/27/13 from NFT toolbox (NFT-2.3_64bit)
% last mod: 10/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hh = waitbar(0,'computing...');
M = size(so,1);
for i = 1 : M
%    waitbar(i/M)
    [dm, Pm, el, in] = utilmesh_dist_mesh_point(so(i,:), C, E);
    dim(i) = dm;
    inm(i) = in;
end
%close(hh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

