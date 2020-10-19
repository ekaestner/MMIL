function return_var = bem_transfer_3shell(bem_input,NVertMax,H_inv_filename,Verbose);

% function return_var = bem_transfer_3shell(bem_input,NVertMax,H_inv_filename,Verbose);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EEG/MEG FORWARD MODEL USING BOUNDARY ELEMENT METHODS
%                      - TRANSFER MATRIX CALCULATION (bem_transfer_1shell.m)
%
% This function computes the "transfer matrix" associated with the voltage
% potential  forward gain matrix for an array of EEG electrodes located on the
% outermost layer of a single layer surface,
%                              -AND/OR-
% the radial magnetic field forward gain matrix for an array of MEG sensors
% located on the  outermost layer of a single layer surface.
% 
% Surface data for each layer is defined via a user specified tesselated grid.
% Each region  of the multilayer surface is assumed to be concentric with
% isotropic conductivity.   EEG sensors are assumed to be located on the surface
% of the single layer at one of the surface tesselation points. Should
% specified EEG sensor points not coincide with a point on the tesselated
% outer surface layer grid, the sensor location will be quantized to the
% closest outer surface tesselation point. Dipole generator(s) are assumed to be
% interior  to the single layer surface
%
% By M.X. Huang, PhD 
% 
% Last Update Feb 9, 2007 -- rename dotprod with dotprod_fun, crossprod
% with crossprod_fun to avoid the conflict with MATLAB files of the same
% names
% 
% 
% INPUTS (Required)
%       bem_input.R_eeg  : EEG sensor locations on scalp (meters)         (Meeg x 3)
%                         (insert dummy value if mode = 2)
%       bem_input.R_meg  : MEG sensor locations on scalp (meters)         (Mmeg x 3)
%                         (insert dummy value if mode = 1)
%       bem_input.O_meg  : MEG sensor orientations (unit-vector)          (Mmeg x 3)
%                         (insert dummy value if mode = 1)
%    bem_input.vertices  : Tessalated Surface layers from INNERMOST TO OUTERMOST. 
%                          Packed as Cell Array where index 1 is innermost surface  (1 x NL)-Cell Array
%                          Individual entries are position row vectors in PCS
%       bem_input.faces  : Tesselated Surface Node connections for each triangle    (1 x NL)-Cell Array
%                          Packed as Cell Array where index 1 is innermost surface
%                          Individual entries are indexes to "vertices"
%       bem_input.sigma  : conductivity from INNERMOST to OUTERMOST        (1 x NL)
%       bem_input.mode   : 1,  compute EEG only
%                          2,  compute MEG only
%                          3,  compute both EEG and MEG
%   bem_input.basis_opt  : 0, constant basis
%                          1, linear
%   bem_input.test_opt   : 0, collocation,
%                          1, Galerkin
%     bem_input.ISA      : 0,  Inhibit Isolated Skull Approach
%                        : 1,  Enable Isolated Skull Approach
%    bem_input.fn_eeg    : EEG transfer matrix filename                   Char String
%                         - EEG Transfer matrices are computed and then saved 
%                          in this "*.mat" file
%
%    bem_input.fn_meg    : MEG transfer matrix filename                   Char String
%                         - MEG Transfer matrices are computed and then saved 
%                         in this "*.mat" file
%
%         WHERE: M = # of sensors; P = # of dipoles;
%               Nn = # of Surface Tesselation Nodes; Nt = # of Surf Tess Triangles 
%
%       H_inv_filename     : file name containing previously calculated H_inv matrices
%                         If file does not exist, H_inv will be calculated and 
%                         the newly calculated H_inv will be saved in a file with this
%                         name.
%
% INPUTS (Optional):
%    NVertMax : Maximum Number of vertices for Surface Tesselation  (for a single layer)
%                   - Surface is re-tesselated if # points exceeds max      
%                   - Default Maximum Number of Triangle Faces is 2292    (1 x NL)
%
%    Verbose : Toggles Verbose mode on/off (on = 1);
%
% OUTPUTS: 
%    return_var: structure variable
%           return_var.basis_opt: string --- constant bem or linear bem
%           return_var.test_opt: string -- collocation, glerkin
%           return_var.mode: mode type, 1 for EEG, 2 for MEG, 3 for both
%           return_var.geometry: geometry matrix
%           return_var.nodes: nodes matrix
%           return_var.cdv: conductivity vector
%           return_var.ISA: 1 for ISA, 0 for no ISA
%           return_var.Te: the EEG transfermation matrix
%           return_var.Tm: the MEG transformation matrix
%           return_var.Te_ISA: the EEG transfermation matrix with ISA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Last modification: M.X. Huang, July 2006
% Adding and modifying return_var to the output, M. Huang
%
% Other people involved in old version: J. Chang, J. Mosher, E. Ermer
