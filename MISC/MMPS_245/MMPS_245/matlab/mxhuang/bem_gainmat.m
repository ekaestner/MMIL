function gain_mat=bem_gainmat(Rq_bemgrid,Xfer,Verbose)
%BEM_GAINMAT - Computes the EEG/MEG forward/gain matrix associated with a set of grid points 
% function gain_mat=bem_gain(Rq_bemgrid,Xfer, Verbose)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EEG/MEG BEM 3-D  GRID FORWARD GAIN CALCULATION UTILITY 
%
% This function computes the EEG/MEG forward matrix (kernel) associated with a set of 
% grid points.
%
% INPUTS (Required):
%          Rq_bemgrid: a structure of specific source locations (Nx3 array)
%
%    Xfer: a structure variable return from bem_transfer
%           Xfer.basis_opt: string --- constant bem or linear bem
%            Xfer.test_opt: string -- collocation, glerkin
%                Xfer.mode: mode type
%            Xfer.geometry: geometry matrix
%               Xfer.nodes: nodes matrix
%                 Xfer.cdv: conductivity vector
%                  Xfer.Te: the EEG transfermation matrix
%                  Xfer.Tm: the MEG transformation matrix
%              Xfer.Te_ISA: the EEG transfermation matrix with ISA
%              Xfer.Tm_ISA: the MEG transformation matrix with ISA
%           Xfer.eegsensor: eegsensor
%               Xfer.R_eeg: R_eeg
%              Xfer.P_sens: P_sens
%               Xfer.R_meg: the MEG sensor locations, no_sensor by 3
%               Xfer.O_meg: the MEG sensor orientations, no_sensor by 3
%
% OUTPUTS: Complete set of output parameters are stored in user specified 
%          "*.mat" file defined by "bem_gaingrid_mfname.mat" 
%          gain_matrix.eeg: the eeg gain matrix (N_eeg_sensor by 3*N_dip)
%                           EEG unit: uV/nAm
%          gain_matrix.meg: the meg gain matrix (N_meg_sensor by 3*N_dip)
%                           MEG unit: fT/nAm
%
%   Other Notes: 3x3 MEG kernel is dotted with sensor orientation to generate more compact
%                1x3 kernel 
%
% By M.X. Huang, PhD, March 2004 
%    
%
%

scale_meg=0.1; % scaling factor converting MEG fields into fT by nAm dipole: 1e-7(mu_0/4pi) x 1e-9 (nAm) x 1e15(fT)
scale_eeg=0.001; % scaling factor converting EEG into uV by nAm dipole: 1e-9 (nAm) x 1e6 (uV)  

if Xfer.ISA==0
    if Xfer.mode==1
        Te_save = Xfer.Te;
        Tm_save = [];  % Dummy value for kernel routine
    elseif Xfer.mode==2
        Tm_save = Xfer.Tm;
        Te_save = [];  % Dummy value for kernel routine
    else
        Tm_save = Xfer.Tm;
        Te_save = Xfer.Te;       
    end
elseif Xfer.ISA==1
    if Xfer.mode==1
        Te_save = Xfer.Te_ISA;
        Tm_save = [];
    elseif Xfer.mode==2
        Tm_save = Xfer.Tm_ISA;
        Te_save = [];
    else
        Tm_save = Xfer.Tm_ISA;
        Te_save = Xfer.Te_ISA;               
    end
else
    if Verbose
        sprintf('Invalid ISA option chosen!!!')
    end
    
end
%
clear Xfer.Te Xfer.Tm Xfer.Te_ISA Xfer.Tm_ISA
%
%%%% THIS PART DETERMINES DIMENSION PARAMETERS AND DIVIDES THE INPUT DIPOLE SET INTO A SET OF BLOCKS %%%% 
%
if Xfer.mode==1  %EEG only
    Meeg = size(Xfer.R_eeg,1);      % Number of Sensors
    S = [];
    R_meg = [];
    O_meg = [];
elseif Xfer.mode==2  %MEG only
    Mmeg = size(Xfer.R_meg,1);
    Xfer.R_eeg = [];
    Xfer.P_wts = [];    
else
    Meeg = size(Xfer.R_eeg,1);      % Number of Sensors
    Mmeg = size(Xfer.R_meg,1);      % Number of Sensors    
end
%
P = size(Rq_bemgrid,1);   % Number of Dipole Locations
%
BLK_SIZE = 50;  % Maximum Number of dipoles to process at a single time
% Set to prevent overstepping memory bounds in Kernel Calculation
blk_tot = floor(P/BLK_SIZE);  % Total Number of blocks
blk_rem = rem(P,BLK_SIZE);    % Residual Number of Dipoles
%
%%%% THIS PART COMPUTES THE GAIN MATRIX KERNEL (PROCESSING A BLOCK AT A TIME) %%%%
%
t0 = clock;  % Start the clock a tickin...
if Verbose
    hwait = waitbar(0,'Creating the gain matrix for selected sources. . .');
end

%

if Xfer.mode==1    % EEG
    Reeg = Xfer.eegsensor(:,2:4);
    Rmeg = [];
    Gmeg_grid = [];
    Geeg_grid = zeros(Meeg,3*P);   % Pre-initialize Output
    for nblk = 1:blk_tot
        n1 = (nblk-1)*3*BLK_SIZE + 1;  % equivalent dipole start index
        n2 = nblk*3*BLK_SIZE;          % equivalent dipole end index
        n3 = (nblk-1)*BLK_SIZE + 1;    % dipole start index
        n4 = nblk*BLK_SIZE;            % dipole end index
        [Geeg_grid(:,n1:n2) dummy] = bem_kernel(Rq_bemgrid(n3:n4,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
        telap = etime(clock,t0);
        trem = telap*(P-n4)/n4;
        if Verbose & ~rem(nblk,10)
            %bst_message_window('overwrite', sprintf('Completed %.0f of %.0f grid points in %3.2f Sec.',n4,P,telap))
            %Expected time to complete: %.2f Sec....',n4,P,telap,trem))
            waitbar(nblk/blk_tot, hwait);
        end
        
    end
    if blk_rem > 0
        n1 = blk_tot*3*BLK_SIZE + 1;   % equivalent dipole start index (for residual dipoles)
        n2 = 3*P;                      % equivalent dipole end index (for residual dipoles)
        n3 = blk_tot*BLK_SIZE + 1;     % dipole start index
        n4 = P;                        % dipole end index
        [Geeg_grid(:,n1:n2) dummy] = bem_kernel(Rq_bemgrid(n3:P,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
    end
    
elseif Xfer.mode==2   % MEG only
    
    Reeg = [];
    Geeg_grid=[];
    Gmeg_grid = zeros(Mmeg,3*P);       % Pre-initialize Output
    Gtemp = zeros(Mmeg,3*BLK_SIZE);
    for nblk = 1:blk_tot
        n1 = (nblk-1)*3*BLK_SIZE + 1;  % equivalent dipole start index
        n2 = nblk*3*BLK_SIZE;          % equivalent dipole end index
        n3 = (nblk-1)*BLK_SIZE + 1;    % dipole start index
        n4 = nblk*BLK_SIZE;            % dipole end index
        Ptemp = BLK_SIZE;              % number dipoles being processed at current time
        [dummy Gtemp] = bem_kernel(Rq_bemgrid(n3:n4,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
        %
        for m=1:Mmeg
            Gmeg_grid(m,n1:n2) = Xfer.O_meg(m,:)*Gtemp(3*m-2:3*m,:);  % Apply Sensor Orient to each lead field
        end
        %                                   
        telap = etime(clock,t0);
        trem = telap*(P-n4)/n4;
        if Verbose & ~rem(nblk,100)
            %bst_message_window('overwrite', sprintf('Completed %.0f of %.0f grid points in %3.1f Sec.',n4,P,telap));
            %Expected time to complete: %.2f Sec....',n4,P,telap,trem))
            waitbar(nblk/blk_tot, hwait);
        end
        
    end
    if blk_rem > 0
        n1 = blk_tot*3*BLK_SIZE + 1;   % equivalent dipole start index (for residual dipoles)
        n2 = 3*P;                      % equivalent dipole end index (for residual dipoles)
        n3 = blk_tot*BLK_SIZE + 1;     % dipole start index
        n4 = P;                        % dipole end index
        Ptemp = n4-n3+1;               % number dipoles being processed at current time
        Gtemp = zeros(Mmeg,3*Ptemp);
        [dummy, Gtemp] = bem_kernel(Rq_bemgrid(n3:P,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
        %
        for m=1:Mmeg 
            Gmeg_grid(m,n1:n2) = Xfer.O_meg(m,:)*Gtemp(3*m-2:3*m,:);   % Apply Sensor Orient to each lead field
        end
    end
    
 elseif Xfer.mode==3   % both EEG and MEG
    Reeg = Xfer.eegsensor(:,2:4);
    Geeg_grid = zeros(Meeg,3*P);   % Pre-initialize Output
    Gmeg_grid = zeros(Mmeg,3*P);       % Pre-initialize Output
    Gtemp = zeros(Mmeg,3*BLK_SIZE);
    for nblk = 1:blk_tot
        n1 = (nblk-1)*3*BLK_SIZE + 1;  % equivalent dipole start index
        n2 = nblk*3*BLK_SIZE;          % equivalent dipole end index
        n3 = (nblk-1)*BLK_SIZE + 1;    % dipole start index
        n4 = nblk*BLK_SIZE;            % dipole end index
        Ptemp = BLK_SIZE;              % number dipoles being processed at current time
        [Geeg_grid(:,n1:n2) Gtemp] = bem_kernel(Rq_bemgrid(n3:n4,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
        %
        for m=1:Mmeg
            Gmeg_grid(m,n1:n2) = Xfer.O_meg(m,:)*Gtemp(3*m-2:3*m,:);  % Apply Sensor Orient to each lead field
        end
        %                                   
        telap = etime(clock,t0);
        trem = telap*(P-n4)/n4;
        if Verbose & ~rem(nblk,100)
            %bst_message_window('overwrite', sprintf('Completed %.0f of %.0f grid points in %3.1f Sec.',n4,P,telap));
            %Expected time to complete: %.2f Sec....',n4,P,telap,trem))
            waitbar(nblk/blk_tot, hwait);
        end
        
    end
    if blk_rem > 0
        n1 = blk_tot*3*BLK_SIZE + 1;   % equivalent dipole start index (for residual dipoles)
        n2 = 3*P;                      % equivalent dipole end index (for residual dipoles)
        n3 = blk_tot*BLK_SIZE + 1;     % dipole start index
        n4 = P;                        % dipole end index
        Ptemp = n4-n3+1;               % number dipoles being processed at current time
        Gtemp = zeros(Mmeg,3*Ptemp);
        [Geeg_grid(:,n1:n2), Gtemp] = bem_kernel(Rq_bemgrid(n3:P,:),Xfer.basis_opt,Xfer.test_opt,Xfer.geometry,Xfer.nodes,Xfer.cdv,Te_save,Xfer.P_wts,Tm_save,Xfer.R_meg);
        %
        for m=1:Mmeg 
            Gmeg_grid(m,n1:n2) = Xfer.O_meg(m,:)*Gtemp(3*m-2:3*m,:);   % Apply Sensor Orient to each lead field
        end
    end   
end

clear Gtemp

% % If required, apply source orientation
% if isstruct(bem_grid_mfname) 
%     if ~isempty(bem_grid_mfname.Orient) % Orientations are actually specified
%         if Verbose & ~isempty(bem_grid_mfname.Orient)
%             bst_message_window('Applying source orientation. . .')  
%             isrc = 0;
%             Gtmp = GBEM_grid; 
%             clear GBEM_grid
%             for src = 1:size(Gtmp,2)/3
%                 isrc = isrc + 1;
%                 GBEM_grid(:,src) = Gtmp(:, 3*(isrc-1)+1:3*isrc) * bem_grid_mfname.Orient(:,src);
%             end
%         end
%     end
% end

telap_fwdgaingrid = etime(clock, t0);
if Verbose
    close(hwait)
    bst_message_window(sprintf('Total time : %3.1f sec.',telap_fwdgaingrid))
end


%
%%%% THIS PART STORES THE FINAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
gain_mat.eeg=Geeg_grid*scale_eeg;  % convert EEG into uV by nAm dipole 
gain_mat.meg=Gmeg_grid*scale_meg;  % convert MEG into fT by nAm dipole

clear Geeg_grid Gmeg_grid




