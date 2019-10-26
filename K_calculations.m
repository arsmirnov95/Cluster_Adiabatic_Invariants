%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                         
%       This is a MATLAB routine to compute the second adiabatic        
%       invariant K. The routine employs Kp index as geomagnetic        
%       input. Computations are based on Tsyganenko-89 model.         
%                                                                                                         
%       To run the code, IRBEM library should be installed on the          
%       computer.                                                                                  
%                                                                                                         
%                                                                                                         
%       Generated by: A.G. Smirnov, E.A. Kronberg, P.W.Daly, N.A.  Aseev, Y.Y. Shprits and A. Kellerman.         
%                                                                                                         
%       In case of any questions contact Artem Smirnov via                 
%       arsmirnov95@gmail.com                                                          
%                                                                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input variables : Kp index;
% x_GSE - satellite position in GSE coordinates (in R_E); 
% t_GSE - timetags for position (matlab date format).
% 

kext      = 4;
options   = [0 0 0 0 0];
sysaxes   = 3; % GSE coordinate system

maginput  = [Kp_interp, zeros(length(t_GSE),24)];
local_PA    = 10:20:90; % local pitch angles. It is sufficient to specify just 1 branch,
%because K is symmetric w.r.t. 90 degrees.

[~,~,Bmirror,~,J,~] = onera_desp_lib_make_lstar_shell_splitting(kext, options, sysaxes, t_GSE, x_GSE(:,1), x_GSE(:,2), x_GSE(:,3), local_PA, maginput);

K = J.*sqrt(Bmirror); % in [ nT^(1/2) R_E].
