function [dx] = KalmanFilterIterativeStandard(Adjust, x_pred)
% Function for Kalman Filter with inner-epoch iteration. 
% The change of the parameters is estimated as pseudo-observation trough 
% an inner-epoch iteration. Start of the iteration is the zero-vector as 
% the change of the parameters is expected to be zero.
% Check [01]: p.244, formulas without Gain Matrix
% 
% INPUT: 	
%   Adjust  containing all adjustment relevant data [struct]	 
%  	x_pred      predicted parameters [| vector]
% OUTPUT:     
%   struct x with   
%       x.x        	vector of adjusted parameters [| vector]
%       x.sigma2   	empirical variance of parameters
%       x.l      	vector of adjusted observations [| vector]
%       x.v     	residual vector [| vector]
%       x.r         redundancy [scalar]
%       x.Qxx       Cofactor Matrix of Parameters
%       x.Sxx       Covariance Matrix of Parameters
%       x.Qvv       Cofactor Matrix of Residuals
%       x.Svv       Covariance Matrix of Residuals
%
% This function belongs to raPPPid, Copyright (c) 2023, M.F. Glaner & J-.B.
% Uwineza
% *************************************************************************


% get variables from Adjust
A = Adjust.A;           % Design matrix; [number of obs x number of parameters]
P_l = Adjust.P;         % weight matrix of observations; square matrix
omc = Adjust.omc;       % vector of observations, observed minus computed [| vector]
P_x = Adjust.P_pred;	% weight matrix of (co)variance-matrix of predicted parameters

%--- Standard KF notation 
dx_k_minus      = Adjust.dx_k_minus;
P_k_minus       = Adjust.P_k_minus; 
% H_k             = Adjust.A;
R_k             = Adjust.Q; 
Phi_k           = Adjust.Transition;
Qd_k            = Adjust.Noise * 3600;   % 3600s = 1h, since the noise is given in **/sqrt(h) units 
% Qd_k            = zeros(size(Adjust.Noise));

%% Check for NaNs
% in observed-minus-computed, indication for missing observations
idx_nan = isnan(omc);
% set rows of missing observations to zero and remove them from adjustment
omc(idx_nan) = 0;
A(idx_nan, :) = 0;


%% --- Standard Iterated Kalman Filter (Ref: Aided Nav, Table 5.7)
% Off-line calculations
H_k = A;
K_k         = (P_k_minus*H_k') / (R_k + H_k*P_k_minus*H_k');  % using / instead of inv
P_k_plus    = (eye(size(P_k_minus)) - K_k*H_k) * P_k_minus;
P_k1_minus  = Phi_k*P_k_plus*Phi_k' + Qd_k;
P_k1_minus  = 1/2 * (P_k1_minus + P_k1_minus');  % to preserve symmetry

% Measurement update 
z_k = omc; 
dx_k_plus   = dx_k_minus + K_k*(z_k - (H_k*dx_k_minus)); 

% Time propagatiKkon 
dx_k1_minus = Phi_k*dx_k_plus; 

dx.P_k1_minus    = P_k1_minus; 
dx.dx_k1_minus   = dx_k1_minus;


%DEBUG
% remove the columns of those parameters which do not contribute to the 
% adjustment which should increase the numerical stability

zero_columns = all(A == 0, 1);                  % check which columns are zero 
zero_columns(Adjust.NO_PARAM+1:end) = 0;          % do not remove ambiguities or estimated ionosphere
A(:,zero_columns) = [];
P_x(:,zero_columns) = [];
P_x(zero_columns,:) = [];
x_pred(zero_columns) = [];

%% save results
idx = ~zero_columns;    	% removed columns have to be considered

dx.x = dx_k1_minus;
% Residuals/Verbesserungen
dx.v         = z_k - (H_k*dx_k_plus);          	  

% Covariance matrix of parameters
dx.Qxx = P_k_plus; 
% DEBUG: Use the standard results
% dx.x = dx_k1_minus; 
dx.res_var = R_k+H_k*P_k_plus*H_k'; 
end

