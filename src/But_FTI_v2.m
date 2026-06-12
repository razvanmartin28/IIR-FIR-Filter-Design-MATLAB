%
%	File BUT_FTI_v2.M
%
%	Function: BUT_FTI_v2
%
%	Synopsis: [B,A] = But_FTI_v2(w_p,w_s,Delta_p,Delta_s,Ts) ; 
%
%   MODIFIED VERSION: Uses the transform s = (1/Ts) * (z-1)/(z+1)
%
function [B,A] = But_FTI_v2(w_p,w_s,Delta_p,Delta_s,Ts)
%
% BEGIN
%
% Constants & Messages 
% ~~~~~~~~~~~~~~~~~~~~
	FN = '<BUT_FTI_v2>: ' ;
	E1 = [FN 'Missing, empty or inconsistent input data => empty outputs. Exit.'] ; 
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
	B = [] ; 
	A = [] ; 
	if (nargin < 3)
	   war_err(E1) ;
	   return ; ; 
	end ; 
	w_p = abs(w_p(1)) ; 
	if (w_p < eps) || (w_p >= (1-eps))
	   war_err(E1) ;
	   return ; ; 
	end ; 
	w_s = abs(w_s(1)) ; 
	if (w_s < eps) || (w_s >= (1-eps))
	   war_err(E1) ;
	   return ; ; 
	end ; 
	Delta_p = abs(Delta_p(1)) ; 
	if (Delta_p < eps) || (Delta_p >= (1-eps))
	   war_err(E1) ;
	   return ; ; 
	end ; 
	if (nargin < 4)
	   Delta_s = Delta_p ;
	end ; 
	Delta_s = abs(Delta_s(1)) ; 
	if (Delta_s < eps) || (Delta_s >= (1-eps))
	   war_err(E1) ;
	   return ; ; 
	end ; 
	if (w_p > w_s)
	   FN = w_p ; 
	   w_p = w_s ; 
	   w_s = FN ; 
	end ; 
	if (nargin < 5)
	   Ts = 2 ; 
	end ; 
	Ts = abs(Ts(1)) ; 
	if (Ts < eps)
	   Ts = 2 ; 
	end ;
%
% Filter design
% ~~~~~~~~~~~~~
% Step #1: Find the parameters of analog Butterworth filter
%
    % MODIFICARE: s = (1/Ts)... deci Omega = (1/Ts)*tan(w/2)
    % A disparut factorul 2 de la numarator.
	w_p = tan(w_p*pi/2)/Ts ; 	% Compute Omega_p.
	w_s = tan(w_s*pi/2)/Ts ; 	% Compute Omega_s.
	FN = 1-Delta_p ; 		% Actually, this is M_p. 
	FN = FN*FN ; 
	FN = (1-FN)/FN ; 
	Delta_s = Delta_s*Delta_s ; 
	M = ceil(log((1-Delta_s)/Delta_s/FN)/log(w_s/w_p)/2) ;	% Minimum order of filter. 
	E1 = w_p/((FN)^(1/2/M)) ; 	% Actualy, this is Omega_c 
					% (the cut-off absolute pulsation).
%
% Step #2: Build the transfer function of analog Butterworth filter
%
	FN = E1*exp(j*(M+(1:2:(2*M)))*pi/M/2) ; % The stable poles of filter. 

%
% Step #3: Build the transfer function of digital Butterworth filter
%
    % MODIFICARE: Transformarea s = (1/Ts)*(z-1)/(z+1) implica z = (1+sTs)/(1-sTs)
    % FN contine polii s_k.
	FN = FN*Ts ; 			% Normalize the poles (s_k * Ts). 
    
    % In formula originala (z = (2+sTs)/(2-sTs)) era E1 = 2-FN.
    % Aici devine E1 = 1-FN.
	E1 = 1-FN ; 
    
    % La numarator si numitor, inlocuim 2 cu 1 conform noii transformari.
	B = real(prod(-FN./E1)*poly(-ones(1,M))) ; 	% Numerator of transfer function. 
	A = real(poly((1+FN)./E1)) ; 			% Denominator of transfer function.
%
% END
%
end