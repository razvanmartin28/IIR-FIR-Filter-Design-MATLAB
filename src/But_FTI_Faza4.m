%
%	File BUT_FTI_FAZA4.M
%
%	Function: BUT_FTI_FAZA4
%
%	Synopsis: [B,A] = But_FTI_Faza4(w_p,w_s,Delta_p,Delta_s,Ts) ; 
%
%	Designs a low-pass discrete IIR filter of Butterworth class, 
%	by solving a design problem with specified tolerances. 
%	The design is realized by means of Tustin's (bilinear) 
%	discretizing transform Method. 
%   
%   MODIFICATION (Phase 4): The filter is designed to have a 
%   non-unitary DC gain: H(0) = 1 + Delta_p.
%
%	Inputs: w_p     = the relative passband upper limit 
%	                  (a number between 0 and 1); 
%	        w_s     = the relative stopband lower limit 
%	                  (a number between 0 and 1, at least equal to wp); 
%	        Delta_p = the tolerance in the passband (a number between 0 and 1); 
%	        Delta_s = the tolerance in the stopband (a number between 0 and 1); 
%	        Ts      = the sampling period required in the Tustin's method
%	                  (by default, Ts=2). 
%
%	Missing, empty or inconsistent inputs return empty or wrong output. 
%
%	Uses:	 WAR_ERR 
%
%	Authors: Bogdan DUMITRESCU & Dan STEFANOIU (Original)
%            Modified by MARTIN RAZVAN-STEFAN for Phase 4
%	Created: March 15, 2010 
%	Revised: January 11, 2026 (Phase 4 Modification)
%
function [B,A] = But_FTI_Faza4(w_p,w_s,Delta_p,Delta_s,Ts)
%
% BEGIN
%
% Constants & Messages 
% ~~~~~~~~~~~~~~~~~~~~
	FN = '<BUT_FTI_FAZA4>: ' ;
	E1 = [FN 'Missing, empty or inconsistent input data => empty outputs. Exit.'] ; 
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
	B = [] ; 
	A = [] ; 
	if (nargin < 3)
	   error(E1) ; % Changed to error() for standard Matlab use
	   return ; 
	end ; 
	w_p = abs(w_p(1)) ; 
	if (w_p < eps) || (w_p >= (1-eps))
	   error(E1) ;
	   return ; 
	end ; 
	w_s = abs(w_s(1)) ; 
	if (w_s < eps) || (w_s >= (1-eps))
	   error(E1) ;
	   return ; 
	end ; 
	Delta_p = abs(Delta_p(1)) ; 
	if (Delta_p < eps) || (Delta_p >= (1-eps))
	   error(E1) ;
	   return ; 
	end ; 
	if (nargin < 4)
	   Delta_s = Delta_p ;
	end ; 
	Delta_s = abs(Delta_s(1)) ; 
	if (Delta_s < eps) || (Delta_s >= (1-eps))
	   error(E1) ;
	   return ; 
	end ; 
	if (w_p > w_s)
	   temp = w_p ; 
	   w_p = w_s ; 
	   w_s = temp ; 
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
	wp_a = 2*tan(w_p*pi/2)/Ts ; 	% Compute Omega_p.
	ws_a = 2*tan(w_s*pi/2)/Ts ; 	% Compute Omega_s.
	
    % Phase 4 Modification: New formula for Order M calculation.
    % Based on conditions: |H(j*wp)| >= 1 and |H(j*ws)| <= Delta_s
    % with |H(0)| = 1 + Delta_p.
    
    num_term = ( (1 + Delta_p + Delta_s) * (1 + Delta_p - Delta_s) ) / ( Delta_p * (Delta_s^2) * (2 + Delta_p) );
    M = ceil( log(num_term) / (2 * log(ws_a/wp_a)) ); % Minimum order of filter.
    
    % Phase 4 Modification: New formula for Omega_c calculation.
    term_c = Delta_p * (2 + Delta_p);
    Omega_c = wp_a / ( term_c^(1/(2*M)) ); 
					
%
% Step #2: Build the transfer function of analog Butterworth filter
%
	Poles = Omega_c*exp(1j*(M+(1:2:(2*M)))*pi/M/2) ; % The stable poles of filter. 
    % We will ensure G(1) = 1 + Delta_p implicitly by carrying the gain.
    
%
% Step #3: Build the transfer function of digital Butterworth filter
%
	Poles = Poles*Ts ; 			% Normalize the poles for bilinear transform. 
	Denom_roots = 2 - Poles ; 
    % We must scale the numerator by (1 + Delta_p) because the analog prototype is scaled.
    
    Gain_Factor = (1 + Delta_p); 
    
	B = real(Gain_Factor * prod(-Poles./Denom_roots) * poly(-ones(1,M))) ; 	% Numerator of transfer function. 
	A = real(poly((2+Poles)./Denom_roots)) ; 			% Denominator of transfer function.
%
% END