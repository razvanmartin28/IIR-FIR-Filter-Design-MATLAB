% CAUTA_ORDIN_MIN_CHEBY1 Gaseste ordinul minim pentru Chebyshev Tip 1
% 
%	Fisier: cauta_ordin_min_Cheby1.m
%
%	Functia: cauta_ordin_min_Cheby1
%
%	Apelare: [M_final, b_final, a_final] = cauta_ordin_min_Cheby1(M_start, Rp_dB, freq_p, idx_s, Delta_s, N_points) 
%
%   Functia gaseste ordinul minim si coeficientii unui filtru 
%   Chebyshev Tip 1 (oscilatii in banda de trecere) plecand de la 
%   un ordin de inceput. Verifica conditia in banda de stopare.
%
% Inputs:
%   M_start   - Ordinul de la care incepe cautarea
%   Rp_dB     - Ripple in banda de trecere (dB)
%   freq_p    - Frecventa limita banda trecere normalizata (0..1)
%   idx_s     - Indexul minim pentru banda de stopare in vectorul de frecvente
%   Delta_s   - Toleranta banda de stopare (liniar)
%   N_points  - Numarul de puncte pentru freqz
%
% Outputs:
%   M_final   - Ordinul minim gasit
%   b_final   - Coeficientii numaratorului
%   a_final   - Coeficientii numitorului
%
% Author: MARTIN RAZVAN-STEFAN
%	Created: January 6, 2026 

function [M_final, b_final, a_final] = cauta_ordin_min_Cheby1(M_start, Rp_dB, freq_p, idx_s, Delta_s, N_points)
    M_curr = M_start;
    found = false;
    M_max_limit = 50; % Protectie (filtrele IIR de ordin mare sunt instabile)
    
    while ~found
        % Proiectare Cheby1 (impune Rp si wp)
        [b_curr, a_curr] = cheby1(M_curr, Rp_dB, freq_p);
        
        % Verificare
        [H_test, ~] = freqz(b_curr, a_curr, N_points);
        H_abs = abs(H_test);
        
        % Conditii: Cheby1 garanteaza banda de trecere prin Rp.
        % Trebuie verificata doar banda de stopare.
        stop_ok = all(H_abs(idx_s:end) <= Delta_s);
        
        if stop_ok
            found = true;
            M_final = M_curr;
            b_final = b_curr;
            a_final = a_curr;
        else
            M_curr = M_curr + 1;
            if M_curr > M_max_limit
                warning('S-a atins limita maxima (%d) pentru Cheby1.', M_max_limit);
                M_final = M_curr;
                b_final = b_curr; a_final = a_curr;
                break; 
            end
        end
    end
    fprintf('Ordin minim Cheby1 gasit: %d\n', M_final);
end