% CAUTA_ORDIN_MIN_CHEBY2 Gaseste ordinul minim pentru Chebyshev Tip 2
% 
%	Fisier: cauta_ordin_min_Cheby2.m
%
%	Functia: cauta_ordin_min_Cheby2
%
%	Apelare: [M_final, b_final, a_final] = cauta_ordin_min_Cheby2(M_start, Rs_dB, freq_s, idx_p, Delta_p, N_points) 
%
%   Functia gaseste ordinul minim si coeficientii unui filtru 
%   Chebyshev Tip 2 (oscilatii in banda de stopare, monoton in trecere)
%   Verifica conditia in banda de trecere.
%
% Inputs:
%   M_start   - Ordinul de la care incepe cautarea
%   Rs_dB     - Atenuare minima in banda de stopare (dB)
%   freq_s    - Frecventa limita banda stopare normalizata (0..1)
%   idx_p     - Indexul maxim pentru banda de trecere in vectorul de frecvente
%   Delta_p   - Toleranta banda de trecere (liniar)
%   N_points  - Numarul de puncte pentru freqz
%
% Outputs:
%   M_final   - Ordinul minim gasit
%   b_final   - Coeficientii numaratorului
%   a_final   - Coeficientii numitorului
%
% Author: MARTIN RAZVAN-STEFAN
%	Created: January 6, 2026 

function [M_final, b_final, a_final] = cauta_ordin_min_Cheby2(M_start, Rs_dB, freq_s, idx_p, Delta_p, N_points)
    M_curr = M_start;
    found = false;
    M_max_limit = 50; 
    
    while ~found
        % Proiectare Cheby2 (impune Rs si ws)
        [b_curr, a_curr] = cheby2(M_curr, Rs_dB, freq_s);
        
        % Verificare
        [H_test, ~] = freqz(b_curr, a_curr, N_points);
        H_abs = abs(H_test);
        
        % Conditii: Cheby2 garanteaza banda de stopare prin Rs.
        % Trebuie verificata banda de trecere [1-Dp, 1+Dp].
        pass_ok = all(H_abs(1:idx_p) >= (1 - Delta_p)) && ...
                  all(H_abs(1:idx_p) <= (1 + Delta_p));
        
        if pass_ok
            found = true;
            M_final = M_curr;
            b_final = b_curr;
            a_final = a_curr;
        else
            M_curr = M_curr + 1;
            if M_curr > M_max_limit
                warning('S-a atins limita maxima (%d) pentru Cheby2.', M_max_limit);
                M_final = M_curr;
                b_final = b_curr; a_final = a_curr;
                break; 
            end
        end
    end
    fprintf('Ordin minim Cheby2 gasit: %d\n', M_final);
end