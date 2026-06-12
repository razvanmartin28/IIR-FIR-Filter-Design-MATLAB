% CAUTA_ORDIN_MIN_ELIPTIC Gaseste ordinul minim pentru filtrul Eliptic (Cauer)
% 
%	Fisier: cauta_ordin_min_Eliptic.m
%
%	Functia: cauta_ordin_min_Eliptic
%
%	Apelare: [M_final, b_final, a_final] = cauta_ordin_min_Eliptic(M_start, Rp_dB, Rs_dB, Wp_norm, idx_s, Delta_s, N_points)
%
%   Functia gaseste ordinul minim si coeficientii unui filtru 
%   Eliptic (Cauer) care satisface tolerantele impuse, plecand de la 
%   un ordin de inceput. Verifica conditia in banda de stopare.
%
% Inputs:
%   M_start   - Ordinul de la care incepe cautarea (de obicei 1)
%   Rp_dB     - Ripple in banda de trecere (dB)
%   Rs_dB     - Atenuare minima in banda de stopare (dB)
%   Wp_norm   - Frecventa limita banda trecere normalizata (0..1)
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
function [M_final, b_final, a_final] = cauta_ordin_min_Eliptic(M_start, Rp_dB, Rs_dB, Wp_norm, idx_s, Delta_s, N_points)
    M_curr = M_start;
    found = false;
    M_max_limit = 50; % Protectie
    
    while ~found
        % Proiectare Eliptic (Cauer)
        [b_curr, a_curr] = ellip(M_curr, Rp_dB, Rs_dB, Wp_norm);
        
        % Verificare
        [H_test, ~] = freqz(b_curr, a_curr, N_points);
        H_abs = abs(H_test);
        
        % Conditii: Eliptic garanteaza banda de trecere prin Rp si banda de stopare prin Rs.
        % Totusi, functia ellip nu accepta ws ca input, ci wp. 
        % Deci tranzitia este dictata de ordin. Trebuie sa verificam daca
        % atenuarea Rs este atinsa DEJA la frecventa ws (idx_s).
        
        stop_ok = all(H_abs(idx_s:end) <= Delta_s);
        
        if stop_ok
            found = true;
            M_final = M_curr;
            b_final = b_curr;
            a_final = a_curr;
        else
            M_curr = M_curr + 1;
            if M_curr > M_max_limit
                warning('S-a atins limita maxima (%d) pentru Eliptic.', M_max_limit);
                M_final = M_curr;
                b_final = b_curr; a_final = a_curr;
                break; 
            end
        end
    end
    fprintf('Ordin minim Eliptic gasit: %d\n', M_final);
end