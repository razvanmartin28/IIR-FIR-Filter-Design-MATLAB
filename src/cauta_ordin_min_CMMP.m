% CAUTA_ORDIN_MIN_CMMP Gaseste ordinul minim pentru metoda CMMP (firls)
%	Fisier: cauta_ordin_min_CMMP.m
%
%	Functia: cauta_ordin_min_CMMP
%
%	Apelare: [M_final, b_final] = cauta_ordin_min_CMMP(M_start, F_ls, A_ls, idx_p, idx_s, Delta_p, Delta_s, N_points) 
%   
%   Functia gaseste ordinul minim si coeficientii unui filtru 
%   cu ajutorul metodei celor mai mici patrate plecand de la
%   un ordin de inceput si avand datele de proiectare disponibile 
%
% Inputs:
%   M_start   - Ordinul de la care incepe cautarea
%   F_ls      - Vectorul de frecvente normalizate [0 wp ws 1]
%   A_ls      - Vectorul de amplitudini dorite [1 1 0 0]
%   idx_p     - Indexul maxim pentru banda de trecere in vectorul de frecvente
%   idx_s     - Indexul minim pentru banda de stopare in vectorul de frecvente
%   Delta_p   - Toleranta banda de trecere
%   Delta_s   - Toleranta banda de stopare
%   N_points  - Numarul de puncte pentru freqz
%
% Outputs:
%   M_final   - Ordinul minim gasit
%   b_final   - Coeficientii filtrului rezultat
%
% Author: MARTIN RAZVAN-STEFAN
%	Created: January 6, 2026 

function [M_final, b_final] = cauta_ordin_min_CMMP(M_start, F_ls, A_ls, idx_p, idx_s, Delta_p, Delta_s, N_points)
    M_curr = M_start;
    found = false;
    M_max_limit = 500; % Protectie la bucla infinita
    while ~found
        % Proiectare (apel specific firls)
        b_curr = firls(M_curr, F_ls, A_ls); 
        
        % Verificare
        [H_test, ~] = freqz(b_curr, 1, N_points);
        H_abs = abs(H_test);
        
        % Conditii
        pass_ok = all(H_abs(1:idx_p) >= (1 - Delta_p)) && ...
                  all(H_abs(1:idx_p) <= (1 + Delta_p));
        stop_ok = all(H_abs(idx_s:end) <= Delta_s);
        
        if pass_ok && stop_ok
            found = true;
            M_final = M_curr;
            b_final = b_curr;
        else
            M_curr = M_curr + 1;
            if M_curr > M_max_limit
                warning('S-a atins limita maxima de iteratii (%d) pentru FIRLS.', M_max_limit);
                M_final = M_curr;
                b_final = b_curr; % Returnam ultimul calculat
                break; 
            end
        end
    end
    fprintf('Ordin minim FIRLS gasit: %d\n', M_final);
end