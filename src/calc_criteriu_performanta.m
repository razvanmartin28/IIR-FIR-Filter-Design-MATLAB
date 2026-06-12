% CALC_CRITERIU_PERFORMANTA calculeaza un scor de performanta pentru un filtru
%
%   Fisier: calc_criteriu_performanta.m
%   Functia: calc_criteriu_performanta
%
%   Criteriul este un COST (valoare mai mica = filtru mai bun).
%   Se normalizeaza erorile pentru a avea aceeasi pondere numerica cu Ordinul,
%   apoi se aplica procentele dorite.
%
%   Formula: C = 0.6*M + 0.3*(E_mag_scalata) + 0.1*(E_phs_scalata)
%   Ponderi tinta: 60% M, 30% Eroare Mag, 10% Eroare Faza
%
% Inputs:
%   M           - Ordinul filtrului
%   H           - Raspunsul in frecventa (vector complex)
%   W           - Vectorul de frecvente normalizate (0..pi)
%   omega_p     - Pulsatia limita banda de trecere
%   omega_s     - Pulsatia limita banda de stopare
%   nume_filtru - String cu numele filtrului (ex: 'Butterworth')
%
% Outputs:
%   C       - Valoarea criteriului (Costul)
%
% Author: MARTIN RAZVAN-STEFAN
% Created: January 2026
function C = calc_criteriu_performanta(M, H, W, omega_p, omega_s, nume_filtru)
    
    % COMPONENTA MAGNITUDINE
    H_abs = abs(H);
    H_ideal = zeros(size(W));
    
    H_ideal(W <= omega_p) = 1;
    % Tranzitia nu se penalizeaza (eroare 0)
    idx_trans = (W > omega_p) & (W < omega_s);
    H_ideal(idx_trans) = H_abs(idx_trans);
    
    err_mag_vect = H_abs - H_ideal;
    rms_mag = sqrt(mean(err_mag_vect.^2));
    
    % COMPONENTA FAZA
    idx_pass = W <= omega_p;
    w_pass = W(idx_pass);
    phi_pass = unwrap(angle(H(idx_pass)));
    
    if length(w_pass) > 1
        p = polyfit(w_pass, phi_pass, 1); 
        phi_ideal_linear = polyval(p, w_pass);
        err_phs_vect = phi_pass - phi_ideal_linear;
        rms_phs = sqrt(mean(err_phs_vect.^2));
    else
        rms_phs = 0;
    end
    
    % SCALARE SI CALCUL COST FINAL
    Scale_Mag = 1140; 
    Scale_Phs = 24;   
    
    Termen_M   = M;
    Termen_Mag = rms_mag * Scale_Mag;
    Termen_Phs = rms_phs * Scale_Phs;
    
    % Aplicarea explicita a procentelor (60% / 30% / 10%)
    C = 0.6 * Termen_M + 0.3 * Termen_Mag + 0.1 * Termen_Phs;

    % Afisare pentru verificare cu numele filtrului
    fprintf('    %s Ordin M = %d, Total Cost: %.4f = M: %.2f (60%%) + Mag: %.2f (30%%) + Faza: %.2f (10%%) \n', ...
            nume_filtru, M, C, 0.6*Termen_M, 0.3*Termen_Mag, 0.1*Termen_Phs);
end