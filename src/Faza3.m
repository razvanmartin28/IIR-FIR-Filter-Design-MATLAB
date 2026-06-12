
% FAZA 3: Concurs de Proiectare - Comparatie Globala (Ordin Minim)

clear; close all; clc;

% Date de Proiectare ---
ns = 9;
ng = 2;

% Calcul Specificatii
Delta_p = 5 / 100; % 0.05
Atten_dB = 30; 
Delta_s = 10^(-Atten_dB/20);  % ~0.0316
Ts = 2.5784; 
omega_p = PS_PRJ_3_Faza_3(ng, ns);
omega_s = omega_p + pi/33;

% Parametri auxiliari dB
Rp_dB = -20*log10(1 - Delta_p);
Rs_dB = -20*log10(Delta_s);

% Frecvente normalizate
Wp_n = omega_p / pi;
Ws_n = omega_s / pi;

% Grila verificare
N_pts = 5000;
[~, W] = freqz(1, 1, N_pts);
idx_p = find(W <= omega_p, 1, 'last');
idx_s = find(W >= omega_s, 1, 'first');

% AFISARE SPECIFICATII
fprintf('Specificatii:\n');
fprintf('wp = %.4f rad\n', omega_p);
fprintf('ws = %.4f rad\n', omega_s);
fprintf('Delta_p = %.4f\n', Delta_p);
fprintf('Delta_s = %.4f\n', Delta_s);
fprintf('Ts = %.4f s\n', Ts);


% Proiectarea celor 4 Filtre (Ordin Minim)

% Butterworth
[B_butt, A_butt] = But_FTI(Wp_n, Ws_n, Delta_p, Delta_s, Ts);

% Calculul parametrilor M si omega_c
wp_analog = 2*tan(omega_p/2)/Ts;
ws_analog = 2*tan(omega_s/2)/Ts;
Mp = 1 - Delta_p;
term_p = (1 - Mp^2)/Mp^2;
term_s = (1 - Delta_s^2)/Delta_s^2;

% Calcul Ordin M
M_butt = ceil(log(term_s/term_p) / (2*log(ws_analog/wp_analog)));

% Eliptic
[M_ellip, B_el, A_el] = cauta_ordin_min_Eliptic(1, Rp_dB, Rs_dB, Wp_n, idx_s, Delta_s, N_pts);

% Cebîsev Tip 1
[M_c1, B_c1, A_c1] = cauta_ordin_min_Cheby1(1, Rp_dB, Wp_n, idx_s, Delta_s, N_pts);

% Cebîsev Tip 2
[M_c2, B_c2, A_c2] = cauta_ordin_min_Cheby2(1, Rs_dB, Ws_n, idx_p, Delta_p, N_pts);


% Calcul Raspunsuri si Criterii Performanta

[H_butt_v, ~] = freqz(B_butt, A_butt, N_pts);
[H_c1_v, ~]   = freqz(B_c1, A_c1, N_pts);
[H_c2_v, ~]   = freqz(B_c2, A_c2, N_pts);
[H_el_v, ~]   = freqz(B_el, A_el, N_pts);

fprintf('\nRezultate Performanta C = 0.6*M + 0.3*ErrMag + 0.1*ErrFaza:\n');
C_butt = calc_criteriu_performanta(M_butt, H_butt_v, W, omega_p, omega_s, "Butterworth:");
C_c1   = calc_criteriu_performanta(M_c1, H_c1_v, W, omega_p, omega_s,     "Cheby1:     ");
C_c2   = calc_criteriu_performanta(M_c2, H_c2_v, W, omega_p, omega_s,     "Cheby2:     ");
C_el   = calc_criteriu_performanta(M_ellip, H_el_v, W, omega_p, omega_s,  "Cauer:      ");


%%
% Structura pentru sortare
Filtres(1) = struct('H', H_butt_v, 'M', M_butt, 'Name', "Butterworth", 'C', C_butt);
Filtres(2) = struct('H', H_c1_v, 'M', M_c1, 'Name', "Cheby1", 'C', C_c1);
Filtres(3) = struct('H', H_c2_v, 'M', M_c2, 'Name', "Cheby2", 'C', C_c2);
Filtres(4) = struct('H', H_el_v, 'M', M_ellip, 'Name', "Cauer", 'C', C_el);

% Sortare Crescatoare dupa C (Cost mic = Performanta mare)
[~, sort_idx] = sort([Filtres.C], 'ascend');
SortedFilters = Filtres(sort_idx);

% Variabile finale pentru salvare
H1 = SortedFilters(1).H; M1 = SortedFilters(1).M; H1_name = SortedFilters(1).Name;
H2 = SortedFilters(2).H; M2 = SortedFilters(2).M; H2_name = SortedFilters(2).Name;
H3 = SortedFilters(3).H; M3 = SortedFilters(3).M; H3_name = SortedFilters(3).Name;
H4 = SortedFilters(4).H; M4 = SortedFilters(4).M; H4_name = SortedFilters(4).Name;

save('Rezultate_Concurs_Proiectare.mat', 'omega_p', ...
     'H1', 'H2', 'H3', 'H4', ...
     'M1', 'M2', 'M3', 'M4', ...
     'H1_name', 'H2_name', 'H3_name', 'H4_name');



% Generare Grafice (Matrice 2x4)
figure('Name', 'Faza 3 - Clasament Filtre', 'Color', 'w', 'Position', [20 50 1500 800]);

% Limite plotare linii
lim_pass_high = 20*log10(1 + Delta_p);
lim_pass_low  = 20*log10(1 - Delta_p);
lim_stop      = 20*log10(Delta_s);

for k = 1:4
    curr_H = SortedFilters(k).H;
    curr_H_dB = 20*log10(abs(curr_H) + eps);
    curr_Phi = unwrap(angle(curr_H));
    curr_M = SortedFilters(k).M;
    curr_Name = SortedFilters(k).Name;
    curr_C = SortedFilters(k).C;
    
    % SPECTRE 
    subplot(2, 4, k);
    plot(W/pi, curr_H_dB, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
    
    % Linii specificatii
    h_wp = xline(omega_p/pi, '--r', 'DisplayName', '\omega_p');
    h_ws = xline(omega_s/pi, '--m', 'DisplayName', '\omega_s');
    yline(lim_pass_high, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'bottom', 'FontSize', 6);
    yline(lim_pass_low, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'top', 'FontSize', 6);
    yline(lim_stop, 'r', 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'top', 'FontSize', 6);
    
    info_str = sprintf('M = %d\nCost=%.2f', curr_M, curr_C);
    text(0.05, -50, info_str, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');
    
    title(sprintf('Locul %d: %s', k, curr_Name));
    xlim([0 1]); ylim([-60 5]);
    if k==1, ylabel('Amplitudine (dB)'); end
    legend([h_wp, h_ws], 'Location', 'best', 'FontSize', 8);
    
    % FAZE
    subplot(2, 4, 4 + k);
    plot(W/pi, curr_Phi, 'k', 'LineWidth', 1.2); grid on;
    title(['Faza ' char(curr_Name)]);
    xlim([0 1]);
    xlabel('Frecventa (\times\pi)');
    if k==1, ylabel('Faza (rad)'); end
end

print('-dpng', '-r300', '../Figuri/faza3_clasament.png');