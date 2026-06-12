% Autor: Martin Razvan-Stefan
%Subpunctul a
clear; clc; close all;
% Datele specifice
ns = 9;
ng = 2;
% [omega_p,omega_s,Delta_p,Ts] = PS_PRJ_3_Faza_1a(ng,ns);
% Salvam un set de date obtinut cu ajutorul functiei
omega_s = 1.6122;
omega_p = 1.3721;
Delta_p = 0.0776;
Delta_s = Delta_p;
Ts = 2.5784;

% Afisare date initiale
fprintf('Date de proiectare:\n');
fprintf('wp = %.4f rad (%.2f pi)\n', omega_p, omega_p/pi);
fprintf('ws = %.4f rad (%.2f pi)\n', omega_s, omega_s/pi);
fprintf('Delta_p = %.4f\n', Delta_p);
fprintf('Ts = %.2f s\n', Ts);

% Proiectarea filtrului folosind But_FTI
[B, A] = But_FTI(omega_p/pi, omega_s/pi, Delta_p, Delta_s, Ts);

% Calculul parametrilor M si omega_c 
% Folosim formulele 6.21, 6.26 si 6.27 din suportul de laborator
wp_analog = 2*tan(omega_p/2)/Ts;
ws_analog = 2*tan(omega_s/2)/Ts;
Mp = 1 - Delta_p;
term_p = (1 - Mp^2)/Mp^2;
term_s = (1 - Delta_s^2)/Delta_s^2;

% Calcul Ordin M
M = ceil(log(term_s/term_p) / (2*log(ws_analog/wp_analog)));

% Calcul Omega_c (analogic)
Omega_c = wp_analog / (term_p^(1/(2*M)));

% Calcul omega_c (discret)
omega_c = 2*atan(Omega_c*Ts/2);

fprintf('\nRezultate obtinute:\n');
fprintf('Ordinul filtrului M = %d\n', M);
fprintf('Pulsatia de taiere wc = %.4f rad (%.4f pi)\n', omega_c, omega_c/pi);

% Trasarea caracteristicilor
N_points = 5000;
[H, W] = freqz(B, A, N_points);
H_dB = 20*log10(abs(H));
Phi = unwrap(angle(H));

% FIGURA 1: Spectrul Filtrului Butterworth
figure('Name', 'Faza 1a - Spectru Amplitudine', 'Color', 'w');

plot(W/pi, H_dB, 'b', 'LineWidth', 1.5); hold on;
grid on;

% Linii verticale (pulsatiile)
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_wc = xline(omega_c/pi, '-.k', 'LineWidth', 1, 'DisplayName', '\omega_c');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');

% Linii orizontale 
% Toleranta de trecere
lim_pass_low = 20*log10(1-Delta_p);
lim_pass_high = 20*log10(1+Delta_p);

% Linia de sus (1 + Delta_p)
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', ...
    'LabelVerticalAlignment', 'top');

% Linia de jos (1 - Delta_p)
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', ...
    'LabelVerticalAlignment', 'bottom');

% Toleranta de stop (Delta_s)
lim_stop = 20*log10(Delta_s);
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', ...
    'LabelVerticalAlignment', 'bottom');

% Adnotari Text (Ts si M)
info_txt = sprintf('M = %d\nTs = %.2f s', M, Ts);
text(0.05, -15, info_txt, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

title('Spectrul Filtrului Butterworth');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Amplitudine (dB)');

legend([h_wp, h_wc, h_ws], 'Location', 'SouthWest');

xlim([0 1]);
ylim([-50 5]); 

% Salvare Figura 1
print('-dpng', '-r300', '../Figuri/spectru_faza1a.png');

% FIGURA 2: Faza Filtrului Butterworth
figure('Name', 'Faza 1a - Faza', 'Color', 'w');
plot(W/pi, Phi, 'r', 'LineWidth', 1.5);
grid on;

title('Faza Filtrului Butterworth');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Faza (rad)');
xlim([0 1]);

% Salvare Figura 2
print('-dpng', '-r300', '../Figuri/faza_faza1a.png');

%% Subpunctul b - Transformarea Biliniara Modificata

% Proiectarea filtrului folosind But_FTI_v2
% Folosim aceeasi paramterii ca la a
[B2, A2] = But_FTI_v2(omega_p/pi, omega_s/pi, Delta_p, Delta_s, Ts);

% Calculul parametrilor M si omega_c pentru noua transformare
% Formulele se adapteaza: dispare factorul 2 din definitia frecventei analogice
wp_analog_v2 = tan(omega_p/2)/Ts;
ws_analog_v2 = tan(omega_s/2)/Ts;

Mp = 1 - Delta_p;
term_p = (1 - Mp^2)/Mp^2;
term_s = (1 - Delta_s^2)/Delta_s^2;

% Calcul Ordin M (ar trebui sa fie identic)
M2 = ceil(log(term_s/term_p) / (2*log(ws_analog_v2/wp_analog_v2)));

% Calcul Omega_c (analogic)
Omega_c2 = wp_analog_v2 / (term_p^(1/(2*M2)));

% Calcul omega_c (discret) - formula inversa
% Daca Omega = (1/Ts)*tan(w/2) => w = 2*atan(Omega*Ts)
omega_c2 = 2*atan(Omega_c2*Ts); 

fprintf('Rezultate obtinute:\n');
fprintf('Ordinul filtrului M = %d\n', M2);
fprintf('Pulsatia de taiere wc = %.4f rad (%.4f pi)\n', omega_c2, omega_c2/pi);

% Trasarea caracteristicilor
[H2, W2] = freqz(B2, A2, N_points);
H2_dB = 20*log10(abs(H2));
Phi2 = unwrap(angle(H2));

% FIGURA 3: Spectrul Filtrului Butterworth (V2)
figure('Name', 'Faza 1b - Spectru Amplitudine', 'Color', 'w');
plot(W2/pi, H2_dB, 'b', 'LineWidth', 1.5); hold on;
grid on;

% Linii verticale (pulsatiile)
h_wp2 = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_wc2 = xline(omega_c2/pi, '-.k', 'LineWidth', 1, 'DisplayName', '\omega_c');
h_ws2 = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');

% Linii orizontale 
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', ...
    'LabelVerticalAlignment', 'top');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', ...
    'LabelVerticalAlignment', 'bottom');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', ...
    'LabelVerticalAlignment', 'bottom');

% Adnotari Text
info_txt2 = sprintf('M = %d\nTs = %.2f s', M2, Ts);
text(0.05, -15, info_txt2, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

title('Spectrul Filtrului Butterworth (V2)');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Amplitudine (dB)');
legend([h_wp2, h_wc2, h_ws2], 'Location', 'SouthWest');
xlim([0 1]);
ylim([-50 5]); 

% Salvare Figura 3
print('-dpng', '-r300', '../Figuri/spectru_faza1b.png');

% FIGURA 4: Faza Filtrului Butterworth (V2)
figure('Name', 'Faza 1b - Faza', 'Color', 'w');
plot(W2/pi, Phi2, 'r', 'LineWidth', 1.5);
grid on;

title('Faza Filtrului Butterworth (V2)');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Faza (rad)');
xlim([0 1]);

% Salvare Figura 4
print('-dpng', '-r300', '../Figuri/faza_faza1b.png');

% 4. Calculul si afisarea erorilor
% Calculam diferenta fata de filtrul de la punctul a (H si Phi calculate anterior)
Err_Spectru = abs(H_dB - H2_dB); % Diferenta in dB
Err_Faza = abs(Phi - Phi2);      % Diferenta de faza

% Calcul norma (Euclidiana / Frobenius pentru vectori e la fel)
norm_err_spectru = norm(H - H2); % Norma diferentei vectorilor complexi H
norm_err_faza = norm(Err_Faza);

fprintf('Norma diferentei functiilor de transfer: %e\n', norm_err_spectru);
fprintf('Norma diferentei fazelor: %e\n', norm_err_faza);

% FIGURA 5: Grafice de eroare
figure('Name', 'Faza 1b - Erori', 'Color', 'w');

% Subplot 1: Eroare Spectre
subplot(2,1,1);
plot(W/pi, Err_Spectru, 'k', 'LineWidth', 1.5);
grid on;
title('Eroare Spectre |H_{dB} - H2_{dB}|');
xlabel('Frecventa normalizata');
ylabel('Eroare (dB)');
% Afisare norma pe grafic
text(0.1, max(Err_Spectru)*0.8, sprintf('Norma ||H-H2|| = %e', norm_err_spectru), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

% Subplot 2: Eroare Faze
subplot(2,1,2);
plot(W/pi, Err_Faza, 'm', 'LineWidth', 1.5);
grid on;
title('Eroare Faze |\Phi - \Phi2|');
xlabel('Frecventa normalizata');
ylabel('Eroare (rad)');
text(0.1, max(Err_Faza)*0.8, sprintf('Norma ||Phi-Phi2|| = %e', norm_err_faza), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k');

% Salvare Figura 5
print('-dpng', '-r300', '../Figuri/erori_faza1b.png');

%% Subpunctul c - Analiza dependentei de Ts 

% Definim seturile de factori
factors_low = [0.1 0.25 0.5 0.75];
factors_high = [1.25 1.75 2.25 3];
sets = {factors_low, factors_high};
fig_titles = {'Faza 1c - Ts mai mic decat referinta', 'Faza 1c - Ts mai mare decat referinta'};
file_names = {'fig_faza1c_low.png', 'fig_faza1c_high.png'};

for fig_idx = 1:2
    current_factors = sets{fig_idx};
    
    h_fig = figure('Name', fig_titles{fig_idx}, 'Color', 'w', ...
           'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]); 
    
    for col = 1:4
        % 1. Calcul noul Ts
        factor = current_factors(col);
        Ts_new = factor * Ts; 
        
        % 2. Proiectare
        [B_new, A_new] = But_FTI(omega_p/pi, omega_s/pi, Delta_p, Delta_s, Ts_new);
        
        % 3. Raspuns in frecventa
        [H_new, ~] = freqz(B_new, A_new, N_points);
        H_dB_new = 20*log10(abs(H_new) + eps); % Adaugam eps pentru a evita log(0)
        Phi_new = unwrap(angle(H_new));
        
        % 4. Calcul erori (CORECTIE)
        
        % Eroarea liniara (pentru norma corecta)
        Err_Lin = abs(H) - abs(H_new);
        norm_err_spectru = norm(Err_Lin); % Norma pe diferenta liniara
        
        % Eroarea de faza (ignoram punctele unde magnitudinea e aproape 0)
        % Diferenta pe cercul unitate elimina salturile de 2*pi
        Err_Faza_Raw = abs(angle(exp(1j*(Phi - Phi_new))));
        % Punem 0 eroare acolo unde semnalul e prea mic (stopband adanc)
        mask_valid = abs(H) > 1e-10; 
        norm_err_faza = norm(Err_Faza_Raw(mask_valid));
        
        % --- PLOTARE ---
        
        % Row 1: Spectrul (dB)
        subplot(4, 4, col); 
        plot(W/pi, H_dB_new, 'b'); grid on;
        title(sprintf('Ts = %.2f s', Ts_new));
        if col==1, ylabel('Ampl. (dB)'); end
        ylim([-100 10]); xlim([0 1]);
        
        % Row 2: Eroarea (Afisam diferenta liniara care e stabila)
        subplot(4, 4, 4 + col);
        plot(W/pi, abs(Err_Lin), 'k'); grid on;
        title(['Err Mag (Norm=' sprintf('%.1e', norm_err_spectru) ')'], 'FontSize', 8);
        if col==1, ylabel('Err (liniar)'); end
        xlim([0 1]);
        
        % Row 3: Faza
        subplot(4, 4, 8 + col);
        plot(W/pi, Phi_new, 'r'); grid on;
        if col==1, ylabel('Faza (rad)'); end
        xlim([0 1]);
        
        % Row 4: Eroarea de Faza
        subplot(4, 4, 12 + col);
        plot(W/pi, Err_Faza_Raw, 'm'); grid on;
        title(['Err Faza (Norm=' sprintf('%.1e', norm_err_faza) ')'], 'FontSize', 8);
        if col==1, ylabel('Err (rad)'); end
        xlim([0 1]);
    end
    
    % Salvare figura
    print(h_fig, '-dpng', '-r300', file_names{fig_idx});
end

%% Subpunctul d: Studiu comparativ - Variatia tolerantelor

% Definirea combinatiilor
factors = [0.5, 1, 1.5, 2];

% Generam toate cele 16 combinatii
combs = [];
for f_p = factors
    for f_s = factors
        combs = [combs; f_p*Delta_p, f_s*Delta_p];
    end
end
% combs este o matrice 16x2. Col 1: Dp, Col 2: Ds.

% Generarea graficelor (2 figuri a cate 8 cazuri)
% Fig 1: Cazurile 1-8, Fig 2: Cazurile 9-16
titles_d = {'Spectrul si Faza filtrului cu toleranta variata (comb 1-8)', 'Spectrul si Faza filtrului cu toleranta variata (comb 9-16)'};
file_names_d = {'../Figuri/fig_faza1d_set1.png', '../Figuri/fig_faza1d_set2.png'};

case_idx = 1;

for fig_k = 1:2
    h_fig = figure('Name', titles_d{fig_k}, 'Color', 'w', ...
           'Units', 'normalized', 'Position', [0.05 0.05 0.7 0.85]);
    
    % Iteram prin 8 cazuri per figura
    for k = 1:8
        if case_idx > 16, break; end
        
        % Extragere tolerante curente
        curr_Dp = combs(case_idx, 1);
        curr_Ds = combs(case_idx, 2);
        
        % Proiectare filtru cu pulsatiile de la subpunctul a
        [B_d, A_d] = But_FTI(omega_p/pi, omega_s/pi, curr_Dp, curr_Ds, Ts);

        % Recalculam M pentru a-l afisa
        wp_a = 2*tan(omega_p/2)/Ts;
        ws_a = 2*tan(omega_s/2)/Ts;
        Mp_curr = 1 - curr_Dp;
        term_p = (1 - Mp_curr^2)/Mp_curr^2;
        term_s = (1 - curr_Ds^2)/curr_Ds^2;
        M_curr = ceil(log(term_s/term_p) / (2*log(ws_a/wp_a)));
        
        % Raspuns in frecventa
        [H_d, W_d] = freqz(B_d, A_d, N_points);
        H_dB_d = 20*log10(abs(H_d)+eps);
        Phi_d = unwrap(angle(H_d));
        
        % --- PLOTARE ---
        % Primele 4 cazuri pe randurile 1-2, urmatoarele 4 pe randurile 3-4.
        if k <= 4
            row_offset_spec = 0; % Randul 1
            row_offset_phas = 4; % Randul 2
            col_pos = k;
        else
            row_offset_spec = 8;  % Randul 3
            row_offset_phas = 12; % Randul 4
            col_pos = k - 4;
        end
        
        % Plot Spectru
        subplot(4, 4, row_offset_spec + col_pos);
        plot(W_d/pi, H_dB_d, 'b'); grid on;
         % Afisam Delta_p si Delta_s pe graficul de amplitudine
        title(sprintf('\\Delta_p=%.3f, \\Delta_s=%.3f', curr_Dp, curr_Ds), 'FontSize', 8);
        ylim([-60 5]); xlim([0 1]);
        if col_pos == 1, ylabel('Ampl (dB)'); end
        
        % Plot Faza
        subplot(4, 4, row_offset_phas + col_pos);
        plot(W_d/pi, Phi_d, 'r'); grid on;
        % Afisam Ordinul M pe graficul de faza
        text(0.1, 0, sprintf('M = %d', M_curr), ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 8);
        xlim([0 1]);
        if col_pos == 1, ylabel('Faza (rad)'); end
        
        case_idx = case_idx + 1;
    end
    
    % Salvare
    print(h_fig, '-dpng', '-r300', file_names_d{fig_k});
end

%% Subpunctul e: Comparatie cu filtre FIR (fir1 si firls)

% Pregatirea datelor
freq_c = omega_c / pi; 
freq_p = omega_p / pi;
freq_s = omega_s / pi;
F_ls = [0 freq_p freq_s 1];
A_ls = [1 1 0 0];

% Limite pentru verificare
N_points = 5000;
[~, W_check] = freqz(1, 1, N_points); 
idx_p = find(W_check <= omega_p, 1, 'last');
idx_s = find(W_check >= omega_s, 1, 'first');

% Metoda Ferestrei - Apelare metoda
[M_fir1, b_fir1] = cauta_ordin_min_Fereastra(M-1, freq_c, idx_p, idx_s, Delta_p, Delta_s, N_points);

% Metoda Celor mai Mici Patrate - Apelare metoda
[M_firls, b_firls] = cauta_ordin_min_CMMP(M-1, F_ls, A_ls, idx_p, idx_s, Delta_p, Delta_s, N_points);


% Generare Grafice Finale

% Recalculare raspunsuri pentru plotare
[H_fir1, W] = freqz(b_fir1, 1, N_points);
H_dB_fir1 = 20*log10(abs(H_fir1)+eps);
Phi_fir1 = unwrap(angle(H_fir1));

[H_firls, ~] = freqz(b_firls, 1, N_points);
H_dB_firls = 20*log10(abs(H_firls)+eps);
Phi_firls = unwrap(angle(H_firls));

% SPECTRE (Stanga: FIR1, Dreapta: FIRLS)
figure('Name', 'Faza 1e - Spectre Comparate', 'Color', 'w', 'Position', [100 100 1200 500]);

% Spectru Fereastra
subplot(1, 2, 1);
p_sig = plot(W/pi, H_dB_fir1, 'b', 'LineWidth', 1.5, 'DisplayName', 'H_{FIR1}'); hold on; grid on;
% Linii Verticale
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_wc = xline(omega_c/pi, '-.k', 'LineWidth', 1, 'DisplayName', '\omega_c');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');

% Linii Orizontale
% Toleranta Trecere (Sus si Jos)
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
% Toleranta Stopare
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);

text(0.05, -50, sprintf('Ordin M = %d', M_fir1), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontSize', 10, 'FontWeight', 'bold');

title('Spectru Metoda Ferestrei');
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Amplitudine (dB)');
xlim([0 1]); ylim([-80 5]);
legend([h_wp, h_wc, h_ws], 'Location', 'southwest');

% Spectru CMMP
subplot(1, 2, 2);
p_sig2 = plot(W/pi, H_dB_firls, 'b', 'LineWidth', 1.5, 'DisplayName', 'H_{FIRLS}'); hold on; grid on;

% Linii Verticale
h_wp2 = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_wc2 = xline(omega_c/pi, '-.k', 'LineWidth', 1, 'DisplayName', '\omega_c');
h_ws2 = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');

% Linii Orizontale
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);

text(0.05, -50, sprintf('Ordin M = %d', M_firls), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontSize', 10, 'FontWeight', 'bold');

title('Spectru Metoda CMMP');
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Amplitudine (dB)');
xlim([0 1]); ylim([-80 5]);
legend([h_wp2, h_wc2, h_ws2], 'Location', 'southwest');

print('-dpng', '-r300', '../Figuri/spectre_sidebyside_faza1e.png');


% FAZE (Stanga: FIR1, Dreapta: FIRLS)
figure('Name', 'Faza 1e - Faze Comparate', 'Color', 'w', 'Position', [100 100 1200 500]);

% Faza Fereastra
subplot(1, 2, 1);
plot(W/pi, Phi_fir1, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda Ferestrei (M=' num2str(M_fir1) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)');
xlim([0 1]);

% Faza CMMP
subplot(1, 2, 2);
plot(W/pi, Phi_firls, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda CMMP (M=' num2str(M_firls) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)');
xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faze_sidebyside_faza1e.png');


% CALCUL ERORI (Magnitudine si Faza) 

% Erori Magnitudine (Liniar)
Err_Mag_FIR1 = abs(abs(H) - abs(H_fir1));
norm_mag_fir1 = norm(Err_Mag_FIR1);

Err_Mag_FIRLS = abs(abs(H) - abs(H_firls));
norm_mag_firls = norm(Err_Mag_FIRLS);

% Erori Faza (Radiani)
Err_Phs_FIR1 = abs(Phi - Phi_fir1);
norm_phs_fir1 = norm(Err_Phs_FIR1);

Err_Phs_FIRLS = abs(Phi - Phi_firls);
norm_phs_firls = norm(Err_Phs_FIRLS);

% Matrice 2x2 cu erorile
figure('Name', 'Faza 1e - Erori Absolute (Spectru si Faza)', 'Color', 'w', 'Position', [100 100 1000 700]);

% Metoda Ferestrei
% Eroare Magnitudine
subplot(2, 2, 1);
plot(W/pi, Err_Mag_FIR1, 'k', 'LineWidth', 1.2); grid on;
title('Eroare Spectru: Butterworth vs Fereastra');
ylabel('Eroare (liniar)'); xlim([0 1]);
text(0.05, max(Err_Mag_FIR1)*0.85, sprintf('Norma = %.2f', norm_mag_fir1), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

% Eroare Faza
subplot(2, 2, 3);
plot(W/pi, Err_Phs_FIR1, 'm', 'LineWidth', 1.2); grid on;
title('Eroare Faza: Butterworth vs Fereastra');
ylabel('Eroare (rad)'); xlim([0 1]);
text(0.05, max(Err_Phs_FIR1)*0.85, sprintf('Norma = %.2f', norm_phs_fir1), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

% Metoda CMMP
% Eroare Magnitudine
subplot(2, 2, 2);
plot(W/pi, Err_Mag_FIRLS, 'k', 'LineWidth', 1.2); grid on;
title('Eroare Spectru: Butterworth vs CMMP');
ylabel('Eroare (liniar)'); xlabel('Frecventa normalizata (\times\pi)'); xlim([0 1]);
text(0.05, max(Err_Mag_FIRLS)*0.85, sprintf('Norma = %.2f', norm_mag_firls), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

% Eroare Faza
subplot(2, 2, 4);
plot(W/pi, Err_Phs_FIRLS, 'm', 'LineWidth', 1.2); grid on;
title('Eroare Faza: Butterworth vs CMMP');
ylabel('Eroare (rad)'); xlabel('Frecventa normalizata (\times\pi)'); xlim([0 1]);
text(0.05, max(Err_Phs_FIRLS)*0.85, sprintf('Norma = %.2f', norm_phs_firls), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

print('-dpng', '-r300', '../Figuri/erori_complete_faza1e.png');
