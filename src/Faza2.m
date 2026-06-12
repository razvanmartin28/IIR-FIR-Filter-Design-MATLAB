clear; close all; clc;
%subpunctul a - filtre Cauer
% Specificatiile de la Faza 1a, dar Delta_s este DUBLAT.
omega_p = 1.3721;
omega_s = 1.6122;
Delta_p = 0.0776;
Delta_s = 2 * Delta_p; % Cerinta specifica Faza 2a
Ts = 2.5784;

% Afisare date initiale
fprintf('Date de proiectare:\n');
fprintf('wp = %.4f rad (%.2f pi)\n', omega_p, omega_p/pi);
fprintf('ws = %.4f rad (%.2f pi)\n', omega_s, omega_s/pi);
fprintf('Delta_p = %.4f\n', Delta_p);
fprintf('Ts = %.2f s\n', Ts);

% Recalculam Butterworth pentru noile specificatii
% Ne asiguram ca avem ordinul corect de comparatie.

% Calcul frecvente analogice predistorsionate
wp_a = 2*tan(omega_p/2)/Ts; 
ws_a = 2*tan(omega_s/2)/Ts;
Mp = 1 - Delta_p;
term_p = (1 - Mp^2)/Mp^2; 
term_s = (1 - Delta_s^2)/Delta_s^2;
M_but = ceil(log(term_s/term_p) / (2*log(ws_a/wp_a)));
fprintf('Ordin Butterworth necesar: %d\n', M_but);

% Proiectare efectiva Butterworth
[B_but, A_but] = But_FTI(omega_p/pi, omega_s/pi, Delta_p, Delta_s, Ts);


% Proiectare Filtru Eliptic (Cauer) - Ordin Minim
% Parametrii pentru ellip (in dB)
Rp_dB = -20*log10(1 - Delta_p); % Ripple in passband
Rs_dB = -20*log10(Delta_s);     % Attenuation in stopband
Wp_norm = omega_p / pi;         % Frecventa normalizata

% Pregatire verificare
N = 5000;
[~, W_check] = freqz(1, 1, N);
idx_s = find(W_check >= omega_s, 1, 'first');

% Filtru Eliptic (Cauer) - Apel Functie
[M_ellip, b_el, a_el] = cauta_ordin_min_Eliptic(1, Rp_dB, Rs_dB, Wp_norm, idx_s, Delta_s, N);


% Generare Grafice

% Calcul raspunsuri in frecventa
[H_but, W] = freqz(B_but, A_but, N);
H_dB_but = 20*log10(abs(H_but)+eps);
Phi_but = unwrap(angle(H_but));

[H_ellip, ~] = freqz(b_el, a_el, N);
H_dB_ellip = 20*log10(abs(H_ellip)+eps);
Phi_ellip = unwrap(angle(H_ellip));

% Limite grafice
lim_pass_low = 20*log10(1-Delta_p);
lim_pass_high = 20*log10(1+Delta_p);
lim_stop = 20*log10(Delta_s);

% SPECTRE COMPARATE (Eliptic vs Butterworth)
figure('Name', 'Faza 2a - Spectre', 'Color', 'w', 'Position', [100 100 1200 500]);

% Subplot Stanga: Eliptic
subplot(1, 2, 1);
plot(W/pi, H_dB_ellip, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
% Linii verticale
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
% Linii orizontale
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', 'LabelVerticalAlignment', 'bottom');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', 'LabelVerticalAlignment', 'top');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', 'LabelVerticalAlignment', 'bottom');
text(0.05, -50, sprintf('Ordin M = %d', M_ellip), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Filtru Eliptic (Cauer)');
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Amplitudine (dB)');
xlim([0 1]); ylim([-80 5]);
legend([h_wp, h_ws], 'Location', 'southwest');

% Subplot Dreapta: Butterworth
subplot(1, 2, 2);
plot(W/pi, H_dB_but, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
% Linii verticale
h_wp2 = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws2 = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
% Linii orizontale
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', 'LabelVerticalAlignment', 'bottom');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', 'LabelVerticalAlignment', 'top');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', 'LabelVerticalAlignment', 'bottom');
text(0.05, -50, sprintf('Ordin M = %d', M_but), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Filtru Butterworth');
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Amplitudine (dB)');
xlim([0 1]); ylim([-80 5]);
legend([h_wp2, h_ws2], 'Location', 'southwest');

print('-dpng', '-r300', '../Figuri/spectre_faza2a.png');


% FAZE COMPARATE (Eliptic vs Butterworth)
figure('Name', 'Faza 2a - Faze', 'Color', 'w', 'Position', [100 100 1200 500]);

% Subplot Stanga: Eliptic
subplot(1, 2, 1);
plot(W/pi, Phi_ellip, 'b', 'LineWidth', 1.2); grid on;
title(['Faza Filtru Eliptic (M=' num2str(M_ellip) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)');
xlim([0 1]);

% Subplot Dreapta: Butterworth
subplot(1, 2, 2);
plot(W/pi, Phi_but, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2); grid on;
title(['Faza Filtru Butterworth (M=' num2str(M_but) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)');
xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faze_faza2a.png');


% EROARE SPECTRU
Err_Mag = abs(abs(H_but) - abs(H_ellip));
norm_mag = norm(Err_Mag);

figure('Name', 'Faza 2a - Eroare Spectru', 'Color', 'w');
plot(W/pi, Err_Mag, 'k', 'LineWidth', 1.2); grid on;
title('Eroare spectru: Butterworth vs. Eliptic (Cauer)');
ylabel('Eroare (liniar)'); xlabel('Frecventa normalizata (\times\pi)');
xlim([0 1]);
text(0.1, max(Err_Mag)*0.85, sprintf('Norma = %.2f', norm_mag), ...
     'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

print('-dpng', '-r300', '../Figuri/erori_faza2a.png');

%% Subpunctul b: Solutii FIR si Comparatie Globala (Cauer, Butt, FIR1, FIRLS)

% Calculam omega_c bazat pe Butterworth
Om_p = (2/Ts) * tan(omega_p/2);
Om_s = (2/Ts) * tan(omega_s/2);
% Calcul omega_c analogic
eps_p = sqrt((1/(1-Delta_p)^2) - 1);
Om_c = Om_p / (eps_p^(1/M_but)); 
% Revenire in discret
omega_c_new = 2 * atan((Om_c * Ts)/2);
freq_c = omega_c_new / pi;

% Vectori pentru firls
freq_p = omega_p / pi;
freq_s = omega_s / pi;
F_ls = [0 freq_p freq_s 1];
A_ls = [1 1 0 0];

% Grila verificare
N = 5000;
[~, W_check] = freqz(1, 1, N);
idx_p = find(W_check <= omega_p, 1, 'last');
idx_s = find(W_check >= omega_s, 1, 'first');


% Metoda Ferestrei
[M_fir1, b_fir1] = cauta_ordin_min_Fereastra(M_but-1, freq_c, idx_p, idx_s, Delta_p, Delta_s, N);

% Metoda CMMP
[M_firls, b_firls] = cauta_ordin_min_CMMP(M_but-1, F_ls, A_ls, idx_p, idx_s, Delta_p, Delta_s, N);


% Calcul Raspunsuri in Frecventa 
[H_fir1, ~] = freqz(b_fir1, 1, N);
H_dB_fir1 = 20*log10(abs(H_fir1)+eps);
Phi_fir1 = unwrap(angle(H_fir1));

[H_firls, ~] = freqz(b_firls, 1, N);
H_dB_firls = 20*log10(abs(H_firls)+eps);
Phi_firls = unwrap(angle(H_firls));


% GRAFICE INDIVIDUALE FIR (Spectru + Faza pe aceeasi figura)

% Metoda Fereastra
figure('Name', 'Faza 2b - FIR Fereastra', 'Color', 'w', 'Position', [100 100 800 600]);

% Subplot 1: Spectru
subplot(2,1,1);
plot(W/pi, H_dB_fir1, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere', 'LabelVerticalAlignment', 'bottom');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere', 'LabelVerticalAlignment', 'top');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom');
text(0.05, -50, sprintf('Ordin M = %d', M_fir1), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Metoda Ferestrei'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-80 5]);
legend([h_wp, h_ws], 'Location', 'southwest');

% Subplot 2: Faza
subplot(2,1,2);
plot(W/pi, Phi_fir1, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda Ferestrei (M=' num2str(M_fir1) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)'); xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faza2b_fir1_individual.png');


% Metoda CMMP
figure('Name', 'Faza 2b - FIR CMMP', 'Color', 'w', 'Position', [100 100 800 600]);

% Subplot 1: Spectru
subplot(2,1,1);
plot(W/pi, H_dB_firls, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere', 'LabelVerticalAlignment', 'bottom');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere', 'LabelVerticalAlignment', 'top');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom');
text(0.05, -50, sprintf('Ordin M = %d', M_firls), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Metoda CMMP'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-80 5]);
legend([h_wp, h_ws], 'Location', 'southwest');

% Subplot 2: Faza
subplot(2,1,2);
plot(W/pi, Phi_firls, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda CMMP (M=' num2str(M_firls) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)'); xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faza2b_firls_individual.png');


% MATRICE 2x2 SPECTRE (TOATE 4 FILTRELE)
figure('Name', 'Faza 2b - Matrice Spectre', 'Color', 'w', 'Position', [50 50 1200 900]);

% 1. BUTTERWORTH
subplot(2,2,1);
plot(W/pi, H_dB_but, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'DisplayName', '\omega_p'); h_ws = xline(omega_s/pi, '--m', 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
yline(lim_pass_low, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'top', 'FontSize', 7);
yline(lim_stop, 'r', 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
text(0.05, -50, sprintf('M = %d', M_but), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');
title('Spectru Filtru Butterworth (IIR)'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-80 5]); legend([h_wp, h_ws], 'Location', 'southwest');

% 2. ELIPTIC
subplot(2,2,2);
plot(W/pi, H_dB_ellip, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'DisplayName', '\omega_p'); h_ws = xline(omega_s/pi, '--m', 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
yline(lim_pass_low, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'top', 'FontSize', 7);
yline(lim_stop, 'r', 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
text(0.05, -50, sprintf('M = %d', M_ellip), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');
title('Spectru Filtru Eliptic/Cauer (IIR)'); xlim([0 1]); ylim([-80 5]); legend([h_wp, h_ws], 'Location', 'southwest');

% 3. FEREASTRA
subplot(2,2,3);
plot(W/pi, H_dB_fir1, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'DisplayName', '\omega_p'); h_ws = xline(omega_s/pi, '--m', 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
yline(lim_pass_low, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'top', 'FontSize', 7);
yline(lim_stop, 'r', 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
text(0.05, -50, sprintf('M = %d', M_fir1), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');
title('Spectru Metoda Ferestrei (FIR)'); xlabel('Frecventa (\times\pi)'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-80 5]); legend([h_wp, h_ws], 'Location', 'southwest');

% 4. CMMP
subplot(2,2,4);
plot(W/pi, H_dB_firls, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'DisplayName', '\omega_p'); h_ws = xline(omega_s/pi, '--m', 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
yline(lim_pass_low, 'g', 'Label', 'Tol. Pass', 'LabelVerticalAlignment', 'top', 'FontSize', 7);
yline(lim_stop, 'r', 'Label', 'Tol. Stop', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
text(0.05, -50, sprintf('M = %d', M_firls), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');
title('Filtru Metoda CMMP (FIR)'); xlabel('Frecventa (\times\pi)'); xlim([0 1]); ylim([-80 5]); legend([h_wp, h_ws], 'Location', 'southwest');

print('-dpng', '-r300', '../Figuri/faza2b_matrice_spectre.png');


% MATRICE 2x2 FAZE (TOATE 4 FILTRELE)
figure('Name', 'Faza 2b - Matrice Faze', 'Color', 'w', 'Position', [50 50 1200 900]);

% 1. BUTTERWORTH
subplot(2,2,1);
plot(W/pi, Phi_but, 'k', 'LineWidth', 1.2); grid on;
title(['Faza Filtru Butterworth (M=' num2str(M_but) ')']); ylabel('Faza (rad)'); xlim([0 1]);

% 2. ELIPTIC
subplot(2,2,2);
plot(W/pi, Phi_ellip, 'b', 'LineWidth', 1.2); grid on;
title(['Faza Filtru Eliptic (M=' num2str(M_ellip) ')']); xlim([0 1]);

% 3. FEREASTRA
subplot(2,2,3);
plot(W/pi, Phi_fir1, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda Fereastra (M=' num2str(M_fir1) ')']); xlabel('Frecventa (\times\pi)'); ylabel('Faza (rad)'); xlim([0 1]);

% 4. CMMP
subplot(2,2,4);
plot(W/pi, Phi_firls, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Metoda CMMP (M=' num2str(M_firls) ')']); xlabel('Frecventa (\times\pi)'); xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faza2b_matrice_faze.png');

%% Subpunctul c: Filtrele Cebîsev si Comparatie Globala (6 Filtre).

% Cebîsev Tip 1 (cheby1) - Cautare ordin minim
% Input: freq_p (normalizat).
[M_c1, b_c1, a_c1] = cauta_ordin_min_Cheby1(1, Rp_dB, freq_p, idx_s, Delta_s, N);

% Cebîsev Tip 2 (cheby2) - Cautare ordin minim
% Input: freq_s (normalizat).
[M_c2, b_c2, a_c2] = cauta_ordin_min_Cheby2(1, Rs_dB, freq_s, idx_p, Delta_p, N);

% Calcul Raspunsuri in Frecventa (Cebîsev)
[H_c1, ~] = freqz(b_c1, a_c1, N);
H_dB_c1 = 20*log10(abs(H_c1)+eps);
Phi_c1 = unwrap(angle(H_c1));

[H_c2, ~] = freqz(b_c2, a_c2, N);
H_dB_c2 = 20*log10(abs(H_c2)+eps);
Phi_c2 = unwrap(angle(H_c2));

% GRAFICE CEBÎSEV (Spectru + Faza)

% CHEBY 1
figure('Name', 'Faza 2c - Cheby1 Individual', 'Color', 'w', 'Position', [100 100 800 600]);
% Spectru
subplot(2,1,1);
plot(W/pi, H_dB_c1, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', 'LabelVerticalAlignment', 'top');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', 'LabelVerticalAlignment', 'bottom');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', 'LabelVerticalAlignment', 'bottom');
text(0.05, -40, sprintf('Ordin M = %d', M_c1), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Cebîsev Tip 1'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-70 10]);
legend([h_wp, h_ws], 'Location', 'southwest');

% Faza
subplot(2,1,2);
plot(W/pi, Phi_c1, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Cebîsev Tip 1 (M=' num2str(M_c1) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)'); xlim([0 1]);
print('-dpng', '-r300', '../Figuri/faza2c_cheby1.png');

% CHEBY 2 
figure('Name', 'Faza 2c - Cheby2 Individual', 'Color', 'w', 'Position', [150 150 800 600]);
% Spectru
subplot(2,1,1);
plot(W/pi, H_dB_c2, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on; grid on;
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', 'LabelVerticalAlignment', 'top');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', 'LabelVerticalAlignment', 'bottom');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', 'LabelVerticalAlignment', 'bottom');
text(0.05, -40, sprintf('Ordin M = %d', M_c2), 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
title('Spectru Cebîsev Tip 2'); ylabel('Ampl (dB)'); xlim([0 1]); ylim([-70 10]);
legend([h_wp, h_ws], 'Location', 'southwest');

% Faza
subplot(2,1,2);
plot(W/pi, Phi_c2, 'r', 'LineWidth', 1.2); grid on;
title(['Faza Cebîsev Tip 2 (M=' num2str(M_c2) ')']);
xlabel('Frecventa normalizata (\times\pi)'); ylabel('Faza (rad)'); xlim([0 1]);
print('-dpng', '-r300', '../Figuri/faza2c_cheby2.png');


% MATRICE 2x3 SPECTRE (TOATE 6 FILTRELE)

% Colectare date intr-o structura pentru iterare usoara
Titles = {'1. Butterworth (IIR)', '2. Metoda Ferestrei (FIR)', '3. Cebîsev Tip 1 (IIR)', ...
          '4. Eliptic (IIR)', '5. Metoda CMMP (FIR)', '6. Cebîsev Tip 2 (IIR)'};
M_vals = [M_but, M_fir1, M_c1, M_ellip, M_firls, M_c2];
H_dB_vals = {H_dB_but, H_dB_fir1, H_dB_c1, H_dB_ellip, H_dB_firls, H_dB_c2};
Phi_vals  = {Phi_but, Phi_fir1, Phi_c1, Phi_ellip, Phi_firls, Phi_c2};

figure('Name', 'Faza 2c - Matrice Spectre 6 Filtre', 'Color', 'w', 'Position', [50 50 1400 900]);

for k = 1:6
    subplot(2, 3, k);
    % Plot Semnal
    plot(W/pi, H_dB_vals{k}, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
    hold on; grid on;
    
    % Linii Verticale
    h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
    h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');
    
    % Linii Orizontale
    yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1+\Delta_p)', ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
    yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1-\Delta_p)', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
    yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
    % Casuta Ordin
    text(0.05, -60, sprintf('Ordin M = %d', M_vals(k)), ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, 'FontWeight', 'bold');
    
    title(Titles{k});
    xlim([0 1]); ylim([-65 10]);
    
    % Axe labels doar pe margini pentru curatenie
    if k >= 4, xlabel('Frecventa normalizata (\times\pi)'); end
    if mod(k,3)==1, ylabel('Amplitudine (dB)'); end
    
    legend([h_wp, h_ws], 'Location', 'best');
end

print('-dpng', '-r300', '../Figuri/faza2c_matrice_spectre.png');

% MATRICE 2x3 FAZE (TOATE 6 FILTRELE)

figure('Name', 'Faza 2c - Matrice Faze 6 Filtre', 'Color', 'w', 'Position', [50 50 1400 900]);

for k = 1:6
    subplot(2, 3, k);
    
    % Plot Faza (FIR cu rosu, IIR cu albastru/negru pt variatie)
    color_line = 'b'; 
    if k == 2 || k == 5, color_line = 'r'; end % FIR sunt pe pozitiile 2 si 5
    
    plot(W/pi, Phi_vals{k}, color_line, 'LineWidth', 1.2); 
    grid on;
    
    title([Titles{k} ' (M=' num2str(M_vals(k)) ')']);
    xlim([0 1]);
    
    if k >= 4, xlabel('Frecventa normalizata (\times\pi)'); end
    if mod(k,3)==1, ylabel('Faza (rad)'); end
end

print('-dpng', '-r300', '../Figuri/faza2c_matrice_faze.png');