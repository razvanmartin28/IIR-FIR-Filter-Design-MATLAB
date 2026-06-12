% FAZA 4: Supliment - Proiectare Butterworth cu castig neunitar
% Autor: Martin Razvan-Stefan
clear; clc; close all;

% --- Datele specifice (Faza 1a) ---
ns = 9;
ng = 2;

% Valorile specifice salvate de la Faza 1a
omega_s = 1.6122;
omega_p = 1.3721;
Delta_p = 0.0776;
Delta_s = Delta_p;
Ts = 2.5784;

% Afisare date initiale
fprintf('Date de proiectare (H(0) = 1 + Delta_p):\n');
fprintf('  wp = %.4f rad (%.2f pi)\n', omega_p, omega_p/pi);
fprintf('  ws = %.4f rad (%.2f pi)\n', omega_s, omega_s/pi);
fprintf('  Delta_p = %.4f\n', Delta_p);
fprintf('  Ts = %.2f s\n', Ts);

% Proiectarea filtrului folosind But_FTI_Faza4
[B, A] = But_FTI_Faza4(omega_p/pi, omega_s/pi, Delta_p, Delta_s, Ts);


% Calculul parametrilor M si omega_c 
wp_analog = 2*tan(omega_p/2)/Ts;
ws_analog = 2*tan(omega_s/2)/Ts;
num_term_M = ( (1 + Delta_p + Delta_s) * (1 + Delta_p - Delta_s) ) / ( Delta_p * (Delta_s^2) * (2 + Delta_p) );
M = ceil( log(num_term_M) / (2 * log(ws_analog/wp_analog)) );

% Calcul Omega_c (analogic)
term_c = Delta_p * (2 + Delta_p);
Omega_c = wp_analog / ( term_c^(1/(2*M)) );

% Calcul omega_c (discret)
omega_c = 2*atan(Omega_c*Ts/2);

fprintf('\nRezultate obtinute:\n');
fprintf('  Ordinul filtrului M = %d\n', M);
fprintf('  Pulsatia de taiere wc = %.4f rad (%.4f pi)\n', omega_c, omega_c/pi);
fprintf('  Castig la DC (H(0) teoretic) = %.4f\n', 1 + Delta_p);


% Trasarea caracteristicilor
N_points = 5000;
[H, W] = freqz(B, A, N_points);

% Calcul amplitudine liniara pentru a verifica H(0)
H_abs = abs(H);
fprintf('  Castig la DC (H(0) masurat)  = %.4f\n', H_abs(1));

H_dB = (H_abs);
Phi = unwrap(angle(H));

%Spectrul Filtrului Butterworth
figure('Name', 'Faza 4 - Spectru Amplitudine', 'Color', 'w', 'Position', [100 100 800 600]);
plot(W/pi, H_dB, 'b', 'LineWidth', 1.5); hold on;
grid on;

% Linii verticale (pulsatiile)
h_wp = xline(omega_p/pi, '--r', 'LineWidth', 1, 'DisplayName', '\omega_p');
h_wc = xline(omega_c/pi, '-.k', 'LineWidth', 1, 'DisplayName', '\omega_c');
h_ws = xline(omega_s/pi, '--m', 'LineWidth', 1, 'DisplayName', '\omega_s');

% Linii orizontale 
% Toleranta de trecere superioara: 1 + 2*Delta_p
lim_pass_high = (1 + 2*Delta_p);
% Toleranta de trecere inferioara: 1 
lim_pass_low = 1; 
% Toleranta de stop: Delta_s
lim_stop = (Delta_s);

% Plotare limite
yline(lim_pass_high, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1 + 2*\Delta_p)', 'LabelVerticalAlignment', 'top');
yline(lim_pass_low, 'g', 'LineWidth', 1, 'Label', 'Tol. Trecere (1)', 'LabelVerticalAlignment', 'bottom');
yline(lim_stop, 'r', 'LineWidth', 1, 'Label', 'Tol. Stop (\Delta_s)', 'LabelVerticalAlignment', 'bottom');

info_txt = sprintf('M = %d\nH(0) = 1 + \\Delta_p', M);
text(0.05, -30, info_txt, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

title('Spectrul Filtrului Butterworth (Câștig Neunitar)');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Amplitudine (dB)');
legend([h_wp, h_wc, h_ws], 'Location', 'best');
xlim([0 1]);
ylim([-50 5]);

print('-dpng', '-r300', '../Figuri/spectru_faza4.png');

% Faza Filtrului Butterworth 
figure('Name', 'Faza 4 - Faza', 'Color', 'w');
plot(W/pi, Phi, 'r', 'LineWidth', 1.5);
grid on;
title('Faza Filtrului Butterworth (Câștig neunitar)');
xlabel('Frecventa normalizata (\times\pi rad/sample)');
ylabel('Faza (rad)');
xlim([0 1]);

print('-dpng', '-r300', '../Figuri/faza_faza4.png');

