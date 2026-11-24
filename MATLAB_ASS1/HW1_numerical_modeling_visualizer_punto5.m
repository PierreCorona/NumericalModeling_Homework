%% Generatore Animazione: Onde Stazionarie Sovrapposte
clear; clc; close all;

% Parametri spaziali e temporali (Senza smorzamento per vedere il moto puro)
x = linspace(0, 1, 1000);
t_end = 2.0; % Un periodo completo della fondamentale
dt = 0.02;   % Passo temporale dell'animazione
t_vec = 0:dt:t_end;

% Creazione figura
fig = figure('Color', 'w', 'Position', [100, 100, 800, 700]);
filename = 'animazione_onde.gif';

for i = 1:length(t_vec)
    t = t_vec(i);

    % Le equazioni fisiche delle onde stazionarie nel tempo
    % u(x,t) = sin(k*x) * cos(w*t)
    u_lenta = sin(pi*x) * cos(pi*t);               % Fondamentale (T = 2.0s)
    u_veloce = (1/8)*sin(8*pi*x) * cos(8*pi*t);    % Increspatura (T = 0.25s)
    u_tot = u_lenta + u_veloce;                    % Combinazione fisica

    % Plot 1: Onda Lenta (La pancia principale)
    subplot(3,1,1);
    plot(x, u_lenta, 'b', 'LineWidth', 2);
    axis([0 1 -1.2 1.2]); grid on;
    title(sprintf('1. Onda Fondamentale Lenta (T = 2.0s) | Tempo: %.2f s', t));
    ylabel('u_1');

    % Plot 2: Onda Veloce (Le increspature)
    subplot(3,1,2);
    plot(x, u_veloce, 'r', 'LineWidth', 1.5);
    axis([0 1 -0.2 0.2]); grid on;
    title('2. Increspature Veloci (T = 0.25s)');
    ylabel('u_2');

    % Plot 3: Onda Combinata (La vera corda)
    subplot(3,1,3);
    plot(x, u_tot, 'k', 'LineWidth', 2);
    axis([0 1 -1.2 1.2]); grid on;
    title('3. Onda Combinata sulla Corda (u_1 + u_2)');
    xlabel('Posizione x'); ylabel('u_{tot}');

    drawnow;

    % Salvataggio del frame per esportare la GIF
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', dt);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', dt);
    end
end

fprintf('Finito! Animazione salvata nella cartella corrente come %s\n', filename);