%% NUMERICAL MODELING AND SIMULATION FOR ACOUSTICS - Homework 1 (Point 5)
clear; clc; close all;

% 1. Parametri fisici e di discretizzazione
L = 1; 
rho = 1; 
mu = 1; 
T = 4;
Nx = 160; 
dt = 0.0025;

h = L / Nx;
Nt = round(T / dt);

% Vettore spaziale e nodi interni
x = linspace(0, L, Nx+1)';
x_int = x(2:end-1);
Nh = length(x_int);

% 2. Assemblaggio Matrici FEM P1
e = ones(Nh, 1);
M = (rho * h / 6) * spdiags([e, 4*e, e], -1:1, Nh, Nh);
K = (mu / h)      * spdiags([-e, 2*e, -e], -1:1, Nh, Nh);

% 3. Condizioni iniziali
% u_0(x) = sin(pi*x) + 1/8 * sin(8*pi*x)
u0_int = sin(pi * x_int) + (1/8) * sin(8 * pi * x_int);
v0_int = zeros(Nh, 1);

% 4. Definizione dei 4 casi di smorzamento (alpha, eta)
cases = struct(...
    'name', {'Undamped', 'Telegraph damping', 'Kelvin-Voigt damping', 'Combined damping'}, ...
    'alpha', {0, 0.8, 0, 0.8}, ...
    'eta', {0, 0, 0.02, 0.02} ...
);
num_cases = length(cases);

% Variabili per il salvataggio dei risultati
t_vec = linspace(0, T, Nt+1);
Energy_norm = zeros(Nt+1, num_cases);

% Istanti di tempo per i profili spaziali: t = 0.5, 1, 2
t_plot = [0.5, 1, 2];
steps_plot = round(t_plot / dt);
U_profiles = zeros(Nx+1, length(t_plot), num_cases); 

fprintf('Avvio simulazioni per i 4 casi...\n');

% 5. Ciclo principale sui 4 casi
for c = 1:num_cases
    alpha = cases(c).alpha;
    eta = cases(c).eta;
    
    fprintf('Calcolo Caso %d: %s (alpha=%.1f, eta=%.2f)...\n', c, cases(c).name, alpha, eta);
    
    % Setup Solutore Strang Splitting
    D = alpha * M + eta * K;
    B = M \ D;
    Q = expm(- (dt / 2) * B);
    
    A_lhs = M + (dt^2 / 4) * K;
    A_rhs = M - (dt^2 / 4) * K;
    
    % Inizializzazione per il caso corrente
    U = u0_int;
    V = v0_int;
    
    % Calcolo Energia Iniziale E_h(0)
    E0 = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
    Energy_norm(1, c) = 1.0; % E_h(0)/E_h(0)
    
    % Integrazione temporale
    for n = 1:Nt
        % Step solver
        Va = Q * V;
        rhs = A_rhs * U + dt * (M * Va);
        U_next = A_lhs \ rhs;
        Vb = (2 / dt) * (U_next - U) - Va;
        V = Q * Vb;
        U = U_next;
        
        % Calcolo Energia al tempo n
        En = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
        Energy_norm(n+1, c) = En / E0;
        
        % Salvataggio profili spaziali se siamo agli step richiesti
        [is_plot_step, idx] = ismember(n, steps_plot);
        if is_plot_step
            % Aggiungo gli zeri ai bordi (condizioni di Dirichlet)
            U_profiles(:, idx, c) = [0; U; 0];
        end
    end
end
fprintf('Simulazioni completate. Generazione grafici...\n');

% =========================================================================
% 6. GRAFICI
% =========================================================================

% Grafico 1: Energia Normalizzata
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
colors = {'k', 'b', 'r', 'g'};
for c = 1:num_cases
    plot(t_vec, Energy_norm(:, c), 'LineWidth', 2, 'Color', colors{c}, 'DisplayName', cases(c).name);
end
xlabel('Tempo (t)');
ylabel('Energia Normalizzata E_h(t) / E_h(0)');
title('Decadimento dell''Energia per i 4 casi di smorzamento');
legend('Location', 'northeast');
set(gca, 'FontSize', 12);
ylim([0 1.05]);

% Grafico 2: Profili Spaziali (Subplots per ogni caso)
figure('Color', 'w', 'Position', [150, 150, 1000, 800]);
t_colors = {'#D95319', '#EDB120', '#7E2F8E'}; % Colori per t=0.5, t=1, t=2
for c = 1:num_cases
    subplot(2, 2, c);
    hold on; grid on;
    
    % Plot condizione iniziale (t=0) come riferimento
    plot(x, [0; u0_int; 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 't = 0 (IC)');
    
    % Plot istanti richiesti
    for idx = 1:length(t_plot)
        plot(x, U_profiles(:, idx, c), 'LineWidth', 1.5, 'Color', t_colors{idx}, ...
             'DisplayName', sprintf('t = %.1f', t_plot(idx)));
    end
    
    title(cases(c).name);
    xlabel('Spazio (x)');
    ylabel('Ampiezza u_h(x,t)');
    xlim([0 1]); ylim([-1.2 1.2]);
    if c == 1
        legend('Location', 'southwest', 'FontSize', 9);
    end
end
sgtitle('Profili spaziali dell''onda nel tempo u_h(x,t)');