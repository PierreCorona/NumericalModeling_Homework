%% NUMERICAL MODELING AND SIMULATION FOR ACOUSTICS - Homework 1 (Point 6 - Newmark)
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

% Parametri per il caso Combined Damping
alpha = 0.8;
eta = 0.02;
D = alpha * M + eta * K;

% 3. Condizioni iniziali
u0_int = sin(pi * x_int) + (1/8) * sin(8 * pi * x_int);
v0_int = zeros(Nh, 1);

% Vettori per il salvataggio dei risultati
t_vec = linspace(0, T, Nt+1);
Energy_Split = zeros(Nt+1, 1);
Energy_Newmark = zeros(Nt+1, 1);

% Istanti di tempo per i profili spaziali: t = 0.5, 1, 2
t_plot = [0.5, 1, 2];
steps_plot = round(t_plot / dt);
U_profiles_Split = zeros(Nx+1, length(t_plot)); 
U_profiles_Newmark = zeros(Nx+1, length(t_plot)); 

% =========================================================================
% SOLUTORE 1: STRANG SPLITTING (Dal Punto 5)
% =========================================================================
fprintf('Esecuzione Solutore Strang Splitting...\n');
B = M \ D;
Q = expm(- (dt / 2) * B);
A_lhs = M + (dt^2 / 4) * K;
A_rhs = M - (dt^2 / 4) * K;

U = u0_int;
V = v0_int;
E0 = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
Energy_Split(1) = 1.0;

for n = 1:Nt
    % Step solver
    Va = Q * V;
    rhs = A_rhs * U + dt * (M * Va);
    U_next = A_lhs \ rhs;
    Vb = (2 / dt) * (U_next - U) - Va;
    V = Q * Vb;
    U = U_next;
    
    % Calcolo Energia
    En = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
    Energy_Split(n+1) = En / E0;
    
    % Salvataggio profili
    [is_plot_step, idx] = ismember(n, steps_plot);
    if is_plot_step
        U_profiles_Split(:, idx) = [0; U; 0];
    end
end

% =========================================================================
% SOLUTORE 2: METODO DI NEWMARK (gamma = 1/2, beta = 1/4)
% =========================================================================
fprintf('Esecuzione Solutore Newmark (Average Acceleration)...\n');
gamma_N = 1/2;
beta_N = 1/4;

% Calcolo accelerazione iniziale a0
% Dalla EDO: M*a0 + D*v0 + K*u0 = 0  => a0 = M \ (-D*v0 - K*u0)
A0 = M \ (-D * v0_int - K * u0_int);

U_nm = u0_int;
V_nm = v0_int;
A_nm = A0;
Energy_Newmark(1) = 1.0;

% Matrice di iterazione del sistema lineare di Newmark (costante nel tempo)
LHS_Newmark = M + (gamma_N * dt) * D + (beta_N * dt^2) * K;

for n = 1:Nt
    % 1. Vettori predittori (per spostamento e velocità)
    U_pred = U_nm + dt * V_nm + ((0.5 - beta_N) * dt^2) * A_nm;
    V_pred = V_nm + ((1 - gamma_N) * dt) * A_nm;
    
    % 2. Risoluzione del sistema lineare per la nuova accelerazione
    RHS = - K * U_pred - D * V_pred;
    A_next = LHS_Newmark \ RHS;
    
    % 3. Correttore (Aggiornamento finale spostamento e velocità)
    U_next = U_pred + (beta_N * dt^2) * A_next;
    V_next = V_pred + (gamma_N * dt) * A_next;
    
    % 4. Calcolo Energia
    En_nm = 0.5 * (V_next' * M * V_next) + 0.5 * (U_next' * K * U_next);
    Energy_Newmark(n+1) = En_nm / E0;
    
    % Aggiornamento variabili per step successivo
    U_nm = U_next;
    V_nm = V_next;
    A_nm = A_next;
    
    % Salvataggio profili
    [is_plot_step, idx] = ismember(n, steps_plot);
    if is_plot_step
        U_profiles_Newmark(:, idx) = [0; U_nm; 0];
    end
end
fprintf('Simulazioni completate. Generazione grafici...\n');

% =========================================================================
% GRAFICI DI CONFRONTO
% =========================================================================
% Grafico 1: Confronto Decadimento Energia
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
plot(t_vec, Energy_Split, 'b', 'LineWidth', 2, 'DisplayName', 'Strang Splitting');
hold on; grid on;
plot(t_vec, Energy_Newmark, 'r--', 'LineWidth', 2, 'DisplayName', 'Newmark (\gamma=1/2, \beta=1/4)');
xlabel('Tempo (t)');
ylabel('Energia Normalizzata E_h(t) / E_h(0)');
title('Confronto Decadimento Energia (Combined Damping)');
legend('Location', 'northeast');
set(gca, 'FontSize', 12);
ylim([0 1.05]);

% Grafico 2: Confronto Profili Spaziali
figure('Color', 'w', 'Position', [150, 150, 1000, 400]);
t_colors = {'#D95319', '#EDB120', '#7E2F8E'};

subplot(1, 2, 1);
hold on; grid on;
plot(x, [0; u0_int; 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 't = 0 (IC)');
for idx = 1:length(t_plot)
    plot(x, U_profiles_Split(:, idx), 'LineWidth', 1.5, 'Color', t_colors{idx}, 'DisplayName', sprintf('t = %.1f', t_plot(idx)));
end
title('Strang Splitting');
xlabel('Spazio (x)'); ylabel('Ampiezza u_h(x,t)');
xlim([0 1]); ylim([-1.2 1.2]); legend('Location', 'southwest', 'FontSize', 9);

subplot(1, 2, 2);
hold on; grid on;
plot(x, [0; u0_int; 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 't = 0 (IC)');
for idx = 1:length(t_plot)
    plot(x, U_profiles_Newmark(:, idx), 'LineWidth', 1.5, 'Color', t_colors{idx}, 'DisplayName', sprintf('t = %.1f', t_plot(idx)));
end
title('Newmark (\gamma=1/2, \beta=1/4)');
xlabel('Spazio (x)'); ylabel('Ampiezza u_h(x,t)');
xlim([0 1]); ylim([-1.2 1.2]);
sgtitle('Confronto Profili Spaziali (Combined Damping)');

% Grafico 3: Differenza numerica di Energia
figure('Color', 'w', 'Position', [200, 200, 800, 400]);
plot(t_vec, Energy_Newmark - Energy_Split, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Tempo (t)');
ylabel('Differenza E_{Newmark} - E_{Splitting}');
title('Differenza Assoluta di Energia tra i due Solutori');
set(gca, 'FontSize', 12);