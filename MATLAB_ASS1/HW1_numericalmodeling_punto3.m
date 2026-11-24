%% NUMERICAL MODELING AND SIMULATION FOR ACOUSTICS - Homework 1 (Point 3)
clear; clc; close all;

%% 1. Parametri fisici e di discretizzazione
L = 1.0;            % Lunghezza del dominio spaziale (0, L)
T = 1;            % Tempo di analisi, ovvero il tempo finale di simulazione (0, T] ==> vediamo mezzo periodo dell'onda (che ha periodo T_onda= 2 sec)
rho = 1.0;          % Densità del mezzo (rho > 0)
mu = 1.0;           % Modulo elastico (mu > 0)

%dati presi dal punto 4

Nel = 100;          % Numero di elementi finiti arbitrario
h = L / Nel;        % Passo della griglia spaziale (mesh uniforme)
x = linspace(0, L, Nel+1)'; % Tutti i nodi della mesh
x_int = x(2:end-1); % Nodi interni (esclusi i bordi di Dirichlet x=0 e x=L)
Nh = length(x_int); % Numero di gradi di libertà (nodi interni)

dt = 0.005;         % Passo di integrazione temporale (Delta t) scelto arbitrariamente
Nt = round(T / dt); % Numero totale di step temporali (Nt = T / dt)
t = (0:Nt) * dt;    % Vettore dei tempi

%% 2. Condizioni iniziali (interpolazione nodale u0 e v0)
% Esempio: primo modo di vibrare con velocità iniziale nulla
u0_func = @(x) sin(pi * x / L);
v0_func = @(x) zeros(size(x));

u0 = u0_func(x_int); % Vettore colonna u0 valutato nei nodi interni
% il seno vale esattamente 0 agli estremi (x=0 e x=L)
v0 = v0_func(x_int); % Vettore colonna v0 valutato nei nodi interni
%v0=0: il sistema viene fatto partire da fermo

%% 3. Assemblaggio delle matrici FEM P1 (Massa M e Rigidezza K)
% Per elementi finiti lineari a tratti su griglia 1D uniforme:
% M_ij = int_0^L (rho * phi_j * phi_i) dx -> tridiagonale con [1, 4, 1] * (rho*h/6)
% K_ij = int_0^L (mu * phi_j' * phi_i') dx -> tridiagonale con [-1, 2, -1] * (mu/h)
% tridiagonale = tre elementi per riga

e = ones(Nh, 1);
M = (rho * h / 6) * spdiags([e, 4*e, e], -1:1, Nh, Nh);
K = (mu / h)      * spdiags([-e, 2*e, -e], -1:1, Nh, Nh);
%costruzione delle matrici di massa e stifness secondo la loro definizione
%(moltiplicate per le funzioni base)

% Conversione in matrici full per l'uso diretto di expm su sistemi piccoli
M_full = full(M);
K_full = full(K);
%in MATLAB la funzione expm accetta in ingresso solo matrici in formato full

%% 4. Simulazione Caso 1: Non smorzato (Undamped: alpha = 0, eta = 0)
alpha_undamped = 0;
eta_undamped   = 0;

[U_undamped, V_undamped, Eh_undamped] = solve_wave_splitting( ...
    M_full, K_full, u0, v0, dt, Nt, alpha_undamped, eta_undamped);

%% 5. Simulazione Caso 2: Smorzato (Damped: alpha > 0, eta > 0)
alpha_damped = 0.5;
eta_damped   = 0.002;

[U_damped, V_damped, Eh_damped] = solve_wave_splitting( ...
    M_full, K_full, u0, v0, dt, Nt, alpha_damped, eta_damped);

%% 6. Verifica numerica e Grafici dell'Energia Discreta Eh^n
figure('Color', 'w', 'Position', [100, 100, 900, 400]);

% Grafico 1: Confronto energia nel tempo
subplot(1, 2, 1);
plot(t, Eh_undamped, 'b-', 'LineWidth', 2, 'DisplayName', 'Undamped (\alpha = \eta = 0)');
hold on;
plot(t, Eh_damped, 'r--', 'LineWidth', 2, 'DisplayName', 'Damped (\alpha > 0, \eta > 0)');
grid on;
xlabel('Tempo t [s]');
ylabel('Energia discreta E_h^n');
title('Evoluzione dell''energia E_h^n');
legend('Location', 'best');

% Grafico 2: Errore relativo di conservazione nel caso non smorzato
subplot(1, 2, 2);
rel_err_energy = abs(Eh_undamped - Eh_undamped(1)) / Eh_undamped(1);
semilogy(t, rel_err_energy + eps, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Tempo t [s]');
ylabel('|E_h^n - E_h^0| / E_h^0');
title('Errore relativo conservazione (\alpha = \eta = 0)');

% Controllo numerico a terminale
max_energy_variation_undamped = max(abs(Eh_undamped - Eh_undamped(1)));
max_energy_increase_damped    = max(diff(Eh_damped));

fprintf('--- RISULTATI VERIFICA ENERGIA ---\n');
fprintf('Caso Undamped: Variazione massima di E_h = %.3e (Conservata!)\n', max_energy_variation_undamped);
fprintf('Caso Damped:   Incremento massimo tra due step (E_h^{n+1} - E_h^n) = %.3e (Non crescente!)\n', max_energy_increase_damped);


%% ========================================================================
%  FUNZIONE DEL SOLUTORE (SPLITTING SCHEME)
%  ========================================================================
function [U_hist, V_hist, Eh] = solve_wave_splitting(M, K, u0, v0, dt, Nt, alpha, eta)
    
    Nh = length(u0);
    
    % 1. Definizione della matrice di smorzamento D = alpha*M + eta*K
    D = alpha * M + eta * K;
    %alpha e eta definiscono lo smorzamento e formano la matrice
    
    % 2. Calcolo di B risolvendo il sistema lineare M*B = D
    B = M \ D;
    
    % 3. Calcolo dell'operatore esponenziale Q = exp(-dt/2 * M^{-1}*D)
    Q = expm(- (dt / 2) * B);
    
    % 4. Pre-calcolo delle matrici del passo centrale (propagazione dell'onda)
    %    Costanti nel tempo: calcolate fuori dal ciclo per efficienza
    A_lhs = M + (dt^2 / 4) * K; % Matrice a sinistra dell'uguale: (M + dt^2/4 * K)
    A_rhs = M - (dt^2 / 4) * K; % Matrice a destra dell'uguale:   (M - dt^2/4 * K)
    
    % 5. Inizializzazione dei vettori di stato e dello storico
    U = u0; % U^0 = u0
    V = v0; % V^0 = v0
    
    U_hist = zeros(Nh, Nt + 1);
    V_hist = zeros(Nh, Nt + 1);
    Eh     = zeros(1, Nt + 1);
    
    U_hist(:, 1) = U;
    V_hist(:, 1) = V;
    
    % Valutazione dell'energia iniziale E_h^0 (indice 1 in MATLAB)
    Eh(1) = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
    
    % 6. Ciclo temporale per n = 0, 1, ..., Nt - 1
    for n = 1:Nt
        % Passo 1: Primo mezzo step di smorzamento (Damping)
        Va = Q * V;
        
        % Passo 2: Aggiornamento spostamento U^{n+1} (Wave propagation)
        rhs = A_rhs * U + dt * (M * Va);
        U_next = A_lhs \ rhs;
        
        % Passo 3: Aggiornamento velocità intermedia V^b (Wave propagation)
        Vb = (2 / dt) * (U_next - U) - Va;
        
        % Passo 4: Secondo mezzo step di smorzamento per ottenere V^{n+1} (Damping)
        V_next = Q * Vb;
        
        % Aggiornamento variabili per lo step successivo
        U = U_next;
        V = V_next;
        
        % Salvataggio soluzioni e calcolo dell'energia discreta E_h^n
        U_hist(:, n + 1) = U;
        V_hist(:, n + 1) = V;
        Eh(n + 1) = 0.5 * (V' * M * V) + 0.5 * (U' * K * U);
    end
end