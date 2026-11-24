%% NUMERICAL MODELING AND SIMULATION FOR ACOUSTICS - Homework 1 (Point 4)
clear; clc; close all;

% 1. Parametri fisici e costanti analitiche
L = 1; T = 1;
rho = 1; mu = 1;
alpha = 0.8; eta = 0.02;

% Parametri ausiliari per la soluzione esatta
d_exact = alpha + eta * pi^2;
beta = sqrt(pi^2 - d_exact^2 / 4);

% Soluzione esatta u_ex(x,t) definita come anonymous function
uex_func = @(x, t) sin(pi * x) .* exp(-d_exact * t / 2) .* ...
                   (cos(beta * t) + (d_exact / (2 * beta)) * sin(beta * t));

% 2. Setup delle risoluzioni richieste dal problema
N_list = [20, 40, 80, 160]; %faccio un ciclo for ogni volta per un numero di elementi che compongono una mesh diverso
num_runs = length(N_list); %l'errore viene calcolato esclusivamente nell'istante finale t=T. Non è una somma o un accumulo degli errori calcolati negli istanti temporali precedenti. Prendiamo l'errore al punto t=T perché in quell'istante si trova l'errore totale di quell'analisi perché gli errori in ogni istante temporale vengono accumulati.


% Vettori per salvare i risultati
h_list = zeros(num_runs, 1); %in h_list stai salvando la lunghezza della mesh spaziale
errors_L2 = zeros(num_runs, 1); %

fprintf('Avvio analisi di convergenza...\n');
fprintf('----------------------------------------------------------\n');
fprintf('   Nx=Nt   |      h      |    L2 Error    |  Conv. Rate  |\n');
fprintf('----------------------------------------------------------\n');

% 3. Ciclo sulle 4 configurazioni di mesh
for i = 1:num_runs
    N = N_list(i);
    
    % Parametri di discretizzazione
    h = L / N; 
    dt = T / N;
    h_list(i) = h;
    
    % Creazione della mesh (nodi interni)
    x = linspace(0, L, N+1)';
    x_int = x(2:end-1); %Tuttavia, poiché le condizioni al contorno ci dicono che agli estremi la soluzione vale sempre zero in ogni istante di tempo, noi conosciamo già il risultato sui bordi. I valori veramente "incogniti" che il risolutore deve trovare sono solo quelli interni.
    Nh = length(x_int); %Nh calcola la lunghezza del vettore x_int, indicando il numero esatto di variabili incognite da risolvere a ogni passo temporale.
    
    % Assemblaggio Matrici FEM P1
    e = ones(Nh, 1);
    M = (rho * h / 6) * spdiags([e, 4*e, e], -1:1, Nh, Nh); %matrice tridiagonale: il comando (-1,1) distribuisce i valori sulla sotto-diagonale, diagonale e sopra-diagonale rispettivamente
    K = (mu / h)      * spdiags([-e, 2*e, -e], -1:1, Nh, Nh);
    
    % Condizioni iniziali per la simulazione numerica
    u0 = sin(pi * x_int);
    v0 = zeros(Nh, 1);
    
    % 4. Esecuzione del solutore 
    % Per N_t step temporali arriviamo esattamente a t = T
    U_final = solve_wave_splitting_final_state(full(M), full(K), u0, v0, dt, N, alpha, eta);
    
    % 5. Calcolo della soluzione esatta al tempo finale T
    U_exact_final = uex_func(x_int, T);
    
    % 6. Calcolo dell'errore in norma L2 pesato con la matrice di massa
    err_vec = U_exact_final - U_final; %errore calcolato all'istante t=T
    errors_L2(i) = sqrt(err_vec' * M * err_vec);
    
    % Calcolo del convergence rate (salto il primo step perché non ho un precedente)
    if i == 1 
        rate = NaN; %Al primissimo giro del ciclo (quando $N = 20$), il computer calcola l'errore per la prima volta. Per capire di quanto è migliorato l'errore, avresti bisogno di confrontarlo con un'analisi precedente, che però non esiste.
        fprintf(' %7d   |  %1.4e |  %1.4e  |      -       |\n', N, h, errors_L2(i));
    else
        rate = log2(errors_L2(i-1) / errors_L2(i)); %il numero di elementi raddoppia a ogni step. Questo significa che il passo spaziale h (che è L/N) viene esattamente dimezzato a ogni iterazione.
        fprintf(' %7d   |  %1.4e |  %1.4e  |    %1.4f    |\n', N, h, errors_L2(i), rate); %Se dimezzo h, il nuovo errore dovrebbe diventare 2^p volte più piccolo
        %l'approssimazione di una curva continua tramite segmenti dritti genera un errore proporzionale al quadrato della larghezza del segmento: quindi il rate di convergenza rimane sempre costante a 2 (circa) 
    end
end
fprintf('----------------------------------------------------------\n');

% 7. Plot Log-Log dell'errore contro h
figure('Color', 'w', 'Position', [200, 200, 600, 500]);
loglog(h_list, errors_L2, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '||u_{ex} - u_h||_{L^2}');
hold on;
grid on;
%perchè loglog? Se disegnassi questa curva su un grafico normale, vedresti una parabola schiacciata verso lo zero. 
% Ma se applichi il logaritmo a destra e a sinistra, la formula diventa:
% log(E) = p*log(h) + log(C): Questa è l'equazione di una retta, in cui la pendenza (l'inclinazione) è esattamente l'ordine di convergenza p.

% Creazione di un triangolo di pendenza teorica O(h^2) per confronto visivo
h_ref = [h_list(2), h_list(3)];
err_ref = [errors_L2(2), errors_L2(2) * (h_ref(2)/h_ref(1))^2] * 1.2; % Shiftato leggermente in alto per visibilità
loglog(h_ref, err_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h^2) Reference');

xlabel('Dimensione maglia h (Log)');
ylabel('Errore L^2 al tempo T=1 (Log)');
title('Analisi di convergenza spaziale-temporale');
legend('Location', 'northwest');
set(gca, 'FontSize', 12);


%% ========================================================================
%  FUNZIONE DEL SOLUTORE: Restituisce SOLO lo stato finale per efficienza
%  ========================================================================
function U = solve_wave_splitting_final_state(M, K, u0, v0, dt, Nt, alpha, eta)
    % 1. Definizione D e B
    D = alpha * M + eta * K;
    B = M \ D;
    
    % 2. Calcolo matrice esponenziale Q
    Q = expm(- (dt / 2) * B);
    
    % 3. Matrici del passo centrale Crank-Nicolson
    A_lhs = M + (dt^2 / 4) * K;
    A_rhs = M - (dt^2 / 4) * K;
    
    % 4. Inizializzazione
    U = u0;
    V = v0;
    
    % 5. Integrazione temporale
    for n = 1:Nt
        Va = Q * V;
        rhs = A_rhs * U + dt * (M * Va);
        U_next = A_lhs \ rhs;
        Vb = (2 / dt) * (U_next - U) - Va;
        V = Q * Vb;
        U = U_next; %all'uscita della funzione hai il displacement dell'onda calcolata al tempo T contente principalmente gli errori di amplitude dissipation (dovuto allo smorzamento) accumulati negli istanti precedenti)
    end
end