%% Izhikevich network - 6 conditions of gE/gI weights

clear; clc; close all;

% Network setup 
N  = 100;
Ne = round(0.8*N);
Ni = round(0.2*N);
dt = 1; T  = 1000; nt = T/dt;

% Izhikevich parameters 
a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)];
b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)];
c = -65*ones(N,1);
d = [8*ones(Ne,1); 2*ones(Ni,1)];

% Connectivity
S = [0.5*rand(N,Ne), -1*rand(N,Ni)];

% Synapse parameters (same for all)
hE = zeros(N,1);
hI = zeros(N,1);
tauE = 5; tauI = 10;
EE = 0; EI = -70;

% 6 conditions: [gE, gI]
conditions = [
    0.05, 0.10;  % Low E
    0.10, 0.10;  % Baseline  
    0.15, 0.10;  % High E
    0.10, 0.20;  % Baseline E, med I
    0.10, 0.40;  % Baseline E, high I
    0.15, 0.40   % High E+I
];

% Create figure with 6 subplots
figure('Color','k', 'Position', [100 100 1200 800]);

for cond = 1:6
    
    % Set gE, gI for this condition
    gE = conditions(cond,1);
    gI = conditions(cond,2);
    
    % Reset states
    V = -65*ones(N,1); u = b.*V;
    hE = zeros(N,1); hI = zeros(N,1);
    firings = [];
    
    % Simulation loop 
    for t = 1:nt
        I_ext = [5*randn(Ne,1); 2*randn(Ni,1)];
        IE = gE .* hE .* (EE - V);
        II = gI .* hI .* (EI - V);
        I  = I_ext + IE + II;
        
        dV = 0.04*V.^2 + 5*V + 140 - u + I;
        du = a.*(b.*V - u);
        V  = V + dt*dV; u  = u + dt*du;
        
        fired = find(V >= 30);
        Sspike = zeros(N,1); Sspike(fired) = 1;
        SE = S(:,1:Ne) * Sspike(1:Ne);
        SI = S(:,Ne+1:end) * Sspike(Ne+1:end);
        hE = hE + dt*(-hE + SE)/tauE;
        hI = hI + dt*(-hI + SI)/tauI;
        
        if ~isempty(fired)
            firings = [firings; t*ones(length(fired),1), fired];
            V(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
        end
    end
    
    % Plot subplot
    subplot(2,3,cond);
    if ~isempty(firings)
        plot(firings(:,1), firings(:,2), '.', 'Color',[0 1 1]);
    end
    ylim([0 N]);
    title(sprintf('gE=%.2f, gI=%.2f', gE, gI), 'Color', 'w', 'FontSize', 10);
    xlabel('Time (ms)');        
    ylabel('Neuron');           
end

sgtitle('Effect of synaptic weights on synchronisation', 'Color', 'w', 'FontSize', 14);
