%% Izhikevich network - tauE/tauI variation (fixed gE=0.15, gI=0.10)

clear; clc; close all;

% Network setup 
N  = 100; Ne = round(0.8*N); Ni = round(0.2*N);
dt = 1; T  = 1000; nt = T/dt;

% Fixed Izhikevich + synapse parameters
a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)]; b = 0.20*ones(N,1);
c = -65*ones(N,1); d = [8*ones(Ne,1); 2*ones(Ni,1)];
S = [0.5*rand(N,Ne), -1*rand(N,Ni)];
gE = 0.15; gI = 0.10; EE = 0; EI = -70;  % FIXED

% 6 conditions: [tauE, tauI]
conditions = [
    2,  5;   % Fast E/I
    5, 10;   % Baseline  
   10, 20;   % Slow E/I
    1,  2;   % Very fast
   20, 40;   % Very slow
    3, 15    % Mixed
];

% Create figure
figure('Color','k', 'Position', [100 100 1200 800]);

for cond = 1:6
    tauE = conditions(cond,1); tauI = conditions(cond,2);
    
    % Reset states
    V = -65*ones(N,1); u = b.*V; hE = zeros(N,1); hI = zeros(N,1); firings = [];
    
    % Simulation loop 
    for t = 1:nt
        I_ext = [5*randn(Ne,1); 2*randn(Ni,1)];
        IE = gE .* hE .* (EE - V); II = gI .* hI .* (EI - V); I = I_ext + IE + II;
        dV = 0.04*V.^2 + 5*V + 140 - u + I; du = a.*(b.*V - u);
        V = V + dt*dV; u = u + dt*du;
        
        fired = find(V >= 30); Sspike = zeros(N,1); Sspike(fired) = 1;
        SE = S(:,1:Ne)*Sspike(1:Ne); SI = S(:,Ne+1:end)*Sspike(Ne+1:end);
        hE = hE + dt*(-hE + SE)/tauE;    % ← tauE changes!
        hI = hI + dt*(-hI + SI)/tauI;    % ← tauI changes!
        
        if ~isempty(fired)
            firings = [firings; t*ones(length(fired),1), fired];
            V(fired) = c(fired); u(fired) = u(fired) + d(fired);
        end
    end
    
    % Plot subplot
    subplot(2,3,cond); cla;
    if ~isempty(firings), plot(firings(:,1), firings(:,2), '.', 'Color',[0 1 1]); end
    ylim([0 N]); title(sprintf('τ_E=%.0f, τ_I=%.0f', tauE, tauI), 'Color','w','FontSize',10);
    xlabel('Time (ms)'); ylabel('Neuron');
end

sgtitle('Effect of synaptic time constants on synchrony (gE=0.15, gI=0.10)', 'Color','w','FontSize',14);
