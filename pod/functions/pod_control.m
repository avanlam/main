function [ ] = pod_control(n_osc)
% POD_CONTROL orchestrates the entire optimisation problem applied to the
% modified Goodwin model, visualising the loss function and the gradient 
% descent

%n_osc = 100;

%% Control Parameters
n_samples = 7;                 	% 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
n_test_samples = 10;                % 2. Size of the test sample: total number of function evaluation will be N * (numDim + 2)
n_trials = 1;                       % 3. Number of realisations for the natural frequencies in the base sample
n_trials_samples = 50;              % 4. Number of realisations for the natural frequencies in the test sample
n_latent_pod=15;                    % 5. Reduced number of oscillators, up to N
n_samples_deim=10;                  % 6. Reduced number of oscillators, up to N
seedNum = 123456789;               	% 7. Seed for random number generator; put [] for autmatic randomization 
funPath = '../syst';                % 8. Folder address: that includes model file, factor space
funFile = 'eval_circ';            	% 9. Model/function file: MATLAB m-file without .m extension
colours = colour_pairs('springXL');
%% System Parameters
final_time=500;             % 1. Final time
n_steps=2*final_time;         % 2. Number of time steps
dt=final_time/n_steps;      % 3. Size of time step
n_asympt=n_steps/5;         % 4. Asymptotic time steps (used to compute QoIs)

x0=repmat([0.0998 0.2468 2.0151 0.0339].', n_osc, 1); 	% 5. Initial point on the stable limit cycle, for each neuron
a=0.8; b=1.2; tau=a+(b-a)*rand(n_osc,1);     % 6. Intrinsic periods: uniform distribution

params.nu=[0.7 0.35 NaN 0.35 NaN 0.35 NaN 1]';          % 7. Fixed parameters (Michaelis-Menten + Light)
params.Ki=[1 1 NaN 1 NaN 1 NaN 1]';
params.k=[NaN NaN 0.7 NaN 0.7 NaN 0.35]';
params.nuc=0.4;
params.Kc=1;
params.omega=2*pi/24;
params.K=0.3; params.L0=0.01;
%% Store POD Specifications
rom.n_osc = n_osc; 
rom.seedNum = seedNum;
rom.funFile = funFile;
rom.funPath = funPath;
rom.funFile = funFile;
rom.n_trials = n_trials;
time.final_time = final_time;
time.n_steps = n_steps;
time.n_asympt = n_asympt;
time.dt = dt;
time.x0 = x0;
%% Randomization
if isempty( seedNum ) == false
    rand('state',seedNum); 
end

%% Generate the original base sample as the single basal value
addpath(funPath);
factors = read_factorSpace('_pod'); % read in the factor space
lb = factors.lb; ub = factors.ub; 

dim_sample = length(lb);  % number of factors
rom.dim_sample = dim_sample;
rom.n_samples = 1;

tic;
base_sample_rom = 0.5.* ( ub' - lb') + lb';
chrono = toc;
     
fprintf(['.. Standard base sample generated in ', num2str(chrono),' seconds \n']);

%% Generate the parametric base sample from the unit hypercube

prom = rom;
prom.n_samples = n_samples;

tic;
p = sobolset(dim_sample); 
% ,'Skip',1e3,'Leap',1e2 %skip the first 1000 values, and then retain every 101st point
% p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
base_sample_prom = p(1:n_samples,:);

for j = 1 : n_samples
    base_sample_prom(j, :) = base_sample_prom(j, :) .* ( ub' - lb') + lb';
end
chrono = toc;
base_sample_prom(end,1) = 28;

fprintf(['.. Parametric base sample generated in ', num2str(chrono),' seconds \n']);

%% Generate the test sample as a grid
tic;
x_all=linspace(lb(1),ub(1),n_test_samples);
y_all=linspace(lb(2),ub(2),n_test_samples);
z_all=linspace(lb(3),ub(3),n_test_samples);

[x, y, z] = meshgrid(x_all, y_all, z_all);
test_sample = [ reshape(x, [n_test_samples^3, 1]), reshape(y, [n_test_samples^3, 1]), reshape(z, [n_test_samples^3, 1])];
chrono = toc;

fprintf(['.. Test sample generated in ', num2str(chrono),' seconds \n']);
%% Construct training set

tic;
data_rom = train_pod(rom, time, base_sample_rom, params, tau); chrono1 = toc;
% data_rom = training{1}; data_rhs = training{2};

tic;
data_prom = train_pod(prom, time, base_sample_prom, params, tau); chrono2 = toc;
% data_prom = training{1}; data_rhs = training{2};

fprintf(['.. Training sets constructed in respectively ', num2str(chrono1), ' and ', num2str(chrono2), ' seconds \n']);
%% Calculate the SVD on the snapshot matrix

[U, S, U_pod, U_old] = svd_pod(data_rom, n_osc, n_latent_pod);
[V, T, V_pod, V_old] = svd_pod(data_prom, n_osc, n_latent_pod);

fprintf('.. SVD done \n');

%% Visualise surface response (1D)

loss_fcn_om=@(orderp,omega) (orderp-1).^2 + 7e1*(2*pi/24-omega).^2;
loss_fcn_L0=@(orderp,L0) (orderp-1).^2 + 1.75e3*L0.^2;

% Oversample to visualize surface response
om_sample=linspace(lb(1),ub(1),20);
L0_sample=linspace(lb(2),ub(2),20);

for j=1:length(om_sample)
    
    params.omega = 2*pi/om_sample(j);
        
    [t_loss,x_loss]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
    [synchrony_loss(:,j),sync_param_loss(j),avg_gene_conc_loss(:,j),spectral_ampl_fact_loss(j),estimated_period_loss(j),order_param_X_loss(:,j)]=...
        compute_QoIs(t_loss,x_loss,n_osc,time.n_asympt,params.omega,params.L0);
    
    %lossfcn(j)=loss_fcn(estimated_period_loss(j),L0vec(j));
    lossfcn_om(j)=loss_fcn_om(order_param_X_loss(end,j),params.omega);
end
params.omega=2*pi/24;
order_om = order_param_X_loss;

for j=1:length(L0_sample)
    
    params.L0 = L0_sample(j);
        
    [t_loss,x_loss]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
    [synchrony_loss(:,j),sync_param_loss(j),avg_gene_conc_loss(:,j),spectral_ampl_fact_loss(j),estimated_period_loss(j),order_param_X_loss(:,j)]=...
        compute_QoIs(t_loss,x_loss,n_osc,time.n_asympt,params.omega,params.L0);
    
    %lossfcn(j)=loss_fcn(estimated_period_loss(j),L0vec(j));
    lossfcn_L0(j)=loss_fcn_L0(order_param_X_loss(end,j),params.L0);
end
params.L0=0.01;

%visualize function
FIG = figure;
t = tiledlayout(1,2);

    nexttile;
    plot(om_sample,loss_fcn_om(1,2*pi./om_sample), 'Color', colours{1}, 'LineWidth', 5); hold on;
    plot(om_sample,loss_fcn_om(order_om(end,:),2*pi./24), 'Color', colours{2}, 'LineWidth', 5);
    plot(om_sample,lossfcn_om, 'Color', colours{5}, 'LineWidth', 5);
    xlabel('$\Omega$','interpreter','latex'); 
    ylabel('$\mathcal{L}_1(\Omega)$','interpreter','latex');

    nexttile;
    ax1 = plot(L0_sample,loss_fcn_L0(1,L0_sample), 'Color', colours{1}, 'LineWidth', 5); hold on;
    ax2 = plot(L0_sample,loss_fcn_L0(order_param_X_loss(end,:),0), 'Color', colours{2}, 'LineWidth', 5);
    ax3 = plot(L0_sample,lossfcn_L0, 'Color', colours{5}, 'LineWidth', 5);
    xlabel('$L_0$','interpreter','latex'); 
    ylabel('$\mathcal{L}_2(L_0)$','interpreter','latex');

FIG.WindowState = 'fullscreen';
lg = legend([ax2,ax1,ax3], 'optimisation', 'regularisation', 'final objective', 'interpreter', 'latex', 'FontSize', 30);
lg.Layout.Tile = 'north';
title(t, '\textbf{Optimisation: objective function}', 'interpreter', 'latex', 'FontSize', 50);
saveas(FIG,strcat('figures/control/opti_loss_', num2str(n_latent_pod)),'epsc')

% % Visualize gradient
% figure;
% plot(L0_sample,1/max(diff(L0_sample))*[lossfcn_om(2)-lossfcn_om(1), 1/2*(lossfcn_om(3:end)-lossfcn_om(1:end-2)), lossfcn_om(end)-lossfcn_om(end-1)]); xlabel('$L_0$','interpreter','latex'); ylabel('$\nabla Loss$','interpreter','latex');
% hold on; 
% yline(0);

%% Control

% Loss function on the order parameter with starting point L0_max
init = [ub(1) ub(2)];
learn_rate=[1e1 8e-6];
epsilon=1e-6;

tol = 1e-4;
Niter=40;

current_fom = init; current_rom = init; current_prom = init;

for k=1:Niter
    disp(k);  
    
    % FOM
    
    tstart=tic;
    
        % Evolve for current value of omega
        params.omega = 2*pi/current_fom(1); params.L0 = 0.01;
        [t_fom,x_fom]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
        [~,~,~,~,~,order_param_X_fom]=...
            compute_QoIs(t_fom,x_fom,n_osc,time.n_asympt,params.omega,params.L0);
        % Evolve for current value of omega + epsilon
        params.omega = 2*pi/(current_fom(1) + epsilon);
        [t_fom_shift,x_fom_shift]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
        [~,~,~,~,~,order_param_X_fom_shift]=...
            compute_QoIs(t_fom_shift,x_fom_shift,n_osc,time.n_asympt,params.omega, params.L0);
    
     time_fom(k,1)=toc(tstart); 
     gradient_fom(1)=1/epsilon*(loss_fcn_om(order_param_X_fom_shift(end),2*pi/(current_fom(1)+epsilon))-loss_fcn_om(order_param_X_fom(end),2*pi/current_fom(1)));
     delta_loss_fom(k,1)=loss_fcn_om(order_param_X_fom(end),2*pi/current_fom(1));
     tic; 
        
        % Evolve for current value of L0
        params.omega = 2*pi/24; params.L0 = current_fom(2);
        [t_fom,x_fom]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
        [~,~,~,~,~,order_param_X_fom]=...
            compute_QoIs(t_fom,x_fom,n_osc,time.n_asympt,params.omega,params.L0);
        % Evolve for current value of L0 + epsilon
        params.L0 = current_fom(2) + epsilon;
        [t_fom_shift,x_fom_shift]=RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 time.final_time], time.x0, time.dt);
        [~,~,~,~,~,order_param_X_fom_shift]=...
            compute_QoIs(t_fom_shift,x_fom_shift,n_osc,time.n_asympt,params.omega,params.L0);
    
    time_fom(k,2)=toc(tstart);
    gradient_fom(2)=1/epsilon*(loss_fcn_L0(order_param_X_fom_shift(end),current_fom(2)+epsilon)-loss_fcn_L0(order_param_X_fom(end),current_fom(2)));
    delta_loss_fom(k,2)=loss_fcn_L0(order_param_X_fom(end),current_fom(2));

    % ROM
    
    tstart=tic;
    
        % Evolve for current value of omega
        params.omega = 2*pi/current_rom(1); params.L0 = 0.01;        
        [t_rom,x_rom_proj]=RK3 (@(t,x) U_pod'*circadian_rhs(t,U_pod*x,params,tau),[0 time.final_time], U_pod'*time.x0, time.dt);
        x_rom=x_rom_proj*U_pod';
        [~,~,~,~,~,order_param_X_rom]=...
            compute_QoIs(t_rom,x_rom,n_osc,time.n_asympt,params.omega,params.L0);
        % Evolve for current value of L0 + epsilon
        params.omega = 2*pi/(current_rom(1) + epsilon);
        [t_rom_shift,x_rom_shift_proj]=RK3 (@(t,x) U_pod'*circadian_rhs(t,U_pod*x,params,tau),[0 time.final_time], U_pod'*time.x0, time.dt);
        x_rom_shift=x_rom_shift_proj*U_pod';
        [~,~,~,~,~,order_param_X_rom_shift]=...
            compute_QoIs(t_rom_shift,x_rom_shift,n_osc,time.n_asympt,params.omega,params.L0);
    
     time_rom(k,1)=toc(tstart); 
     gradient_rom(1)=1/epsilon*(loss_fcn_om(order_param_X_rom_shift(end),2*pi/(current_rom(1)+epsilon))-loss_fcn_om(order_param_X_rom(end),2*pi/current_rom(1)));
     delta_loss_rom(k,1)=loss_fcn_om(order_param_X_rom(end),2*pi/current_rom(1));
     tic; 
        
        % Evolve for current value of L0
        params.omega = 2*pi/24; params.L0 = current_rom(2);
        [t_rom,x_rom_proj]=RK3 (@(t,x) U_pod'*circadian_rhs(t,U_pod*x,params,tau),[0 time.final_time], U_pod'*time.x0, time.dt);
        x_rom=x_rom_proj*U_pod';
        [~,~,~,~,~,order_param_X_rom]=...
            compute_QoIs(t_rom,x_rom,n_osc,time.n_asympt,params.omega,params.L0);
        % Evolve for current value of L0 + epsilon
        params.L0 = current_rom(2) + epsilon;
        [t_rom_shift,x_rom_shift_proj]=RK3 (@(t,x) U_pod'*circadian_rhs(t,U_pod*x,params,tau),[0 time.final_time], U_pod'*time.x0, time.dt);
        x_rom_shift=x_rom_shift_proj*U_pod';
        [~,~,~,~,~,order_param_X_rom_shift]=...
            compute_QoIs(t_rom_shift,x_rom_shift,n_osc,time.n_asympt,params.omega,params.L0);
    
    time_rom(k,2)=toc(tstart);
    gradient_rom(2)=1/epsilon*(loss_fcn_L0(order_param_X_rom_shift(end),current_rom(2)+epsilon)-loss_fcn_L0(order_param_X_rom(end),current_rom(2)));
    delta_loss_rom(k,2)=loss_fcn_L0(order_param_X_rom(end),current_rom(2));
    
    % PROM
    tstart=tic;
    
        % Evolve for current value of omega
        params.omega = 2*pi/current_prom(1); params.L0 = 0.01;        
        [t_prom,x_prom_proj]=RK3 (@(t,x) V_pod'*circadian_rhs(t,V_pod*x,params,tau),[0 time.final_time], V_pod'*time.x0, time.dt);
        x_prom=x_prom_proj*V_pod';
        [~,~,~,~,~,order_param_X_prom]=...
            compute_QoIs(t_prom,x_prom,n_osc,time.n_asympt,params.omega,params.L0);
        % Evolve for current value of omega + epsilon
        params.omega = 2*pi/(current_prom(1) + epsilon); 
        [t_prom_shift,x_prom_proj_shift]=RK3 (@(t,x) V_pod'*circadian_rhs(t,V_pod*x,params,tau),[0 time.final_time], V_pod'*time.x0, time.dt);
        x_prom_shift=x_prom_proj_shift*V_pod';
        [~,~,~,~,~,order_param_X_prom_shift]=...
            compute_QoIs(t_prom_shift,x_prom_shift,n_osc,time.n_asympt,params.omega,params.L0);
    
     time_prom(k,1)=toc(tstart); 
     gradient_prom(1)=1/epsilon*(loss_fcn_om(order_param_X_prom_shift(end),2*pi/(current_prom(1)+epsilon))-loss_fcn_om(order_param_X_prom(end),2*pi/current_prom(1)));
     delta_loss_prom(k,1)=loss_fcn_om(order_param_X_prom(end),2*pi/current_prom(1));
     tic; 
        
        % Evolve for current value of L0
        params.omega = 2*pi/24; params.L0 = current_prom(2);
        [t_prom,x_prom_proj]=RK3 (@(t,x) V_pod'*circadian_rhs(t,V_pod*x,params,tau),[0 time.final_time], V_pod'*time.x0, time.dt);
        x_prom=x_prom_proj*V_pod';
        [~,~,~,~,~,order_param_X_prom]=...
            compute_QoIs(t_prom,x_prom,n_osc,time.n_asympt,params.omega,params.L0);        
        % Evolve for current value of L0 + epsilon
        params.L0 = current_prom(2) + epsilon;
        [t_prom_shift,x_prom_proj_shift]=RK3 (@(t,x) V_pod'*circadian_rhs(t,V_pod*x,params,tau),[0 time.final_time], V_pod'*time.x0, time.dt);
        x_prom_shift=x_prom_proj_shift*V_pod';
        [~,~,~,~,~,order_param_X_prom_shift]=...
            compute_QoIs(t_prom_shift,x_prom_shift,n_osc,time.n_asympt,params.omega,params.L0);
    
    time_prom(k,2)=toc(tstart);
    gradient_prom(2)=1/epsilon*(loss_fcn_L0(order_param_X_prom_shift(end),current_prom(2)+epsilon)-loss_fcn_L0(order_param_X_prom(end),current_prom(2)));
    delta_loss_prom(k,2)=loss_fcn_L0(order_param_X_prom(end),current_prom(2));
    
    iter_fom(k,:)=current_fom;
    iter_rom(k,:)=current_rom;
    iter_prom(k,:)=current_prom;
    
    % Update
    current_fom=current_fom-learn_rate.*gradient_fom;
    current_rom=current_rom-learn_rate.*gradient_rom;
    current_prom=current_prom-learn_rate.*gradient_prom;
    
%     delta = [learn_rate.*gradient_fom, learn_rate.*gradient_rom, learn_rate.*gradient_prom];
%     idx = find(delta<tol);
%     pr = {'FOM: Omega', 'FOM: L0', 'ROM: Omega', 'ROM: L0', 'pROM: Omega', 'pROM: L0'};
%     
%     if ~ isempty(idx)
%         disp(pr{idx})
%     end
    
end



% FIG = figure;
% t = tiledlayout(1,2);
% 
% nexttile;
% plot(1:Niter,time_fom(:,1),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,time_rom(:,1),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,time_prom(:,1),'x','color',colours{3}, 'MarkerSize', 16);
% xlabel('$iter$','interpreter','latex'); ylabel('$Cost$','interpreter','latex');
% 
% nexttile;
% plot(1:Niter,time_fom(:,2),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,time_rom(:,2),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,time_prom(:,2),'x','color',colours{3}, 'MarkerSize', 16);
% xlabel('$iter$','interpreter','latex');
% 
% lg = legend({'FOM','ROM','pROM'},'interpreter','latex', 'FontSize', 30);
% lg.Layout.Tile = 'north';
% title(t, '\textbf{Optimisation: cost}', 'interpreter','latex', 'FontSize', 50)
% 
% FIG.WindowState = 'fullscreen';
% saveas(FIG,strcat('figures/control/opti_cost_', num2str(n_latent_pod)),'epsc')


% FIG = figure;
% t = tiledlayout(2,2);
% 
% nexttile;
% plot(1:Niter,iter_fom(:,1),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,iter_rom(:,1),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,iter_prom(:,1),'x','color',colours{3}, 'MarkerSize', 16);
% ylabel('$\Omega$','interpreter','latex');
% 
% nexttile;
% plot(1:Niter,iter_fom(:,2),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,iter_rom(:,2),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,iter_prom(:,2),'x','color',colours{3}, 'MarkerSize', 16);
% ylabel('$L_0$','interpreter','latex');
% 
% nexttile;
% plot(1:Niter,delta_loss_fom(:,1),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,delta_loss_rom(:,1),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,delta_loss_prom(:,1),'x','color',colours{3}, 'MarkerSize', 16);
% xlabel('$iter$','interpreter','latex'); ylabel('$loss$','interpreter','latex');
% 
% nexttile;
% plot(1:Niter,delta_loss_fom(:,2),'.','color',colours{1}, 'MarkerSize', 40); hold on;
% plot(1:Niter,delta_loss_rom(:,2),'x','color',colours{2}, 'MarkerSize', 16);
% plot(1:Niter,delta_loss_prom(:,2),'x','color',colours{3}, 'MarkerSize', 16);
% xlabel('$iter$','interpreter','latex'); ylabel('$loss$','interpreter','latex');
% 
% lg = legend({'FOM','ROM','pROM'},'interpreter','latex', 'FontSize', 30);
% lg.Layout.Tile = 'north';
% title(t, '\textbf{Optimisation: }', 'interpreter','latex', 'FontSize', 50)
% 
% FIG.WindowState = 'fullscreen';
% saveas(FIG,strcat('figures/control/opti_iter_', num2str(n_latent_pod)),'epsc')

FIG = figure;
t = tiledlayout(1,2);

nexttile;
plot(om_sample,lossfcn_om, 'Color', [0.8 0.8 0.8], 'LineWidth', 4); hold on;
plot(iter_fom(:,1),delta_loss_fom(:,1),'.','color',colours{1}, 'MarkerSize', 40);
plot(iter_rom(:,1),delta_loss_rom(:,1),'x','color',colours{2}, 'MarkerSize', 16);
plot(iter_prom(:,1),delta_loss_prom(:,1),'x','color',colours{3}, 'MarkerSize', 16);
xline(iter_fom(end,1), 'color',colours{1}, 'LineWidth', 4); 
xline(iter_rom(end,1), '--', 'color',colours{2}, 'LineWidth', 4); 
xline(iter_prom(end,1), ':', 'color',colours{3}, 'LineWidth', 4);
xlim([lb(1) ub(1)]);
xlabel('$\Omega$','interpreter','latex'); ylabel('$\mathcal{L}_1(\Omega)$','interpreter','latex');

[~, idx_rom] = min(delta_loss_rom(:,2)); [~, idx_prom] = min(delta_loss_prom(:,2)); 
nexttile;
plot(L0_sample,lossfcn_L0, 'Color', [0.8 0.8 0.8], 'LineWidth', 4); hold on;
ax1 = plot(iter_fom(:,2),delta_loss_fom(:,2),'.','color',colours{1}, 'MarkerSize', 40);
ax2 = plot(iter_rom(:,2),delta_loss_rom(:,2),'x','color',colours{2}, 'MarkerSize', 16);
ax3 = plot(iter_prom(:,2),delta_loss_prom(:,2),'x','color',colours{3}, 'MarkerSize', 16);
ax4 = xline(iter_fom(end,2), 'color',colours{1}, 'LineWidth', 4); 
ax5 = xline(iter_rom(idx_rom,2), '--', 'color',colours{2}, 'LineWidth', 4); 
ax6 = xline(iter_prom(idx_prom,2), ':', 'color',colours{3}, 'LineWidth', 4);
xlim([0.003 ub(2)]);
xlabel('$L_0$','interpreter','latex'); ylabel('$\mathcal{L}_2(L_0)$','interpreter','latex');

lg = legend([ax1 ax2 ax3 ax4 ax5 ax6], {'FOM','ROM','pROM', 'min FOM', 'min ROM', 'min pROM'},'interpreter','latex', 'FontSize', 30, 'NumColumns', 2);
lg.Layout.Tile = 'north';
title(t, '\textbf{Optimisation: loss evolution}', 'interpreter','latex', 'FontSize', 50)

FIG.WindowState = 'fullscreen';
saveas(FIG,strcat('figures/control/opti_deltaloss_', num2str(n_latent_pod)),'epsc')

end

