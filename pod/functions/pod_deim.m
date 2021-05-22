% function [ rom_inp, prom_inp ] = pod_deim(method, num_osc)

method = 'determ';

num_osc = 100;
set(0,'defaultAxesFontSize',25);

%% Control Parameters
n_samples = 20;                 	% 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
n_test_samples = 15;                % 2. Size of the test sample: total number of function evaluation will be N * (numDim + 2)
n_trials = 1;                       % 3. Number of realisations for the natural frequencies in the base sample
n_trials_samples = 50;              % 4. Number of realisations for the natural frequencies in the test sample
n_latent_pod=10;                    % 5. Reduced number of oscillators, up to N
n_samples_deim=10;                  % 6. Reduced number of oscillators, up to N
seedNum = 123456789;               	% 7. Seed for random number generator; put [] for autmatic randomization 
funPath = ['../syst_', method]; 	% 8. Folder address: that includes model file, factor space
funFile = 'eval_circ';            	% 9. Model/function file: MATLAB m-file without .m extension
smplMtd = 'LHS';                  	% 10. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers; if blank, default is LHS
%% Plot Flags
samplePlotFlag = 0;
svdPlotFlag = 0;
simulPlotFlag = 0;
testFreqFlag = 0;
testParamFlag = 1;
testFreqParamFlag = 1;
testReducedDimFlag = 1;
testFullDimFlag = 1;
%% System Parameters
final_time=500;             % 1. Final time
n_steps=2*final_time;       % 2. Number of time steps
dt=final_time/n_steps;      % 3. Size of time step
n_asympt=n_steps/5;         % 4. Asymptotic time steps (used to compute QoIs)

x0=repmat([0.0998 0.2468 2.0151 0.0339].', num_osc, 1); 	% 5. Initial point on the stable limit cycle, for each neuron
a=0.8; b=1.2; tau_all=a+(b-a)*rand(num_osc,n_trials);     % 6. Intrinsic periods: uniform distribution

params.nu=[0.7 0.35 NaN 0.35 NaN 0.35 NaN 1]';          % 7. Fixed parameters (Michaelis-Menten + Light)
params.Ki=[1 1 NaN 1 NaN 1 NaN 1]';
params.k=[NaN NaN 0.7 NaN 0.7 NaN 0.35]';
params.nuc=0.4;
params.Kc=1;
params.omega=2*pi/24;
params.K=0.3; params.L0=0.01;
%% Store POD Specifications
rom_inp.N = num_osc; 
rom_inp.seedNum = seedNum;
rom_inp.funFile = funFile;
rom_inp.funPath = funPath;
rom_inp.funFile = funFile;
rom_inp.smplMtd = smplMtd;
rom_inp.n_trials = n_trials;
sys_inp.final_time = final_time;
sys_inp.n_steps = n_steps;
sys_inp.dt = dt;
sys_inp.x0 = x0;
%% Randomization
if isempty( seedNum ) == false
    rand('state',seedNum); 
end
%% Generate the original base sample as a simple point
addpath(funPath);
factors = read_factorSpace('_pod'); % read in the factor space
lb = factors.lb; ub = factors.ub; 

dim_sample = length(lb);  % number of factors
rom_inp.dim_sample = dim_sample;

rom_inp.n_samples = 1;

base_sample_rom = 0.5.* ( ub' - lb') + lb';

% p = sobolset(numDim * 2,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
% p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
% baseSample = p(1:N,:);
%         

x_all=linspace(lb(1),ub(1),n_test_samples);
y_all=linspace(lb(2),ub(2),n_test_samples);

[xx, yy] = meshgrid(x_all, y_all);
test_sample = [ reshape(xx, [n_test_samples^2, 1]), reshape(yy, [n_test_samples^2, 1])];

fprintf('.. Base and test sample generated \n');
%% Generate the parametric base sample from a unit hypercube

prom_inp = rom_inp;
prom_inp.n_samples = n_samples;

base_sample_prom = sym_LHS( n_samples, dim_sample, 'maxmin', 10);

for j = 1 : n_samples
    base_sample_prom(j, :) = base_sample_prom(j, :) .* ( ub' - lb') + lb';
end

fprintf('.. Parametric base sample generated \n');

if samplePlotFlag
   FIG = figure();
   FIG.WindowState = 'fullscreen'; 

   colours = colour_pairs('duo');
   scatter(test_sample(:,1),test_sample(:,2), 90, colours{1}, 'filled', 'Marker', 's'); hold on
   scatter(base_sample_prom(:,1),base_sample_prom(:,2), 90, colours{2}, 'filled', 'Marker', 's'); hold off
   legend("base sample","test sample",'interpreter','latex', 'FontSize', 35)
   xlim([lb(1),ub(1)]); ylim([lb(2),ub(2)]); axis square
   ylabel('$\Omega$','interpreter','latex', 'FontSize', 40);
   xlabel('$\nu_6$','interpreter','latex', 'FontSize', 40); 
   title('Parametric samples','interpreter','latex', 'FontSize', 45);

   saveas(FIG, strcat('figures_', method, '/samples_', num2str(num_osc)), 'epsc');
end
%% Construct training set

training = train_pod(rom_inp, sys_inp, base_sample_rom, params, tau_all);
data_rom = training{1}; data_rhs = training{2};

training = train_pod(prom_inp, sys_inp, base_sample_prom, params, tau_all);
data_prom = training{1}; data_rhs = training{2};

fprintf('.. Training set constructed \n');

%% Calculate the SVD on the snapshot matrix

[U, S, U_pod, U_old] = svd_pod(data_rom, num_osc, n_latent_pod);
[V, T, V_pod, V_old] = svd_pod(data_prom, num_osc, n_latent_pod);

fprintf('.. SVD done \n');
%% Plot the decay of the SVD values

if svdPlotFlag
    
    [FIG1, FIG2, FIG3] = plot_svd(U, S, U_old, V, T, V_old);
    FIG1.WindowState = 'fullscreen'; 
    saveas(FIG1,strcat('figures_', method, '/subangle_', num2str(num_osc)),'epsc')  
    FIG2.WindowState = 'fullscreen'; 
    saveas(FIG2,strcat('figures_', method, '/svd_', num2str(num_osc)),'epsc')   
    FIG3.WindowState = 'fullscreen'; 
    saveas(FIG3,strcat('figures_', method, '/podvectors_', num2str(num_osc)),'epsc')   
        
end

%% A single simulation = one random realisation of the model (to exemplify the POD)
tau = tau_all(:,1);


if simulPlotFlag
    [ full_output, rom_outp, prom_outp ] = simulate_sample(base_sample_rom, params, tau, final_time, x0, dt, n_asympt, U_pod, V_pod);
    
    [FIG1, FIG2] = plot_one_simul(num_osc, n_asympt, full_output, rom_outp, prom_outp);
    
    FIG1.WindowState = 'fullscreen';
    saveas(FIG1,strcat('figures_', method, '/osc_', num2str(num_osc)),'epsc')
    FIG2.WindowState = 'fullscreen';
    saveas(FIG2,strcat('figures_', method, '/singlesim_', num2str(num_osc)),'epsc')
    
    FIG3 = plot_optimal_projection(full_output, U_pod, V_pod);
    
    FIG3.WindowState = 'fullscreen';
    saveas(FIG3,strcat('figures_', method, '/optiproj_', num2str(num_osc)),'epsc')
    
    % Display times
    fprintf('Elapsed time for full model: %f\n', full_output.elapsed_time);
    fprintf('Elapsed time for ROM model: %f\n', rom_outp.elapsed_time);
    fprintf('Elapsed time for pROM model: %f\n', prom_outp.elapsed_time);
    
end

%% Simulation for the full range of natural frequency, with fixed parameters

% %deim with snapshots
% Pdeim=deim(U{4});
% Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
% PTU=U_pod(Pdeim_ext,:);
% PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
% meanU_v=mean(U{4},1);

% Test set with iid samples (independent of training)
if testFreqFlag
    rng('default'); tau_all=a+(b-a)*rand(num_osc, n_trials_samples);
    
    [ full_outp_freq, rom_outp_freq, prom_outp_freq ] = accumulate_sample([params.nu(6) params.omega], params, tau_all, final_time, x0, dt, n_asympt, U_pod, V_pod);

    fprintf('.. Test data accumulated \n');

    FIG = figure(); FIG.WindowState = 'fullscreen';
    t = tiledlayout(3,1);
    t.TileSpacing = 'tight';
    sgtitle('\textbf{Variation over the natural frequencies of the oscillators}', 'FontSize', 40, 'interpret', 'latex')
    
    colours = colour_pairs('spring');
    % Visualize confidence interval for synchrony variable
    nexttile;
    plot_confidence_interval(full_outp_freq.synchrony_mean-1.96*sqrt(full_outp_freq.synchrony_std/size(tau_all,2)),full_outp_freq.synchrony_mean+1.96*sqrt(full_outp_freq.synchrony_std/size(tau_all,2)),full_outp_freq.time,colours{1}); hold on;
    plot_confidence_interval(rom_outp_freq.synchrony_mean-1.96*sqrt(rom_outp_freq.synchrony_std/size(tau_all,2)),rom_outp_freq.synchrony_mean+1.96*sqrt(rom_outp_freq.synchrony_std/size(tau_all,2)),rom_outp_freq.time,colours{2});
    plot_confidence_interval(prom_outp_freq.synchrony_mean-1.96*sqrt(prom_outp_freq.synchrony_std/size(tau_all,2)),prom_outp_freq.synchrony_mean+1.96*sqrt(prom_outp_freq.synchrony_std/size(tau_all,2)),prom_outp_freq.time,colours{3}); hold off;
    xlabel('$t$','interpreter','latex'); ylabel('$Q$','interpreter','latex');
    title('Synchrony parameter')

    % Visualize confidence interval for avg gene concentration X(t)
    nexttile;
    plot_confidence_interval(full_outp_freq.avg_gene_mean-1.96*sqrt(full_outp_freq.avg_gene_std/size(tau_all,2)),full_outp_freq.avg_gene_mean+1.96*sqrt(full_outp_freq.avg_gene_std/size(tau_all,2)),full_outp_freq.time,colours{1}); hold on;
    plot_confidence_interval(rom_outp_freq.avg_gene_mean-1.96*sqrt(rom_outp_freq.avg_gene_std/size(tau_all,2)),rom_outp_freq.avg_gene_mean+1.96*sqrt(rom_outp_freq.avg_gene_std/size(tau_all,2)),rom_outp_freq.time,colours{2});
    plot_confidence_interval(prom_outp_freq.avg_gene_mean-1.96*sqrt(prom_outp_freq.avg_gene_std/size(tau_all,2)),prom_outp_freq.avg_gene_mean+1.96*sqrt(prom_outp_freq.avg_gene_std/size(tau_all,2)),prom_outp_freq.time,colours{3}); hold off;
    xlabel('$t$','interpreter','latex'); ylabel('$X$','interpreter','latex');
    title('Average gene concentration')

    % Visualize confidence interval for order parameter as fcn of time
    nexttile;
    plot_confidence_interval(full_outp_freq.order_param_mean-1.96*sqrt(full_outp_freq.order_param_std/size(tau_all,2)),full_outp_freq.order_param_mean+1.96*sqrt(full_outp_freq.order_param_std/size(tau_all,2)),full_outp_freq.time(end-n_asympt:end),colours{1}); hold on;
    plot_confidence_interval(rom_outp_freq.order_param_mean-1.96*sqrt(rom_outp_freq.order_param_std/size(tau_all,2)),rom_outp_freq.order_param_mean+1.96*sqrt(rom_outp_freq.order_param_std/size(tau_all,2)),rom_outp_freq.time(end-n_asympt:end),colours{2});
    plot_confidence_interval(prom_outp_freq.order_param_mean-1.96*sqrt(prom_outp_freq.order_param_std/size(tau_all,2)),prom_outp_freq.order_param_mean+1.96*sqrt(prom_outp_freq.order_param_std/size(tau_all,2)),prom_outp_freq.time(end-n_asympt:end),colours{3}); hold off;
    ylim([0 1]);
    xlabel('$t$','interpreter','latex'); ylabel('$R$','interpreter','latex');
    title('Kuramoto order parameter')

    lg = legend({'FOM','ROM','pROM'},'interpreter','latex');
    lg.Layout.Tile = 'north';
    
    saveas(FIG,strcat('figures_', method, '/testF_', num2str(num_osc)),'epsc')   

end

%% A single simulation for varying v_6 and omega

if testParamFlag
    
    tau=tau_all(:,1);
    colours = colour_pairs('spring');

%     [  full_outp_test, rom_outp_test, prom_outp_test, deim_outp_test ] = simulate_sample_grid(x_all, y_all, params, tau, final_time, x0, dt, n_asympt, U, U_pod, V_pod);
% 
%     [  full_outp_x, rom_outp_x, prom_outp_x, deim_outp_x ] = simulate_sample_grid(x_all, base_sample_rom(2), params, tau, final_time, x0, dt, n_asympt, U, U_pod, V_pod);
% 
%     [  full_outp_y, rom_outp_y, prom_outp_y, deim_outp_y ] = simulate_sample_grid(base_sample_rom(1), y_all, params, tau, final_time, x0, dt, n_asympt, U, U_pod, V_pod);    
%     
%     filename = 'test_parameters.mat'; save(filename, 'full_outp_test', 'rom_outp_test', 'prom_outp_test' ...
%         , 'full_outp_x', 'rom_outp_x', 'prom_outp_x', 'full_outp_y', 'rom_outp_y', 'prom_outp_y');

    % Visualize quantities
        % Rho
        sync_param = [full_outp_test.sync_param(:), rom_outp_test.sync_param(:), prom_outp_test.sync_param(:)];
        FIG = plot_grid(sync_param, '\rho', x_all,y_all, [0.9 1], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_sync_', num2str(num_osc)),'epsc')   

        % Spectral amplification factor
        spectral = [full_outp_test.spectral(:), rom_outp_test.spectral(:), prom_outp_test.spectral(:)];
        FIG = plot_grid(spectral, 'S', x_all,y_all, [0 60], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_spectral_', num2str(num_osc)),'epsc')   

        % Estimated period
        period = [full_outp_test.period(:), rom_outp_test.period(:), prom_outp_test.period(:)];
        FIG = plot_grid(period, '\bar{T}', x_all,y_all, [20 28], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_period_', num2str(num_osc)),'epsc')   

        % Order parameter
        order_param = [full_outp_test.order_param(end,:)', rom_outp_test.order_param(end,:)', prom_outp_test.order_param(end,:)'];
        FIG = plot_grid(order_param, 'R', x_all,y_all, [0 1], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_order_', num2str(num_osc)),'epsc')   

        % Error
        err_rom= 1/sqrt(numel(full_outp_test.x)).*reshape(vecnorm(vecnorm(full_outp_test.x-rom_outp_test.x,2, 1), 2, 2), [n_test_samples^2, 1]);
        err_prom= 1/sqrt(numel(full_outp_test.x)).*reshape(vecnorm(vecnorm(full_outp_test.x-prom_outp_test.x,2, 1), 2, 2), [n_test_samples^2, 1]);
        FIG = plot_grid([err_rom err_prom], 'error', x_all, y_all, [0 0.1], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_error_', num2str(num_osc)),'epsc')   

        % 2-D Error
        xx = {full_outp_x.x, rom_outp_x.x, prom_outp_x.x};
        yy = {full_outp_y.x, rom_outp_y.x, prom_outp_y.x};
        FIG = plot_error(x_all, y_all, xx, yy, lb, ub);
        saveas(FIG,strcat('figures_', method, '/testP_error2d_', num2str(num_osc)),'epsc')   

end


%% simulate for all instances of the random frequency AND varying coupling 

% %deim with snapshots
% Pdeim=deim(U{4});
% Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
% PTU=U_pod(Pdeim_ext,:);
% PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
% meanU_v=mean(U{4},1);

rng('default'); tau_all=a+(b-a)*rand(num_osc,n_trials_samples);

if testFreqParamFlag
    
    [ full_outp_test, rom_outp_test, prom_outp_test, deim_outp_test ] = accumulate_sample_grid(x_all, y_all, params, tau_all, final_time, x0, dt, n_asympt, U, U_pod, V_pod);

    % Visualize quantities
        % Rho
        sync_param = [full_outp_test.sync_param_mean(:), rom_outp_test.sync_param_mean(:), prom_outp_test.sync_param_mean(:), deim_outp_test.sync_param_mean(:)];
        FIG = plot_grid(sync_param, '\rho', x_all,y_all, [0.9 1], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testParam_sync_', num2str(num_osc)),'epsc')   

        % Spectral amplification factor
        spectral = [full_outp_test.spectral_mean(:), rom_outp_test.spectral_mean(:), prom_outp_test.spectral_mean(:), deim_outp_test.spectral_mean(:)];
        FIG = plot_grid(spectral, 'S', x_all,y_all, [0 20], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testParam_spectral_', num2str(num_osc)),'epsc')   

        % Estimated period
        period = [full_outp_test.period_mean(:), rom_outp_test.period_mean(:), prom_outp_test.period_mean(:), deim_outp_test.period_mean(:)];
        FIG = plot_grid(period, '\bar{T}', x_all,y_all, [20 30], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testParam_period_', num2str(num_osc)),'epsc')   

        % Order parameter
        order_param = [full_outp_test.order_param_mean(end,:)', rom_outp_test.order_param_mean(end,:)', prom_outp_test.order_param_mean(end,:)', deim_outp_test.order_param_mean(end,:)'];
        FIG = plot_grid(order_param, 'R', x_all,y_all, [0 1], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testParam_order_', num2str(num_osc)),'epsc')   

        % Cost
        cost = [full_outp_test.elapsed_time(:), rom_outp_test.elapsed_time(:), prom_outp_test.elapsed_time(:), deim_outp_test.elapsed_time(:)];
        FIG = plot_grid(cost, 'cost', x_all,y_all, [0 max(cost,[],'all')], lb, ub);
        saveas(FIG,strcat('figures_', method, '/testParam_cost_', num2str(num_osc)),'epsc')   

end

%% behavior varying reduced dimension using single instances
latent_dim=unique(round(logspace(0,log10(num_osc),20)));
 tau=tau_all(:,1);

if testReducedDimFlag
    colours = colour_pairs('spring');

    for k=1:length(latent_dim)
        disp(k);

        [U, S, U_pod] = svd_pod(data_rom, num_osc, latent_dim(k));
        [V, T, V_pod] = svd_pod(data_prom, num_osc, latent_dim(k));

        Pdeim=deim(U{4});
        Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
        PTU=U_pod(Pdeim_ext,:);
        PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
        meanU_v=mean(U{4},1);
        PTtau=tau(Pdeim);

        [full_output, rom_output, prom_output] = simulate_sample([params.nu(6) params.omega], params, tau, final_time, x0, dt, n_asympt, U_pod, V_pod);

        full_output_k.synchrony(:,k) = full_output.synchrony;
        full_output_k.avg_gene(:,k) = full_output.avg_gene;
        full_output_k.sync_param(:,k) = full_output.sync_param;
        full_output_k.spectral(:,k) = full_output.spectral;
        full_output_k.period(:,k) = full_output.period;
        full_output_k.order_param(:,k) = full_output.order_param;
        full_output_k.elapsed_time(:,k) = full_output.elapsed_time;

        rom_output_k.synchrony(:,k) = rom_output.synchrony;
        rom_output_k.avg_gene(:,k) = rom_output.avg_gene;
        rom_output_k.sync_param(:,k) = rom_output.sync_param;
        rom_output_k.spectral(:,k) = rom_output.spectral;
        rom_output_k.period(:,k) = rom_output.period;
        rom_output_k.order_param(:,k) = rom_output.order_param;
        rom_output_k.elapsed_time(:,k) = rom_output.elapsed_time;

        prom_output_k.synchrony(:,k) = prom_output.synchrony;
        prom_output_k.avg_gene(:,k) = prom_output.avg_gene;
        prom_output_k.sync_param(:,k) = prom_output.sync_param;
        prom_output_k.spectral(:,k) = prom_output.spectral;
        prom_output_k.period(:,k) = prom_output.period;
        prom_output_k.order_param(:,k) = prom_output.order_param;
        prom_output_k.elapsed_time(:,k) = prom_output.elapsed_time;

        %for i=1:4
        %    [Urhs_loc{i},S{i},~]=svd(data_rhs(:,i:4:end-4+i).','econ'); Urhs_loc{i}=Urhs_loc{i}(:,1:latent_dim(k)); Srhs{i}=diag(S{i});
        %end
        %clear Urhs;
        %clear PTUinv PTU;
        %Urhs(ids,:)=(blkdiag(Urhs_loc{:}));
        %Pdeim = deim(Urhs_loc{4});
        %ids_pod=[1:4:4*n_latent_pod-3, 2:4:4*n_latent_pod-2, 3:4:4*n_latent_pod-1, 4:4:4*n_latent_pod];
        %ids_deim=[1:4:4*latent_dim(k)-3, 2:4:4*latent_dim(k)-2, 3:4:4*latent_dim(k)-1, 4:4:4*latent_dim(k)];
        %Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
        %for i=1:4, PTUinv_loc{i}=(U{i}'*Urhs_loc{i})*(Urhs_loc{i}(Pdeim,:)\eye(size(Urhs_loc{i},2))); end;  PTUinv(:,ids_deim)=blkdiag(PTUinv_loc{:});
        %for i=1:4, PTU_loc{i}=U{i}(Pdeim,:); end; PTU(ids_deim,:)=blkdiag(PTU_loc{:});
        %meanU_v=mean(U{4},1);
        %PTtau=tau(Pdeim);

    end

    FIG = figure;
    t = tiledlayout(2,1);
    t.TileSpacing = 'tight';

        nexttile;
        semilogx(latent_dim,full_output_k.order_param(end,:),'color',colours{1}, 'Marker', 'x'); hold on;
        semilogx(latent_dim,rom_output_k.order_param(end,:),'color',colours{2}, 'Marker', 'x'); 
        semilogx(latent_dim,rom_output_k.order_param(end,:),'color',colours{3}, 'Marker', 'x'); 
        ylim([0 1]);
        ylabel('$R$','interpreter','latex');
        title('Order parameter','interpreter','latex')

        lg = legend({'FOM','ROM','pROM'},'interpreter','latex');
        lg.Layout.Tile = 'north';

        nexttile;
        loglog(latent_dim,abs(rom_output_k.order_param(end,:)-full_output_k.order_param(end,:))./abs(full_output_k.order_param(end,:))+1e-16,'color',colours{2}); hold on;
        loglog(latent_dim,abs(prom_output_k.order_param(end,:)-full_output_k.order_param(end,:))./abs(full_output_k.order_param(end,:))+1e-16,'color',colours{3}); 
        xticks(1:10:100); ylim([0 1]);
        xlabel('$k/4$','interpreter','latex'); ylabel('$|R_{FOM}-R_{ROM}|/|R_{FOM}|$','interpreter','latex');
        title('Error on order parameter','interpreter','latex')

    title(t, 'Order parameter while varying reduced dimension', 'interpreter','latex', 'FontSize', 40)
    FIG.WindowState = 'fullscreen';
    saveas(FIG,strcat('figures_', method, '/testK_order_', num2str(num_osc)),'epsc')   


    FIG = figure;
    t = tiledlayout(1,2);
    t.TileSpacing = 'tight';

        nexttile;
        loglog(latent_dim,full_output_k.elapsed_time,'color',colours{1}, 'Marker', 'x'); hold on;
        loglog(latent_dim,rom_output_k.elapsed_time,'color',colours{2}, 'Marker', 'x'); 
        loglog(latent_dim,prom_output_k.elapsed_time,'color',colours{3}, 'Marker', 'x'); 
        xlabel('$k/4$','interpreter','latex'); xticks(1:10:100); 
        ylabel('$Cost$','interpreter','latex');
        title('Cost','interpreter','latex')

        lg = legend({'FOM','ROM','pROM'},'interpreter','latex');
        lg.Layout.Tile = 'north';

        nexttile;
        loglog(mean(full_output_k.elapsed_time),1e-16,'o','color',colours{1}); hold on;
        loglog(rom_output_k.elapsed_time,abs(rom_output_k.order_param(end,:)-full_output_k.order_param(end,:))./abs(full_output_k.order_param(end,:))+1e-16,'o--','color',colours{2}); 
        loglog(prom_output_k.elapsed_time,abs(prom_output_k.order_param(end,:)-full_output_k.order_param(end,:))./abs(full_output_k.order_param(end,:))+1e-16,'o--','color',colours{3}); 
        xlabel('Cost','interpreter','latex'); ylabel('$Err$','interpreter','latex');
        title('Error for each cost','interpreter','latex')

    title(t, 'Cost while varying reduced dimension', 'interpreter','latex', 'FontSize', 40)
    FIG.WindowState = 'fullscreen';
    saveas(FIG,strcat('figures_', method, '/testK_cost_', num2str(num_osc)),'epsc')   

end
%% behavior varying dimension of the full problem

if testFullDimFlag
    colours = colour_pairs('spring');

    full_dim=[20; 50; 100; 200];
    reduced_dim=10*ones(size(full_dim));
    rom_inp_n = rom_inp;
    prom_inp_n = prom_inp;

    for n=1:length(full_dim)
        disp(n);

        rom_inp_n.N = full_dim(n); 
        prom_inp_n.N = full_dim(n); 

        %construct training set and svd
        x0=repmat([0.0998 0.2468 2.0151 0.0339].',full_dim(n),1); sys_inp.x0 = x0;
        tau_all=a+(b-a)*rand(full_dim(n),n_trials); 
        %tau_all=sort(tau_all,1);

        training = train_pod(rom_inp_n, sys_inp, base_sample_rom, params, tau_all);
        data_rom = training{1}; data_rhs = training{2};    

        training = train_pod(prom_inp_n, sys_inp, base_sample_prom, params, tau_all);
        data_prom = training{1}; data_rhs = training{2};

        [U, S, U_pod] = svd_pod(data_rom, full_dim(n), n_latent_pod);
        [V, T, V_pod] = svd_pod(data_prom, full_dim(n), n_latent_pod);

        %simulate models for one instance of the parameters
        params.nu(6)=x_all(end);
        params.omega=y_all(end);
        tau=tau_all(:,1);

    %     Pdeim=deim(U{4});
    %     Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
    %     PTU=U_pod(Pdeim_ext,:);
    %     PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
    %     meanU_v=mean(U{4},1);
    %     PTtau=tau(Pdeim);

        [full_output, rom_output, prom_output] = simulate_sample([params.nu(6) params.omega], params, tau, final_time, x0, dt, n_asympt, U_pod, V_pod);

        full_output_n.synchrony(:,n) = full_output.synchrony;
        full_output_n.avg_gene(:,n) = full_output.avg_gene;
        full_output_n.sync_param(:,n) = full_output.sync_param;
        full_output_n.spectral(:,n) = full_output.spectral;
        full_output_n.period(:,n) = full_output.period;
        full_output_n.order_param(:,n) = full_output.order_param;
        full_output_n.elapsed_time(:,n) = full_output.elapsed_time;

        rom_output_n.synchrony(:,n) = rom_output.synchrony;
        rom_output_n.avg_gene(:,n) = rom_output.avg_gene;
        rom_output_n.sync_param(:,n) = rom_output.sync_param;
        rom_output_n.spectral(:,n) = rom_output.spectral;
        rom_output_n.period(:,n) = rom_output.period;
        rom_output_n.order_param(:,n) = rom_output.order_param;
        rom_output_n.elapsed_time(:,n) = rom_output.elapsed_time;

        prom_output_n.synchrony(:,n) = prom_output.synchrony;
        prom_output_n.avg_gene(:,n) = prom_output.avg_gene;
        prom_output_n.sync_param(:,n) = prom_output.sync_param;
        prom_output_n.spectral(:,n) = prom_output.spectral;
        prom_output_n.period(:,n) = prom_output.period;
        prom_output_n.order_param(:,n) = prom_output.order_param;
        prom_output_n.elapsed_time(:,n) = prom_output.elapsed_time;

        %for i=1:4
        %    [Urhs_loc{i},S{i},~]=svd(data_rhs(:,i:4:end-4+i).','econ'); Urhs_loc{i}=Urhs_loc{i}(:,1:latent_dim(k)); Srhs{i}=diag(S{i});
        %end
        %clear Urhs;
        %clear PTUinv PTU;
        %Urhs(ids,:)=(blkdiag(Urhs_loc{:}));
        %Pdeim = deim(Urhs_loc{4});
        %ids_pod=[1:4:4*n_latent_pod-3, 2:4:4*n_latent_pod-2, 3:4:4*n_latent_pod-1, 4:4:4*n_latent_pod];
        %ids_deim=[1:4:4*latent_dim(k)-3, 2:4:4*latent_dim(k)-2, 3:4:4*latent_dim(k)-1, 4:4:4*latent_dim(k)];
        %Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
        %for i=1:4, PTUinv_loc{i}=(U{i}'*Urhs_loc{i})*(Urhs_loc{i}(Pdeim,:)\eye(size(Urhs_loc{i},2))); end;  PTUinv(:,ids_deim)=blkdiag(PTUinv_loc{:});
        %for i=1:4, PTU_loc{i}=U{i}(Pdeim,:); end; PTU(ids_deim,:)=blkdiag(PTU_loc{:});
        %meanU_v=mean(U{4},1);
        %PTtau=tau(Pdeim);

    end

    FIG = figure;
    t = tiledlayout(2,1);
    t.TileSpacing = 'tight';

        nexttile;
        semilogx(full_dim,full_output_n.order_param(end,:),'color',colours{1}); hold on;
        semilogx(full_dim,rom_output_n.order_param(end,:),'color',colours{2}); 
        semilogx(full_dim,prom_output_n.order_param(end,:),'color',colours{3}); 
        ylim([0 1]);
        ylabel('$R$','interpreter','latex'); xlim([min(full_dim) max(full_dim)]);
        nexttile;
        loglog(full_dim,abs(rom_output_n.order_param(end,:)-full_output_n.order_param(end,:))./abs(full_output_n.order_param(end,:))+1e-16,'color',colours{2}); hold on;
        loglog(full_dim,abs(prom_output_n.order_param(end,:)-full_output_n.order_param(end,:))./abs(full_output_n.order_param(end,:))+1e-16,'color',colours{3}); 
        ylim([0 1]);
        xlabel('$n/4$','interpreter','latex'); ylabel('$|R_{FOM}-R_{ROM}|/|R_{FOM}|$','interpreter','latex'); xlim([min(full_dim) max(full_dim)]);

    title(t, 'Order parameter while varying full dimension', 'interpreter','latex', 'FontSize', 40)
    FIG.WindowState = 'fullscreen';
    saveas(FIG,strcat('figures_', method, '/testK_order_', num2str(num_osc)),'epsc')   

    FIG = figure;
    t = tiledlayout(1,2);
    t.TileSpacing = 'tight';

        nexttile;
        loglog(full_dim,full_output_n.elapsed_time,'color',colours{1}); hold on;
        loglog(full_dim,rom_output_n.elapsed_time,'color',colours{2}); 
        loglog(full_dim,prom_output_n.elapsed_time,'color',colours{3}); 
        xlabel('$n/4$','interpreter','latex'); ylabel('$Cost$','interpreter','latex'); xlim([min(full_dim) max(full_dim)]);

        nexttile;
        loglog(full_output_n.elapsed_time,1e-16*ones(size(full_output_n.elapsed_time)),'o--','color',colours{1}); hold on;
        loglog(rom_output_n.elapsed_time,abs(rom_output_n.order_param(end,:)-full_output_n.order_param(end,:))./abs(full_output_n.order_param(end,:))+1e-16,'o--','color',colours{2}); 
        loglog(prom_output_n.elapsed_time,abs(prom_output_n.order_param(end,:)-full_output_n.order_param(end,:))./abs(full_output_n.order_param(end,:))+1e-16,'o--','color',colours{3}); 
        xlabel('$Cost$','interpreter','latex'); ylabel('$Err$','interpreter','latex');

        lg = legend({'FOM','ROM','pROM'},'interpreter','latex');
        lg.Layout.Tile = 'north';

    title(t, 'Cost while varying full dimension', 'interpreter','latex', 'FontSize', 40)
    FIG.WindowState = 'fullscreen';
    saveas(FIG,strcat('figures_', method, '/testN_cost_', num2str(num_osc)),'epsc')   
end

% end