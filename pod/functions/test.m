n_osc = 100;

%% Control Parameters
n_samples = 20;                 	% 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
n_test_samples = 15;                % 2. Size of the test sample: total number of function evaluation will be N * (numDim + 2)
n_trials = 1;                       % 3. Number of realisations for the natural frequencies in the base sample
n_trials_samples = 50;              % 4. Number of realisations for the natural frequencies in the test sample
n_latent_pod=10;                    % 5. Reduced number of oscillators, up to N
n_samples_deim=10;                  % 6. Reduced number of oscillators, up to N
seedNum = 123456789;               	% 7. Seed for random number generator; put [] for autmatic randomization 
funPath = '../syst';                % 8. Folder address: that includes model file, factor space
funFile = 'eval_circ';            	% 9. Model/function file: MATLAB m-file without .m extension
smplMtd = 'LHS';                  	% 10. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers; if blank, default is LHS
%% System Parameters
final_time=500;             % 1. Final time
n_steps=final_time;      % 2. Number of time steps
dt=final_time/n_steps;      % 3. Size of time step
n_asympt=n_steps/5;         % 4. Asymptotic time steps (used to compute QoIs)

x0=repmat([0.0998 0.2468 2.0151 0.0339].', n_osc, 1); 	% 5. Initial point on the stable limit cycle, for each neuron
a=0.8; b=1.2; tau_all=a+(b-a)*rand(n_osc,n_trials);     % 6. Intrinsic periods: uniform distribution

params.nu=[0.7 0.35 NaN 0.35 NaN 0.35 NaN 1]';          % 7. Fixed parameters (Michaelis-Menten + Light)
params.Ki=[1 1 NaN 1 NaN 1 NaN 1]';
params.k=[NaN NaN 0.7 NaN 0.7 NaN 0.35]';
params.nuc=0.4;
params.Kc=1;
params.omega=2*pi/24;
params.K=0.3; params.L0=0.01;
%% Store POD Specifications
rom_inp.n_osc = n_osc; 
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

%% 
addpath(funPath);
factors = read_factorSpace('_pod'); % read in the factor space
lb = factors.lb; ub = factors.ub; 

dim_sample = length(lb);  % number of factors
rom_inp.dim_sample = dim_sample;

prom_inp = rom_inp;
prom_inp.n_samples = 16;

x_all=linspace(lb(1),ub(1),4);
y_all=linspace(lb(2),ub(2),4);
[xx, yy] = meshgrid(x_all, y_all);
base_first = [ reshape(xx, [16, 1]), reshape(yy, [16, 1])];


addpath(funPath);
factors = read_factorSpace('_pod_2'); % read in the factor space
lb = factors.lb; ub = factors.ub; 

dim_sample = length(lb);  % number of factors
rom_inp.dim_sample = dim_sample;

prom_inp = rom_inp;
prom_inp.n_samples = 16;


x_all=linspace(lb(1),ub(1),4);
y_all=linspace(lb(2),ub(2),4);
[xx, yy] = meshgrid(x_all, y_all);
base_second = [ reshape(yy, [16, 1]), reshape(xx, [16, 1])];

fprintf('.. Parametric base sample generated \n');

disp(size(base_first))
disp(size(intersect(base_first, base_second, 'rows')))


%% 


training = train_pod(prom_inp, sys_inp, base_first, params, tau_all);
data_prom_1 = training{1}; data_rhs = training{2};

training = train_pod(prom_inp, sys_inp, base_second, params, tau_all);
data_prom_2 = training{1}; data_rhs = training{2};

fprintf('.. Training generated \n');
disp(size(data_prom_1))
disp(size(intersect(data_prom_1, data_prom_2, 'rows')))

size(setdiff(data_prom_1, data_prom_2, 'rows'))

%% 

[U, S, U_pod, U_old] = svd_pod(data_prom_1, n_osc, n_latent_pod);
[V, T, V_pod, V_old] = svd_pod(data_prom_2, n_osc, n_latent_pod);


