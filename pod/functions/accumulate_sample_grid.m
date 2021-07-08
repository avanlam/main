function [full_outp_test, rom_outp_test, prom_outp_test, deim_outp_test] = accumulate_sample_grid(x_all, y_all, params, tau, time, U, U_pod, V_pod)
%ACCUMULATE_SAMPLE organises the simulation of the modified Goodwin model
% into four different output structures, with the mean and standard
% deviations across different sets of natural frequencies for the
% oscillators described in tau, in the case of the test for a grid of 
% parameter configurations
% INPUT: 
%   *'x_all' and 'y_all' contains the test values for the constitutive 
%   parameters at which the model will be sampled in the
%   simulate_sample_grid function.
%   * the 'params' structure contains the standard or default values for the
%   parameters
%   * the 'time' structure contains the detail of the time integration
%   * the 'U', 'U_pod' and 'V_pod' matrices represent the projection 
%   matrices for respectively the DEIM procedure and the ROM + pROM 
%   POD-procedures

[full_outp_test, rom_outp_test, prom_outp_test] = simulate_sample_grid(x_all, y_all, params, tau(:,1), time.final_time, time.x0, time.dt, time.n_asympt, U, U_pod, V_pod);

for i=2:size(tau,2)

    [full_outp, rom_outp, prom_outp, deim_outp] = simulate_sample_grid(x_all, y_all, params, tau(:,i), time.final_time, time.x0, time.dt, time.n_asympt, U, U_pod, V_pod);

    full_outp_test.sync_param(:,:,i) = full_outp.sync_param;
    full_outp_test.spectral(:,:,i) = full_outp.spectral;
    full_outp_test.period(:,:,i) = full_outp.period;
    full_outp_test.order_param(:,:,:,i) = full_outp.order_param;

    rom_outp_test.sync_param(:,:,i) = rom_outp.sync_param;
    rom_outp_test.spectral(:,:,i) = rom_outp.spectral;
    rom_outp_test.period(:,:,i) = rom_outp.period;
    rom_outp_test.order_param(:,:,:,i) = rom_outp.order_param;

    prom_outp_test.sync_param(:,:,i) = prom_outp.sync_param;
    prom_outp_test.spectral(:,:,i) = prom_outp.spectral;
    prom_outp_test.period(:,:,i) = prom_outp.period;
    prom_outp_test.order_param(:,:,:,i) = prom_outp.order_param;
    
    deim_outp_test.sync_param(:,:,i) = deim_outp.sync_param;
    deim_outp_test.spectral(:,:,i) = deim_outp.spectral;
    deim_outp_test.period(:,:,i) = deim_outp.period;
    deim_outp_test.order_param(:,:,:,i) = deim_outp.order_param;
end

% Compute averages
last = ndims(full_outp_test.sync_param);

full_outp_test.sync_param_mean = mean(full_outp_test.sync_param,last); 
rom_outp_test.sync_param_mean = mean(rom_outp_test.sync_param,last); 
prom_outp_test.sync_param_mean = mean(prom_outp_test.sync_param,last);
deim_outp_test.sync_param_mean = mean(deim_outp_test.sync_param,last);

full_outp_test.spectral_mean = mean(full_outp_test.spectral,last); 
rom_outp_test.spectral_mean = mean(rom_outp_test.spectral,last); 
prom_outp_test.spectral_mean = mean(prom_outp_test.spectral,last);
deim_outp_test.spectral_mean = mean(deim_outp_test.spectral,last);

full_outp_test.period_mean = mean(full_outp_test.period,last); 
rom_outp_test.period_mean = mean(rom_outp_test.period,last); 
prom_outp_test.period_mean = mean(prom_outp_test.period,last);
deim_outp_test.period_mean = mean(deim_outp_test.period,last);

last = ndims(full_outp_test.order_param);
full_outp_test.order_param_mean = mean(full_outp_test.order_param,last); 
rom_outp_test.order_param_mean = mean(rom_outp_test.order_param,last); 
prom_outp_test.order_param_mean = mean(prom_outp_test.order_param,last);
deim_outp_test.order_param_mean = mean(deim_outp_test.order_param,last);

% Compute variances
last = ndims(full_outp_test.sync_param);
full_outp_test.sync_param_std = var(full_outp_test.sync_param,1,last); 
rom_outp_test.sync_param_std = var(rom_outp_test.sync_param,1,last); 
prom_outp_test.sync_param_std = var(prom_outp_test.sync_param,1,last);
deim_outp_test.sync_param_std = var(deim_outp_test.sync_param,1,last);

full_outp_test.spectral_std = var(full_outp_test.spectral,1,last); 
rom_outp_test.spectral_std = var(rom_outp_test.spectral,1,last); 
prom_outp_test.spectral_std = var(prom_outp_test.spectral,1,last);
deim_outp_test.spectral_std = var(deim_outp_test.spectral,1,last);

full_outp_test.period_std = var(full_outp_test.period,1,last); 
rom_outp_test.period_std = var(rom_outp_test.period,1,last); 
prom_outp_test.period_std = var(prom_outp_test.period,1,last);
deim_outp_test.period_std = var(deim_outp_test.period,1,last);

last = ndims(full_outp_test.order_param);
full_outp_test.order_param_std = var(full_outp_test.order_param,1,last); 
rom_outp_test.order_param_std = var(rom_outp_test.order_param,1,last); 
prom_outp_test.order_param_std = var(prom_outp_test.order_param,1,last);
deim_outp_test.order_param_std = var(deim_outp_test.order_param,1,last);

end

