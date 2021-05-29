function [xx, yy, zz] = build_and_test(base_sample, prom, time, tau_all, n_osc, n_latent_pod, test_samples, params, tau, U, U_pod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    data = train_pod(prom, time, base_sample, params, tau_all); 
    [~, ~, W_pod, ~] = svd_pod(data, n_osc, n_latent_pod);
    
    n_test_samples = size(test_samples{1},1);
    
    [  full_x, rom_x, prom_x, ~ ] = simulate_sample_grid(test_samples{1}, [n_test_samples 1 1], params, tau, time, U, U_pod, W_pod);
    [  full_y, rom_y, prom_y, ~ ] = simulate_sample_grid(test_samples{2}, [1 n_test_samples 1], params, tau, time, U, U_pod, W_pod);    
    [  full_z, rom_z, prom_z, ~ ] = simulate_sample_grid(test_samples{3}, [1 1 n_test_samples], params, tau, time, U, U_pod, W_pod);    

    xx = {full_x.x, rom_x.x, prom_x.x};
    yy = {full_y.x, rom_y.x, prom_y.x};
    zz = {full_z.x, rom_z.x, prom_z.x};
end

