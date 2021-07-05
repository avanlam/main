function [full_outp_test, rom_outp_test, prom_outp_test] = accumulate_sample(sample, params, tau, time, U, V)
%ACCUMULATE_SAMPLE Summary of this function goes here
%   Detailed explanation goes here


[full_outp_test, rom_outp_test, prom_outp_test] = simulate_sample(sample, params, tau(:,1), time, U, V);

for i=2:size(tau,2)

    [full_outp, rom_outp, prom_outp] = simulate_sample(sample, params, tau(:,i), time, U, V);

    full_outp_test.synchrony(:,i) = full_outp.synchrony;
    full_outp_test.avg_gene(:,i) = full_outp.avg_gene;
    full_outp_test.sync_param(:,i) = full_outp.sync_param;
    full_outp_test.spectral(:,i) = full_outp.spectral;
    full_outp_test.period(:,i) = full_outp.period;
    full_outp_test.order_param(:,i) = full_outp.order_param;

    rom_outp_test.synchrony(:,i) = rom_outp.synchrony;
    rom_outp_test.avg_gene(:,i) = rom_outp.avg_gene;
    rom_outp_test.sync_param(:,i) = rom_outp.sync_param;
    rom_outp_test.spectral(:,i) = rom_outp.spectral;
    rom_outp_test.period(:,i) = rom_outp.period;
    rom_outp_test.order_param(:,i) = rom_outp.order_param;

    prom_outp_test.synchrony(:,i) = prom_outp.synchrony;
    prom_outp_test.avg_gene(:,i) = prom_outp.avg_gene;
    prom_outp_test.sync_param(:,i) = prom_outp.sync_param;
    prom_outp_test.spectral(:,i) = prom_outp.spectral;
    prom_outp_test.period(:,i) = prom_outp.period;
    prom_outp_test.order_param(:,i) = prom_outp.order_param;
end

% Compute averages
last = ndims(full_outp_test.synchrony);
full_outp_test.synchrony_mean = mean(full_outp_test.synchrony,last); 
rom_outp_test.synchrony_mean = mean(rom_outp_test.synchrony,last); 
prom_outp_test.synchrony_mean = mean(full_outp_test.synchrony,last);

full_outp_test.avg_gene_mean = mean(full_outp_test.avg_gene,last); 
rom_outp_test.avg_gene_mean = mean(rom_outp_test.avg_gene,last); 
prom_outp_test.avg_gene_mean = mean(prom_outp_test.avg_gene,last);

full_outp_test.sync_param_mean = mean(full_outp_test.sync_param,last); 
rom_outp_test.sync_param_mean = mean(rom_outp_test.sync_param,last); 
prom_outp_test.sync_param_mean = mean(prom_outp_test.sync_param,last);

full_outp_test.spectral_mean = mean(full_outp_test.spectral,last); 
rom_outp_test.spectral_mean = mean(rom_outp_test.spectral,last); 
prom_outp_test.spectral_mean = mean(prom_outp_test.spectral,last);

full_outp_test.period_mean = mean(full_outp_test.period,last); 
rom_outp_test.period_mean = mean(rom_outp_test.period,last); 
prom_outp_test.period_mean = mean(prom_outp_test.period,last);

full_outp_test.order_param_mean = mean(full_outp_test.order_param,last); 
rom_outp_test.order_param_mean = mean(rom_outp_test.order_param,last); 
prom_outp_test.order_param_mean = mean(prom_outp_test.order_param,last);

% Compute variances
full_outp_test.synchrony_std = var(full_outp_test.synchrony,1,last); 
rom_outp_test.synchrony_std = var(rom_outp_test.synchrony,1,last);
prom_outp_test.synchrony_std = var(prom_outp_test.synchrony,1,last);

full_outp_test.avg_gene_std = var(full_outp_test.avg_gene,1,last); 
rom_outp_test.avg_gene_std = var(rom_outp_test.avg_gene,1,last); 
prom_outp_test.avg_gene_std = var(prom_outp_test.avg_gene,1,last);

full_outp_test.sync_param_std = var(full_outp_test.sync_param,1,last); 
rom_outp_test.sync_param_std = var(rom_outp_test.sync_param,1,last); 
prom_outp_test.sync_param_std = var(prom_outp_test.sync_param,1,last);

full_outp_test.spectral_std = var(full_outp_test.spectral,1,last); 
rom_outp_test.spectral_std = var(rom_outp_test.spectral,1,last); 
prom_outp_test.spectral_std = var(prom_outp_test.spectral,1,last);

full_outp_test.period_std = var(full_outp_test.period,1,last); 
rom_outp_test.period_std = var(rom_outp_test.period,1,last); 
prom_outp_test.period_std = var(prom_outp_test.period,1,last);

full_outp_test.order_param_std = var(full_outp_test.order_param,1,last); 
rom_outp_test.order_param_std = var(rom_outp_test.order_param,1,last); 
prom_outp_test.order_param_std = var(prom_outp_test.order_param,1,last);

end

