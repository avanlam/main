function [full_outp, rom_outp, prom_outp, deim_outp] = simulate_sample_grid(x_all, y_all, params, tau, final_time, x0, dt, n_asympt, U, U_pod, V_pod)

m = size(x_all, 2);
n = size(y_all, 2);

num_n = numel(tau);

[tmp1,tmp2]=meshgrid(x_all,y_all);

for j = 1 : numel(tmp1)
    params.nu(6) = tmp1(j);
    params.omega = 2*pi/tmp2(j);
    
    % Full order model (FOM)
    tic;
    [time_full,x_full(:,:,j)] = RK3 (@(t,x) circadian_rhs(t,x,params,tau),[0 final_time], x0, dt );
    elapsed_time_full(j)=toc;
    [synchrony_full(:,j),sync_param_full(j),avg_gene_conc_full(:,j),spectral_ampl_fact_full(j),estimated_period_full(j),order_param_X_full(:,j)] = compute_QoIs(time_full,x_full(:,:,j),num_n,n_asympt,params.omega,params.L0);
    
    % Reduced order model (ROM) : POD
    tic;
    [time_rom,x_rom_proj] = RK3 (@(t,x) U_pod'*circadian_rhs(t,U_pod*x,params,tau),[0 final_time], U_pod'*x0, dt);
    x_rom(:,:,j)=x_rom_proj*U_pod';
    elapsed_time_rom(j)=toc;
    [synchrony_rom(:,j),sync_param_rom(j),avg_gene_conc_rom(:,j),spectral_ampl_fact_rom(j),estimated_period_rom(j),order_param_X_rom(:,j)] = compute_QoIs(time_rom,x_rom(:,:,j),num_n,n_asympt,params.omega,params.L0);
        
    % Parametric reduced order model (pROM) : POD
    tic;
    [time_prom,x_prom_proj] = RK3 (@(t,x) V_pod'*circadian_rhs(t,V_pod*x,params,tau),[0 final_time], V_pod'*x0, dt);
    x_prom(:,:,j)=x_prom_proj*V_pod';
    elapsed_time_prom(j)=toc;
    [synchrony_prom(:,j),sync_param_prom(j),avg_gene_conc_prom(:,j),spectral_ampl_fact_prom(j),estimated_period_prom(j),order_param_X_prom(:,j)] = compute_QoIs(time_prom,x_prom(:,:,j),num_n,n_asympt,params.omega,params.L0);
        
    % Reduced order model (ROM) : DEIM with snapshots
    Pdeim=deim(U{4});
    Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
    PTU=U_pod(Pdeim_ext,:);
    PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
    meanU_v=mean(U{4},1);
    PTtau=tau(Pdeim);
    
    tic;
    [time_deim,x_deim_proj]=RK3(@(t,x) PTUinv*circadian_rhs_colloc(t,x,params,PTtau,PTU,meanU_v), [0 final_time], U_pod'*x0, dt);
    x_deim=x_deim_proj*U_pod';
    elapsed_time_deim(j)=toc;
    [synchrony_deim(:,j),sync_param_deim(j),avg_gene_conc_deim(:,j),spectral_ampl_fact_deim(j),estimated_period_deim(j),order_param_X_deim(:,j)] = compute_QoIs(time_deim,x_deim,num_n,n_asympt,params.omega,params.L0);

end

full_outp.time = time_full;
full_outp.x = x_full;
full_outp.synchrony = permute(reshape(synchrony_full, size(synchrony_full, 1), n, m), [1 3 2]);
full_outp.sync_param = reshape(sync_param_full, n, m)';
full_outp.avg_gene = permute(reshape(avg_gene_conc_full, size(avg_gene_conc_full, 1), n, m), [1 3 2]);
full_outp.spectral = reshape(spectral_ampl_fact_full, n, m)';
full_outp.period = reshape(estimated_period_full, n, m)';
full_outp.order_param = permute(reshape(order_param_X_full, size(order_param_X_full, 1), n, m), [1 3 2]);
full_outp.elapsed_time = reshape(elapsed_time_full, n, m)';

rom_outp.time = time_rom;
rom_outp.x = x_rom;
rom_outp.synchrony = permute(reshape(synchrony_rom, size(synchrony_rom, 1), n, m), [1 3 2]);
rom_outp.sync_param = reshape(sync_param_rom, n, m)';
rom_outp.avg_gene = permute(reshape(avg_gene_conc_rom, size(avg_gene_conc_rom, 1), n, m), [1 3 2]);
rom_outp.spectral = reshape(spectral_ampl_fact_rom, n, m)';
rom_outp.period = reshape(estimated_period_rom, n, m)';
rom_outp.order_param = permute(reshape(order_param_X_rom, size(order_param_X_rom, 1), n, m), [1 3 2]);
rom_outp.elapsed_time = reshape(elapsed_time_rom, n, m)';

prom_outp.time = time_prom;
prom_outp.x = x_prom;
prom_outp.synchrony = permute(reshape(synchrony_prom, size(synchrony_prom, 1), n, m), [1 3 2]);
prom_outp.sync_param = reshape(sync_param_prom, n, m)';
prom_outp.avg_gene = permute(reshape(avg_gene_conc_prom, size(avg_gene_conc_prom, 1), n, m), [1 3 2]);
prom_outp.spectral = reshape(spectral_ampl_fact_prom, n, m)';
prom_outp.period = reshape(estimated_period_prom, n, m)';
prom_outp.order_param = permute(reshape(order_param_X_prom, size(order_param_X_prom, 1), n, m), [1 3 2]);
prom_outp.elapsed_time = reshape(elapsed_time_prom, n, m)';

deim_outp.time = time_deim;
deim_outp.x = x_deim;
deim_outp.synchrony = permute(reshape(synchrony_deim, size(synchrony_deim, 1), n, m), [1 3 2]);
deim_outp.sync_param = reshape(sync_param_deim, n, m)';
deim_outp.avg_gene = permute(reshape(avg_gene_conc_deim, size(avg_gene_conc_deim, 1), n, m), [1 3 2]);
deim_outp.spectral = reshape(spectral_ampl_fact_deim, n, m)';
deim_outp.period = reshape(estimated_period_deim, n, m)';
deim_outp.order_param = permute(reshape(order_param_X_deim, size(order_param_X_deim, 1), n, m), [1 3 2]);
deim_outp.elapsed_time = reshape(elapsed_time_deim, n, m)';

end

