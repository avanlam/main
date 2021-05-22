function [full_outp, ROM_outp, pROM_outp] = simulate_sample(sample, params, tau, final_time, x0, dt, n_asympt, U_pod, V_pod)

num_n = numel(tau);

for j = 1 : size(sample, 1)
    params.nu(6) = sample(j, 1);
    params.omega = 2*pi/sample(j, 2);
    
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
%     Pdeim=deim(U{4});
%     Pdeim_ext = 4*(Pdeim-1)+(1:4); Pdeim_ext=sort(Pdeim_ext(:));
%     PTU=U_pod(Pdeim_ext,:);
%     PTUinv=U_pod(Pdeim_ext,:)\eye(size(U_pod,2));
%     meanU_v=mean(U{4},1);
%     PTtau=tau(Pdeim);
    
%     tic;
%     [time_deim,x_deim_proj]=RK3(@(t,x) PTUinv*circadian_rhs_colloc(t,x,params,PTtau,PTU,meanU_v), [0 final_time], U_pod'*x0, dt);
%     x_deim=x_deim_proj*U_pod';
%     elapsed_time_deim(j)=toc;
%     [synchrony_prom(:,j),sync_param_prom(j),avg_gene_conc_prom(:,j),spectral_ampl_fact_prom(j),estimated_period_prom(j),order_param_X_prom(:,j)] = compute_QoIs(time_prom,x_prom,num_n,n_asympt,params.omega,params.L0);
end

full_outp.time = time_full;
full_outp.x = x_full;
full_outp.synchrony = synchrony_full;
full_outp.sync_param = sync_param_full;
full_outp.avg_gene = avg_gene_conc_full;
full_outp.spectral = spectral_ampl_fact_full;
full_outp.period = estimated_period_full;
full_outp.order_param = order_param_X_full;
full_outp.elapsed_time = elapsed_time_full;

ROM_outp.time = time_rom;
ROM_outp.x = x_rom;
ROM_outp.synchrony = synchrony_rom;
ROM_outp.sync_param = sync_param_rom;
ROM_outp.avg_gene = avg_gene_conc_rom;
ROM_outp.spectral = spectral_ampl_fact_rom;
ROM_outp.period = estimated_period_rom;
ROM_outp.order_param = order_param_X_rom;
ROM_outp.elapsed_time = elapsed_time_rom;

pROM_outp.time = time_prom;
pROM_outp.x = x_prom;
pROM_outp.synchrony = synchrony_prom;
pROM_outp.sync_param = sync_param_prom;
pROM_outp.avg_gene = avg_gene_conc_prom;
pROM_outp.spectral = spectral_ampl_fact_prom;
pROM_outp.period = estimated_period_prom;
pROM_outp.order_param = order_param_X_prom;
pROM_outp.elapsed_time = elapsed_time_prom;

end

