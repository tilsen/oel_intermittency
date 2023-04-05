function [P] = model_default_params()

%% simulation params
P(1).name = 'T';
P(end).val = 5;                  
P(end).type = 'simulation';
P(end).desc = 'maximal run time';

P(end+1).name = 'stop_on_completion';
P(end).val = 1e-3;
P(end).type = 'simulation';
P(end).desc = 'stop if the sequence is completed';

P(end+1).name = 'dt';
P(end).val = 1e-3;
P(end).type = 'simulation';
P(end).desc = 'time step';

%%
P(end+1).name = 'Nc';
P(end).val = 100;
P(end).type = 'systems';
P(end).desc = 'number of c-systems';

P(end+1).name = 'Ns';
P(end).val = 3;
P(end).type = 'systems';
P(end).desc = 'number of s-systems';

Ns = P(ismember({P.name},'Ns')).val;
P(end+1).name = 'Sx0';
P(end).val = linspace(.5,.1,Ns);
P(end).type = 'systems';
P(end).desc = 'initial syntactic system activation';

P(end+1).name = 'S_S';
P(end).val = [0 1 0; 1 0 -1; 0 -1 0];
P(end).type = 'systems';
P(end).desc = 's-to-s coupling matrix';

P(end+1).name = 'C_S_prop';
P(end).val = 1/3;
P(end).type = 'systems';
P(end).desc = 'proportion of systems coupled to verbal S-system';

P(end+1).name = 'Ne';
P(end).val = 5;
P(end).type = 'systems';
P(end).desc = 'number of initial energy levels (usually Ns + 2 extra for ground and selection level)';

P(end+1).name = 'f';
P(end).val = 8;
P(end).type = 'systems';
P(end).desc = 'oscillator intrinsic frequency (Hz)';

P(end+1).name = 'eta_Cp';
P(end).val = 0.1;
P(end).type = 'noise amplitude';
P(end).desc = 'c-system phase noise';

P(end+1).name = 'eta_Cx';
P(end).val = 0.1;
P(end).type = 'noise amplitude';
P(end).desc = 'c-system activation noise';

P(end+1).name = 'eta_Sp';
P(end).val = 0.1;
P(end).type = 'noise amplitude';
P(end).desc = 's-system phase noise';

P(end+1).name = 'eta_Sx';
P(end).val = 0.5;
P(end).type = 'noise amplitude';
P(end).desc = 's-system activation  noise';

P(end+1).name = 'eta_Sx0';
P(end).val = 0.1;
P(end).type = 'noise amplitude';
P(end).desc = 's-system activation reset noise';

P(end+1).name = 'chi_SC';
P(end).val = 0;
P(end).type = 'activation coupling strength';
P(end).desc = 's-to-c coupling';

P(end+1).name = 'chi_SS';
P(end).val = 0;
P(end).type = 'activation coupling strength';
P(end).desc = 's-to-s coupling';

P(end+1).name = 'chi_CC';
P(end).val = 0;
P(end).type = 'activation coupling strength';
P(end).desc = 'c-to-c coupling';

P(end+1).name = 'chi_CS';
P(end).val = 4;
P(end).type = 'activation coupling strength';
P(end).desc = 'c-to-s coupling';

P(end+1).name = 'psi_SC';
P(end).val = 10;
P(end).type = 'phase coupling strength';
P(end).desc = 's-to-c coupling';

P(end+1).name = 'psi_SS';
P(end).val = 10;
P(end).type = 'phase coupling strength';
P(end).desc = 's-to-s coupling';

P(end+1).name = 'psi_CC';
P(end).val = 0;
P(end).type = 'phase coupling strength';
P(end).desc = 'c-to-c coupling';

P(end+1).name = 'psi_CS';
P(end).val = 0.05;
P(end).type = 'phase coupling strength';
P(end).desc = 'c-to-s coupling';

P(end+1).name = 'gw_ub';
P(end).val = 1.1;
P(end).type = 'global excitation potential';
P(end).desc = 'upper limit';     

P(end+1).name = 'gw_lb';
P(end).val = 0;
P(end).type = 'global excitation potential';
P(end).desc = 'lower limit'; 

P(end+1).name = 'gain_gw';
P(end).val = 0.1;
P(end).type = 'global excitation potential';
P(end).desc = 'gain'; 

P(end+1).name = 'sigma_VE';
P(end).val = 0.05;
P(end).type = 's-system excitation potential';
P(end).desc = 'width of force regions'; 

P(end+1).name = 'gain_VE';
P(end).val = 1;
P(end).type = 's-system excitation potential';
P(end).desc = 'strength of excitation potential'; 

P(end+1).name = 'gain_Ai';
P(end).val = 1.2;
P(end).type = 'annealer';
P(end).desc = 'gain of annealer in initial organization phase'; 

P(end+1).name = 'gain_Ar';
P(end).val = 6;
P(end).type = 'annealer';
P(end).desc = 'gain of annealer after initial organization phase'; 

P(end+1).name = 'Me_unocc_cost';
P(end).val = 15;
P(end).type = 'excitation monitor';
P(end).desc = 'cost of unoccupied level'; 

P(end+1).name = 'Me_multiocc_cost';
P(end).val = 15;
P(end).type = 'excitation monitor';
P(end).desc = 'cost of multiply occupied level'; 

P(end+1).name = 'Me_decay';
P(end).val = 10;
P(end).type = 'excitation monitor';
P(end).desc = 'decay rate'; 

P(end+1).name = 'sigme_VMe';
P(end).val = 0.05;
P(end).type = 'excitation monitor';
P(end).desc = 'width of excitation monitor potential'; 

P(end+1).name = 'reset_Me';
P(end).val = 0.5;
P(end).type = 'excitation monitor';
P(end).desc = 'threshold for Monitor-induced reset'; 

P(end+1).name = 'reorg_Me';
P(end).val = -0.5;
P(end).type = 'excitation monitor';
P(end).desc = 'threshold for Monitor-induced regorganization'; 

P(end+1).name = 'tau_Fx_dem';
P(end).val = 0.5*ones(1,Ns);
P(end).type = 'feedback';
P(end).desc = 'feedback threshold for syntactic systems'; 

P(end+1).name = 'tau_Sx_sel';
P(end).val = 1;
P(end).type = 'feedback';
P(end).desc = 'threshold for syntactic system feedback growth'; 

P(end+1).name = 'env_num_N';
P(end).val = 2;
P(end).type = 'environment';
P(end).desc = 'number of noun concepts activated by environment'; 

P(end+1).name = 'env_num_V';
P(end).val = 1;
P(end).type = 'environment';
P(end).desc = 'number of verb concepts activated by environment'; 

P = struct2table(P);

end