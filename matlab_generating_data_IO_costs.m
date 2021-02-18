clear all

addpath '~/gams34.1_linux_x64_64_sfx'

cd '~/IO_GNEP'

load("data_for_matlab_original_costs.mat")

for i = 1:num_trials
    load(sprintf("IO_costs_%d.mat",i))
    num_arcs = double(num_arcs);
    num_nodes = double(num_nodes);
    num_players = double(num_players);
    alpha = double(alpha);    
    
    iwgdx(sprintf('IO_traffic_data_gdx_iteration_%d',i),'IO_costs_c_vec',...
        'IO_costs_c_hat_vec','node_arc_incidence','num_arcs','num_nodes',...
        'num_players','alpha')
end