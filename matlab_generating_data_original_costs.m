clear all

addpath '~/gams34.1_linux_x64_64_sfx'

cd '~/IO_GNEP'

rng(5)

load("data_for_matlab_original_costs.mat")

for i = 1:num_trials
    c_vector = unifrnd(double(lowerbound_c),double(upperbound_c),1,double(num_arcs));
    c_hat_vector = unifrnd(double(lowerbound_chat),double(upperbound_chat),1,double(num_arcs));
    num_arcs = double(num_arcs);
    num_nodes = double(num_nodes);
    num_players = double(num_players);
    alpha = double(alpha);
    
    iwgdx(sprintf('traffic_data_gdx_iteration_%d',i),'c_vector',...
        'c_hat_vector','node_arc_incidence','num_arcs','num_nodes',...
        'num_players','alpha')
    save(sprintf('costs_iteration_%d',i),'c_vector','c_hat_vector')

end