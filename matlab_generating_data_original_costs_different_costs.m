clear all

addpath '~/gams34.1_linux_x64_64_sfx'

cd '~/IO_GNEP'

rng(5)

load("data_for_matlab_original_costs.mat")

for i = 1:num_trials
    num_arcs = double(num_arcs);
    num_nodes = double(num_nodes);
    num_players = double(num_players);
    alpha = double(alpha);
    node_arc_incidence = node_arc_incidence;
    
    c_matrix = unifrnd(double(lowerbound_c),double(upperbound_c),...
        double(num_players),double(num_arcs));
    c_hat_matrix = unifrnd(double(lowerbound_chat),double(upperbound_chat),...
        double(num_players),double(num_arcs));

%     c_matrix = ones(num_players,num_arcs)*10;
% 
%     for j = 1:num_players
%         random_c_numbers = randi([1,num_arcs],1,10);
%         c_matrix(j,random_c_numbers) = 2;
%     end
    
    iwgdx(sprintf('traffic_data_gdx_iteration_%d',i),'c_matrix',...
        'c_hat_matrix','node_arc_incidence','num_arcs','num_nodes',...
        'num_players','alpha')
    save(sprintf('costs_iteration_%d',i),'c_matrix','c_hat_matrix')

end