
$gdxin traffic_data_gdx_iteration_%NUMBER%.gdx
$load

parameters
alpha 
$load alpha = alpha
num_arcs 
$load num_arcs = num_arcs
num_nodes
$load num_nodes = num_nodes
num_players
$load num_players = num_players
;

*https://support.gams.com/gams:a_scalar_drives_the_length_of_a_set
*https://www.gams.com/latest/docs/UG_DollarControlOptions.html#DOLLAReval
$eval NUM_ARCS num_arcs 
$eval NUM_NODES num_nodes
$eval NUM_PLAYERS num_players

set
i number of arcs /1*%NUM_ARCS%/
j number of nodes /1*%NUM_NODES%/
p number of players /1*%NUM_PLAYERS%/
;

alias (p,p1)
alias (j,od1,od2,j2,j3);
alias (i,i2);

positive variables
x(p,i)
u(i);

variables
v(p,j)
;

*https://forum.gamsworld.org/viewtopic.php?t=9414
parameter
c(i)
$load c = c_vector
c_hat(i)
$load c_hat = c_hat_vector
f(j)
D(j,i)
$load D = node_arc_incidence
;

equations
x_inequals(p,i),joint_constraints(i),equalities_x(p,j);
x_inequals(p,i) .. c(i)*(sum(p1,x(p1,i))) + c(i)*x(p,i) + c_hat(i) + sum(j,D(j,i)*v(p,j)) + u(i) =g= 0;
joint_constraints(i) .. alpha - (sum(p1,x(p1,i))) =g= 0;
equalities_x(p,j) .. sum(i,D(j,i)*x(p,i)) - f(j) =e= 0;

model traffic_problem /x_inequals.x,joint_constraints.u,equalities_x.v/;

File traffic_results /traffic_results_original_costs_%NUMBER%.csv/;
traffic_results.nd = 9;
put traffic_results;
put "data" /;
*https://www.gams.com/latest/docs/UG_FlowControl.html
*https://www.gams.com/latest/docs/UG_Put.html
loop((od1,od2) $ (od1.val<>od2.val),
loop(j2,f(j2) = 0);
f(od1) = -1;
f(od2) = 1;
*Put all the variables equal to 0
loop((p,i),x.l(p,i) = 0.5);
loop(i,u.l(i) = 0.5);
loop((p,j),v.l(p,j) = 0.5);
option mcp = PATH;
traffic_problem.optfile = 1;
solve traffic_problem using mcp;
*Try a Couple of Different Starting Points if Did not Get to Optimality
if (traffic_problem.modelStat>1,
    loop((p,i),x.l(p,i) = 0.75);
    loop(i,u.l(i) = 0.75);
    loop((p,j),v.l(p,j) = 0.75);
    traffic_problem.optfile = 1;
    solve traffic_problem using mcp;
);
if (traffic_problem.modelStat>1,
    loop((p,i),x.l(p,i) = 1);
    loop(i,u.l(i) = 1);
    loop((p,j),v.l(p,j) = 1);
    traffic_problem.optfile = 1;
    solve traffic_problem using mcp;
);
abort$(traffic_problem.modelStat>1) "DID not solve to optimality";
loop((p,i),put x.l(p,i)/);
put "data" /   
);