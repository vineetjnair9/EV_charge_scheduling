% function EV_scheduler(N,Ts,
%% Define the problem parameters
N = 5; % prediction horizon
Ts = 1; % sample time
C = 4; % number of cars
n_nu = 5; % number of nodes (>= to the # of stations)
capacity_nu = 1*ones(n_nu,1); % vector of vehicle capacity per station
Ec = 5; % Power threshold (switches power to decrementing)
Ecapacity = 15; %kWh (battery capacity)
Pmax = 45; %kW (max charging power)

%% Run the feas_traj function to get the sequences (ordered-sets) of edges
% % to traverse
% % Initial nodes (nu)
% 
% % Define connected graph G = (V,E,A) 
% 
% traj = feas_traj(
% 
% Before developing the above, use a basic problem
traj.max_edges = 3;

%% Set up the MIQP problem (YALMIP)
E = sdpvar(C,N+1);
d = sdpvar(C,N+1);
x = binary(C,N+1,n_nu);
y = binary(C,N+1,n_nu);
% create slack variables for switching the mode
e_md = binary(C,N+1,3);
d_md = binary(C,N+1,traj.max_edges);

% define the dynamical constraints
dynamics = [];
for k = 1:N
    for c = 1:C
        for e_piece = 1:3
            dynamics = [dynamics, sum(e_md(c,k,:)),implies(e_md(c,k,i),...
            E(c,k+1) == E(c,k) + P_charge(E(c,k),Pmax,Ec,Capacity)*Ts,...
