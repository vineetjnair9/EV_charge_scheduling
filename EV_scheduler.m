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
P_drive_test = 30; % kW

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
traj.sequence = {[1 3 4 5],[2 4 5]};
traj.distances = {[25 10 30], [15 30]};

%% Set up the MIQP problem (YALMIP)
E = sdpvar(C,N+1);
d = sdpvar(C,N+1);
x = binary(C,N+1,n_nu);
y = binary(C,N+1);
% create slack variables for switching the mode
% e_md = binary(C,N+1,3);
d_md = binary(C,N+1,traj.max_edges);
gam = binary(C,N+1);
delta_dist = binary(C,N+1);
eps_counter = sdpvar(C,N+1);
current_node = sdpvar(C,N+1);

% define the dynamical constraints
e_dynamics = [];
d_dynamics = [];
eps_dynamics = [];
gamma_logic = [];
edge_logic = [];
for k = 1:N
    for c = 1:C
        % define the switched energy dynamics
        e_dynamics = [e_dynamics,...
            implies(not(gam(c,k)) + y(c,k) == 2,...
            E(c,k+1) == E(c,k) + P_charge(E(c,k),Pmax,Ec,Capacity)*Ts),...
            implies(not(gam(c,k)) + not(y(c,k)) == 2,...
            E(c,k+1) == E(c,k)),...
            implies(gam(c,k) == 1,...
            E(c,k+1) == E(c,k) - P_drive_test*Ts)];
%             implies(gamma == 1, E(c,k) - P_drive(v_ij(c,k))*Ts)
        
        % define the switched distance dynamics
        d_dynamics = [d_dynamics,...
            implies(not(gam(c,k)) + y(c,k) == 2,...
            d(c,k+1) == d(c,k)),...
            implies(not(gam(c,k)) + not(y(c,k)) == 2,...
            d(c,k+1) == d(c,k)),...
            implies(gam(c,k) == 1,...
            d(c,k+1) == d(c,k) + delta_dist(c,k+1)*Ts)];
        
        % define the specific edge counter
        eps_dynamics = [eps_dynamics,...
            implies(eps(c,k) >= current_node(c),...
            eps_counter(c,k+1) == 0),...
            implies(eps(c,k) < current_node(c,k+1),...
            eps_counter(c,k+1) == delta_dist(c)*Ts)];
        
        
    end
end