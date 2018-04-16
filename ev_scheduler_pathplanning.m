% function EV_scheduler(N,Ts,
%% Define the problem parameters
clear, close all, clc
clear('yalmip');
N = 30; % prediction horizon
Ts = 1; % sample time
C = 2; % number of cars
n_nu = 5; % number of nodes (>= to the # of stations)
capacity_nu = 1*ones(n_nu,1); % vector of vehicle capacity per station
Ec = 5; % Power threshold (switches power to decrementing)
Ecapacity = 15; %kWh (battery capacity)
E_min = 0;
Pmax = 45; %kW (max charging power)
P_drive_test = .1; % kW
car_capacity = 10;
margin = .30; % since the equality constraints never hold

%% Run the feas_traj function to get the sequences (ordered-sets) of edges
% % to traverse
% % Initial nodes (nu)
% 
% % Define connected graph G = (V,E,A) 
% 
% traj = feas_traj(C
% 
% Before developing the above, use a basic problem
traj.num_cars = 2; % scalar indicating the number of cars on the network
traj.edges = {4,3}; % cell array indicating the numbers of edges traversed
traj.max_edges = 4; 
traj.sequence = {[1 3 4 5],[2 4 5]}; % cell array of the node sequence
traj.distances = {[25 10 30], [15 30]}; % cell array of the edge weights
traj.cum_distances = {[0 25 35 65], [0 15 45]};

%% Set up the MIQP problem (YALMIP)
E = sdpvar(C,N+1);
d = sdpvar(C,N+1);
x = binvar(C,N+1,n_nu);
y = binvar(C,N+1);
gam = binvar(C,N+1);
delta_dist = binvar(C,N+1);
% agg_time = binvar(N+1);
e_md = binvar(C,N+1,3);
d_md = binvar(C,N+1,3);
edge_selector = binvar(C,N+1,traj.max_edges);
g = binvar(C,N+1);

% define the dynamical constraints
constraints = [];
e_dynamics = [];
d_dynamics = [];
gamma_logic = [];
delta_dynamics = [];
for k = 1:N
    for c = 1:C
        % define the switched energy dynamics
        e_dynamics = [e_dynamics, sum(e_md(c,k,:)) == 1,...
            implies(e_md(c,k,1),[gam(c,k) + y(c,k) == 2,...
            E(c,k+1) == E(c,k) + p_charge(E(c,k),Pmax,Ec,Ecapacity)*Ts]),...
            implies(e_md(c,k,2), [gam(c,k) + not(y(c,k)) == 2,...
            E(c,k+1) == E(c,k)]),...
            implies(e_md(c,k,3), [not(gam(c,k)) == 1,...
            E(c,k+1) == E(c,k) - P_drive_test*Ts])];
% % %             implies(gam(c,k) == 1, E(c,k) - P_drive(v_ij(c,k))*Ts)
%             e_dynamics = [e_dynamics, sum(e_md(c,k,:)) == 1,...
%                 implies(e_md(c,k,1), [E_min <= E(c,k) <= Ecapacity, ...
%                 E(c,k+1) == E(c,k) + p_charge(E(c,k),Pmax,Ec,Ecapacity)*Ts]),...
%                 implies(e_md(c,k,2), [E_min <= E(c,k) <= Ecapacity,...
%                 E(c,k+1) == E(c,k)]),...
%                 implies(e_md(c,k,3), [E_min <= E(c,k) <= Ecapacity,...
%                 E(c,k+1) == E(c,k) - P_drive_test*Ts])];
        
        % define the switched distance dynamics
        d_dynamics = [d_dynamics, sum(d_md(c,k,:)) == 1,...
            implies(d_md(c,k,1), [gam(c,k) + y(c,k) == 2,...
            d(c,k+1) == d(c,k)]),...
            implies(d_md(c,k,2), [gam(c,k) + not(y(c,k)) == 2,...
            d(c,k+1) == d(c,k)]),...
            implies(d_md(c,k,3), [not(gam(c,k)) == 1,...
            d(c,k+1) == d(c,k) + delta_dist(c,k+1)*Ts])];
%         d_dynamics = [d_dynamics, sum(d_md(c,k,:)) == 1,...
%             implies(d_md(c,k,1), [0 <= d(c,k) <= traj.cum_distances{c}(end),...
%             d(c,k+1) == d(c,k)]),...
%             implies(d_md(c,k,2), [0 <= d(c,k) <= traj.cum_distances{c}(end),...
%             d(c,k+1) == d(c,k)]),...
%             implies(d_md(c,k,3), [0 <= d(c,k) <= traj.cum_distances{c}(end),...
%             d(c,k+1) == d(c,k) + 5*Ts])];
        
        % define the logic for deciding to stop
        gamma_logic = [gamma_logic, sum(g(c,k)),...
            implies(g(c,k),[margin*traj.distances{c} <= d(c,k+1) <=(2-margin)*traj.distances{c},...
            gam(c,k)]),...
            implies(g(c,k), [d(c,k+1) ~= traj.distances{c},...
            ~gam(c,k)])];
            
        % define how delta := e_{ij}/t_{ij} switches 
        % NEED TO CHECK ON THE "=" CASES SINCE THEY'RE PREFERABLY STRICT
        n = 1;
        delta_dynamics = [delta_dynamics, sum(edge_selector(c,k,:)) == 1];
        while n < traj.edges{c}
            delta_dynamics = [delta_dynamics,...
                implies(edge_selector(c,k,n), ...
                [traj.cum_distances{c}(n) <= d(c,k+1) <= traj.cum_distances{c}(n+1),...
                delta_dist(c,k+1) == 5])];
%                 delta_dist(c,k+1) == traj.distances{c}(n)/travel_time(c,n)])];
            n = n + 1;
        end
        
        % general constraints
        constraints = [constraints,...
            E_min <= E(c,k) <= Ecapacity,...
            0 <= d(c,k) <= traj.cum_distances{c}(end)];
    end
%         constraints = [constraints,...
%             sum(y(c,k+1)) <= car_capacity];
end

initial_conditions = [];
terminal_constraints = [];
for c = 1:C
    initial_conditions = [initial_conditions,...
        E(c,1) == Ecapacity, d(c,N+1) == 0];
    terminal_constraints = [terminal_constraints, ...
        d(c,N+1) >= .5*traj.cum_distances{c}(end), E_min <= E(c,N+1) <= Ecapacity,...
            0 <= d(c,N+1) <= traj.cum_distances{c}(end)];
end
constraints = [constraints, e_dynamics, d_dynamics, gamma_logic, delta_dynamics...
    initial_conditions, terminal_constraints];

options = sdpsettings('verbose',0,'solver','gurobi');
p = optimize(constraints,[],options);
if p.problem == 1
    error('Infeasible');
else
    p
    figure(1); 
    subplot(121); plot(value(E)','rx');title('Energy vs time');
    subplot(122); 
    plot(value(d)','bo');title('Distance vs time');
end