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
P_drive_test = .5; % kW
P_charge_test = .5; % kW
car_capacity = 10;
d_margin = .90; % since the equality constraints never hold

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
traj.distances = {[0.1 25 10 30 60], [0.1 15 30 60]}; % cell array of the edge weights
traj.cum_distances = {[0 25 35 65], [0 15 45]};

%% Set up the MIQP problem (YALMIP)

E = sdpvar(C,N+1);
d = sdpvar(C,N+1);
y = binvar(C,N+1);
gam = binvar(C,N+1);
delta_dist = binvar(C,N+1);
% agg_time = binvar(N+1);
mode_sel = binvar(C,N+1,3);
mode_encode = binvar(C,N+1,3);
edge_selector = binvar(C,N+1,traj.max_edges);

% define the dynamical constraints
constraints = [];
e_dynamics = [];
gamma_logic = [];
delta_dynamics = [];
selector = [];
station_constraints = [];
for k = 1:N
    for c = 1:C
        % define the switched energy dynamics
        e_dynamics = [e_dynamics, ...
            implies(mode_sel(c,k,1),...
            [E(c,k+1) == E(c,k) + P_charge_test*Ts,...
            d(c,k+1) == d(c,k)]),...
            implies(mode_sel(c,k,2),...
            [E(c,k+1) == E(c,k),...
            d(c,k+1) == d(c,k)]),...
            implies(mode_sel(c,k,3),...
            [E(c,k+1) == E(c,k) - P_drive_test*Ts,...
            d(c,k+1) == d(c,k) + 5*Ts])
            ];
        %             d(c,k+1) == d(c,k) + delta_dist(c,k+1)*Ts])];
        
        selector = [selector, mode_sel(c,k,1) + mode_sel(c,k,2) + mode_sel(c,k,3) == 1,...
            implies(mode_sel(c,k,1), [0 <= gam(c,k) + y(c,k) <= 2,...
            gam(c,k) + y(c,k) == 2]),...
            implies(mode_sel(c,k,2), [0 <= gam(c,k) + y(c,k) <= 2,...
            gam(c,k) + not(y(c,k)) == 2]),...
            implies(mode_sel(c,k,3), [0 <= gam(c,k) + y(c,k) <= 2,...
            not(gam(c,k)) == 1])
            ];
        
        % define the logic for deciding to make a stop at node j
%         distance_set = false;
%         for j = 1:traj.edges{c}-1
%             distance_set = distance_set...
%                 + (2-d_margin)*traj.distances{c}(j) <= d(c,k) &&...
%                 d(c,k) <= d_margin*traj.distances{c}(j+1);
%         end
        if c == 1
            gamma_logic = [gamma_logic,...
                implies(1 - gam(c,k), ((d(c,k) >= traj.distances{c}(1)) +...
                    (d(c,k) <= traj.distances{c}(2)) + (d(c,k) >= traj.distances{c}(2)) + ...
                    (d(c,k) <= traj.distances{c}(3)) + (d(c,k) >= traj.distances{c}(3)) + ...
                    (d(c,k) <= traj.distances{c}(4)) + (d(c,k) >= traj.distances{c}(4)) + ...
                    (d(c,k) <= traj.distances{c}(5)) + (d(c,k) >= traj.distances{c}(5))) == 9),...
                implies(gam(c,k), d(c,k) == traj.distances{c})];
        elseif c == 2
            gamma_logic = [gamma_logic,...
                implies(1 - gam(c,k), ((d(c,k) >= traj.distances{c}(1)) +...
                    (d(c,k) <= traj.distances{c}(2)) + (d(c,k) >= traj.distances{c}(2)) + ...
                    (d(c,k) <= traj.distances{c}(3)) + (d(c,k) >= traj.distances{c}(3)) + ...
                    (d(c,k) <= traj.distances{c}(4)) + (d(c,k) >= traj.distances{c}(4))) == 7),...
                implies(gam(c,k), d(c,k) == traj.distances{c})];
        end
        %
        % define how delta := e_{ij}/t_{ij} switches
        % NEED TO CHECK ON THE "=" CASES SINCE THEY'RE PREFERABLY STRICT
        %         n = 1;
        %         delta_dynamics = [delta_dynamics, sum(edge_selector(c,k,:)) == 1];
        %         while n < traj.edges{c}
        %             delta_dynamics = [delta_dynamics,...
        %                 implies(edge_selector(c,k,n), ...
        %                 [traj.cum_distances{c}(n) <= d(c,k+1) <= traj.cum_distances{c}(n+1),...
        %                 delta_dist(c,k+1) == 5])];
        % %                 delta_dist(c,k+1) == traj.distances{c}(n)/travel_time(c,n)])];
        %             n = n + 1;
        %         end
        
        % station constraints: capacity
        %         station_constraints = [station_constraints,...
        %             implies(y(c,k) ==
        
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
        E(c,1) == Ecapacity, d(c,N+1) == 0,...
        sum(y(c,:))*P_charge_test + sum(~gam(c,:))*P_drive_test + E(c,N+1) - E(c,1),...
        ];
    terminal_constraints = [terminal_constraints, ...
        E_min <= E(c,N+1) <= Ecapacity,...
        0 <= d(c,N+1) <= 1.5*traj.cum_distances{c}(end),...
        mode_sel(c,N+1,1) == 1, sum(mode_sel(c,N+1,:)) == 1];
end
constraints = [constraints, e_dynamics, gamma_logic, delta_dynamics...
    initial_conditions, terminal_constraints, selector];

options = sdpsettings('verbose',0,'solver','gurobi');
p = optimize(constraints,[(100-d(:,N+1))'*(100-d(:,N+1))],options);
%% Show the results
if p.problem == 1
    error('Infeasible');
elseif p.problem ~= 0
    p
    error('The above went wrong');
else
    p
    figure(1);
    subplot(121); plot(value(E)','-x');title('Energy vs time');
    legend('Car 1','Car 2');
    subplot(122);
    plot(value(d)','-x');title('Distance vs time');
    legend('Car 1','Car 2');
end