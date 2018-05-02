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
E_max = 15; %kWh (battery capacity)
E_min = 0;
d_max = 100;
Pmax = 45; %kW (max charging power)
P_drive_test = .5; % kW
P_charge_test = 5; % kW
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
n_nodes = 5;
A = [0 1 1 0 0;
     1 0 0 1 0;
     1 0 0 1 1;
     0 1 1 0 1;
     0 0 1 1 0];
E = [0  30 25 0  0;
     0  0  0  15 0;
     0  0  0  10 45;
     0  0  0  0  30;
     0  0  0  0  0];
E = E + E';


m_x = [E_min; 0];
M_x = [E_max; d_max];
M_theta_var = [1 1 1]';
m_leg = -15; % min and max distances for edges
M_leg = 30;
x0 = [E_max; 0];
%% Set up the MIQP problem (YALMIP)

x = sdpvar(C,N+1,2);
y = binvar(C,N);
gam = binvar(C,N);
theta_var = binvar(C,N,3,2);
delta_dist = binvar(C,N);

% define the dynamical constraints
constraints = [];
d_dynamics = [];
gamma_logic = [];
selector = [];
station_constraints = [];
for k = 1:N
    for c = 1:C
        % define the dynamics in big-M
        constraints = [constraints,...
            (m_x - M_x).*theta_var(c,k,1,1:2) <= - x(c,k+1,:) + x(c,k,:) + [P_charge_test; 0]*Ts;
            (m_x - M_x).*theta_var(c,k,1,1:2) <= x(c,k+1,:)   - x(c,k,:) - [P_charge_test; 0]*Ts;
            (m_x - M_x).*theta_var(c,k,2,1:2) <= - x(c,k+1,:) + x(c,k,:);
            (m_x - M_x).*theta_var(c,k,2,1:2) <= x(c,k+1,:)   - x(c,k,:);
            (m_x - M_x).*theta_var(c,k,3,1:2) <= - x(c,k+1,:) + x(c,k,:) + [-P_drive_test; delta_dist(c,k)];
            (m_x - M_x).*theta_var(c,k,3,1:2) <= x(c,k+1,:)   - x(c,k,:) - [-P_drive_test; delta_dist(c,k)],...
            theta_var(c,k,1,1:2) + theta_var(c,k,2,1:2) + theta_var(c,k,3,1:2) == [2;2],...
            theta_var(c,k,1,1) == theta_var(c,k,1,2),...
            theta_var(c,k,2,1) == theta_var(c,k,2,2),...
            theta_var(c,k,3,1) == theta_var(c,k,3,2),...
            ];
       % define the switching logic, i(t)
        constraints = [constraints,...
            0 <= theta_var(c,k,1,1:2) <= M_theta_var(1)*abs(2 - gam(c,k) - y(c,k));
            0 <= theta_var(c,k,2,1:2) <= M_theta_var(2)*abs(1 - gam(c,k) - y(c,k));
            0 <= theta_var(c,k,3,1:2) <= M_theta_var(3)*abs(0 - gam(c,k) - y(c,k))];
        
%         for j = 1:traj.distances{c}
%             constraints = [constraints,...
%                 0 <= gam(c,k) <= M_gam*(x(c,k,2) - traj.distances{c}(j))
    
        constraints = [constraints,...
            m_x <= x(c,k+1,:) <= M_x];
        
    end
end

initial_conditions = [];
terminal_constraints = [];
for c = 1:C
    initial_conditions = [initial_conditions,...
        ];
    terminal_constraints = [terminal_constraints, ...
        ];
end
constraints = [constraints, ...
    initial_conditions, terminal_constraints, selector];

figure(1); subplot(121);
plot(constraints,x(:,:,1),[],[],sdpsettings('relax',1));
xlabel('Car #');
ylabel('Time step');
zlabel('Energy');

subplot(122);
plot(constraints,x(:,:,2),[],[],sdpsettings('relax',1))
xlabel('Car #');
ylabel('Time step');
zlabel('Distance');
X = intersect(x(:,:,1),x(:,:,2))
% options = sdpsettings('verbose',0,'solver','gurobi');
% p = optimize(constraints,(100-x(:,N+1,2))'*(100-x(:,N+1,2)),options);
%% Show the results
if p.problem == 1
    error('Infeasible');
elseif p.problem ~= 0
    p
    error('The above went wrong');
else
    p
    figure(1); 
    subplot(121); plot(value(x(:,:,1))','-x');title('Energy vs time');
    legend('Car 1','Car 2');
    subplot(122); 
    plot(value(x(:,:,2))','-x');title('Distance vs time'); 
    legend('Car 1','Car 2');
    
    figure(2); 
    subplot(221); plot(reshape(value(theta_var(1,:,:,1)),[],3),'-x');
    title('State: Car 1, State: Energy');
    legend('Charging','Waiting','Driving');
    subplot(222); plot(reshape(value(theta_var(1,:,:,2)),[],3),'-x');
    title('State: Car 1, State: Distance');
    legend('Charging','Waiting','Driving');
    subplot(223); plot(reshape(value(theta_var(2,:,:,1)),[],3),'-x');
    title('State: Car 2, State: Energy');
    legend('Charging','Waiting','Driving');
    subplot(224); plot(reshape(value(theta_var(2,:,:,2)),[],3),'-x');
    title('State: Car 2, State: Distance');
    legend('Charging','Waiting','Driving');
end