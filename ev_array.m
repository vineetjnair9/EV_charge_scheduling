% uses the array form of the variables for clarity
%% Define the problem parameters
clear, close all, clc
clear('yalmip');
N = 15; % prediction horizon
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


m_x = [E_min; -1e-10];
M_x = [E_max; d_max];
M_theta_var = [1 1 1]';
m_leg = -15; % min and max distances for edges
M_leg = 30;
x0 = [E_max; 0];
%% Set up the MIQP problem (YALMIP)

% x{car}(state,k)
x = sdpvar(repmat(2,1,C),repmat(N+1,1,C));
y = binvar(C,N);
gam = binvar(C,N);

% theta_var{car}(mode,k)
theta_var = binvar(repmat(3,1,C),repmat(N+1,1,C));
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
            (m_x - M_x).*theta_var{c}(1,k) <= - x{c}(:,k+1) + x{c}(:,k) + [P_charge_test; 0]*Ts;
            (m_x - M_x).*theta_var{c}(1,k) <=   x{c}(:,k+1) - x{c}(:,k) - [P_charge_test; 0]*Ts;
            (m_x - M_x).*theta_var{c}(2,k) <= - x{c}(:,k+1) + x{c}(:,k);
            (m_x - M_x).*theta_var{c}(2,k) <=   x{c}(:,k+1) - x{c}(:,k);
            (m_x - M_x).*theta_var{c}(3,k) <= - x{c}(:,k+1) + x{c}(:,k) + [-P_drive_test; 5];
            (m_x - M_x).*theta_var{c}(3,k) <=   x{c}(:,k+1) - x{c}(:,k) - [-P_drive_test; 5],...
            theta_var{c}(1,k) + theta_var{c}(2,k) + theta_var{c}(3,k) == 2
            ];
        
            
       % define the switching logic, i(t)
        constraints = [constraints,...
            0 <= theta_var{c}(1,k) <= M_theta_var(1)*abs(2 - gam(c,k) - y(c,k));
            0 <= theta_var{c}(2,k) <= M_theta_var(2)*abs(1 - gam(c,k) - y(c,k));
            0 <= theta_var{c}(3,k) <= M_theta_var(3)*abs(0 - gam(c,k) - y(c,k))];
        
%         for j = 1:traj.distances{c}
%             constraints = [constraints,...
%                 0 <= gam(c,k) <= M_gam*(x(c,k,2) - traj.distances{c}(j))
    
        constraints = [constraints,...
            m_x <= x{c}(:,k+1) <= M_x];
        
    end
end

initial_conditions = [m_x <= x{c}(:,1) <= M_x];
terminal_constraints = [];
for c = 1:C
    initial_conditions = [initial_conditions,...
        x{c}(:,1) == x0];
    k = N+1;
    terminal_constraints = [terminal_constraints, ...
        theta_var{c}(1,k) + theta_var{c}(2,k) + theta_var{c}(3,k) == 2
        ];
end
constraints = [constraints, ...
    initial_conditions, terminal_constraints];

figure(3); subplot(121);
plot(constraints,x{1}(1,:),[],[],sdpsettings('relax',1));
xlabel('Car #');
ylabel('Time step');
zlabel('Energy');

subplot(122);
plot(constraints,x{1}(2,:),[],[],sdpsettings('relax',1))
xlabel('Car #');
ylabel('Time step');
zlabel('Distance');
X = intersect(x{1}(1,:),x{1}(2,:))
options = sdpsettings('verbose',0,'solver','gurobi');
obj = (100-x{1}(2,N+1))^2 +(100-x{2}(2,N+1))^2;
p = optimize(constraints,obj,options);
%% Show the results
if p.problem == 1
    p
    error('Infeasible');
elseif p.problem ~= 0
    p
    error('The above went wrong');
else
    p
    states = cell(C,1);
    for k = 1:N+1
        for c = 1:C
            states{c} = [states{c},value(x{c}(:,k))];
        end         
    end
    
    figure(1);  
    for c = 1
        subplot(121);
        plot(states{c}(1,:),'-x'); title('Energy vs time');
        legend(['Car ' num2str(c)]); hold on
        subplot(122); 
        plot(states{c}(2,:),'-x'); title('Distance vs time'); 
        legend(['Car ' num2str(c)]); hold on
    end
    
    
    figure(2); 
    for k = 1:N+1
         if value(theta_var{1}(1,k)) == false
             v_k = 3;
         elseif value(theta_var{1}(2,k)) == false
             v_k = 2;
         else
             v_k = 1;
         end
             
        scatter(k,v_k,'kx'); hold on
        title('State: Car 1'); grid on
        axis([1 N+1 0 3]);
        
        
    end
    legend('3: Charging','2: Waiting','1: Driving');

    disp('Decision-making logic');
    gamma__y = [value(gam(1,:))',value(y(1,:))']
end