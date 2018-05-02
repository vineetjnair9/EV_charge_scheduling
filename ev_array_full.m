% uses the array form of the variables for clarity
%% Define the problem parameters
clear, close all, clc
clear('yalmip');
N = 10; % prediction horizon
Ts = 1; % sample time
C = 3; % number of cars
n_nu = 5; % number of nodes (>= to the # of stations)
capacity_nu = 1*ones(n_nu,1); % vector of vehicle capacity per station
Ec = 5; % Power threshold (switches power to decrementing)
E_max = 15; %kWh (battery capacity)
E_min = 0;
d_max = 100;
Pmax = 45; %kW (max charging power)
P_drive_test = 3; % kW
P_charge_test = 1; % kW
D_test = 5;


%% Initialization conditions and graph
start.edge = [0 0 0];
start.node = [2 3 3];
n_nodes = 5;
x0 = [E_max; 0];

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

%% m and M to define convex hull
% If m is not defined for a pair it is taken to be zero
m_x = [E_min; -1e-10];
M_x = [E_max; d_max];
M_beta_var = [1 1 1]';
M_edge = max(max(E));
soft = 2;

%% Set up the MIQP problem (YALMIP) using big-M relaxation

% x{car}(state,k)
x = sdpvar(repmat(2,1,C),repmat(N+1,1,C));
% y(c,k)
y = binvar(C,N);
% gam(c,k)
gam = binvar(C,N);
% beta_var{car}(mode,k)
beta_var = binvar(repmat(3,1,C),repmat(N+1,1,C));
% xhi{c}{k}(node_from,node_to)
xhi = binvar(repmat(n_nodes^2,1,C),repmat(N+1,1,C));
% delta_dist(c,k)
delta_dist = binvar(C,N);
% eps_var(c,k)
eps_var = sdpvar(C,N+1);
% current_edge(c,k)
current_edge = sdpvar(C,N+1);

% slack variable
xtra = sdpvar(C,N+1);

% define the dynamical constraints
constraints = [];
states = [];
logic = [];
con = [];
for k = 1:N
    for c = 1:C
        
        % define the dynamics in big-M
        states = [states,...
            (m_x - M_x).*(1-beta_var{c}(1,k)) <= - x{c}(:,k+1) + x{c}(:,k) + [P_charge_test; 0]*Ts;
            (m_x - M_x).*(1-beta_var{c}(1,k)) <=   x{c}(:,k+1) - x{c}(:,k) - [P_charge_test; 0]*Ts;
            (m_x - M_x).*(1-beta_var{c}(2,k)) <= - x{c}(:,k+1) + x{c}(:,k);
            (m_x - M_x).*(1-beta_var{c}(2,k)) <=   x{c}(:,k+1) - x{c}(:,k);
            (m_x - M_x).*(1-beta_var{c}(3,k)) <= - x{c}(:,k+1) + x{c}(:,k) + [-P_drive_test; D_test];
            (m_x - M_x).*(1-beta_var{c}(3,k)) <=   x{c}(:,k+1) - x{c}(:,k) - [-P_drive_test; D_test],...
            beta_var{c}(1,k) + beta_var{c}(2,k) + beta_var{c}(3,k) == 1
            m_x <= x{c}(:,k+1) <= M_x,...
            ];
            
       % define the switching logic, i(t)
        logic = [logic,...
            0 <= (1-beta_var{c}(1,k)) <= M_beta_var(1)*(2 - gam(c,k) - y(c,k));
            0 <= (1-beta_var{c}(2,k)) <= M_beta_var(2)*abs(1 - gam(c,k) - y(c,k));
            0 <= (1-beta_var{c}(3,k)) <= M_beta_var(3)*(gam(c,k))];
        
        constraints = [constraints,...
            current_edge(:,k+1) == 10];
        
        % add the logic for gamma and epsilon
        logic = [logic,...
            (0-M_edge)*(1-gam(c,k)) <= -current_edge(c,k) + eps_var(c,k),...
            (0-M_edge)*(gam(c,k)) <= 0.99*current_edge(c,k) - eps_var(c,k),...
            ];
        logic = [logic,...
            (0-soft*M_edge)*(1-gam(c,k)) <= -eps_var(c,k+1),...
            (0-soft*M_edge)*(1-gam(c,k)) <=  eps_var(c,k+1),...
            (0-soft*M_edge)*(gam(c,k)) <= -eps_var(c,k+1) + eps_var(c,k) + D_test*Ts + xtra(c,k),...
            (0-soft*M_edge)*(gam(c,k)) <=  eps_var(c,k+1) - eps_var(c,k) - D_test*Ts + xtra(c,k),...
            ];
        
        % create the switching conditions for edges
        for i = 1:n_nodes
            next = (i-1)*n_nodes;
            % Consider the cases where a switch must be made due to reaching a node
            % There are n of these cases since this can occur at each node in a
            % network.
            xhi_idx = i:n_nodes:n_nodes^2;
%             disp([i xhi_idx next+1:next+n_nodes])
%             con = [con,...
%                 (0-M_edge)*(2-gam(c,k)-sum(xhi{c}(xhi_idx,k)))...
%                 - xhi{c}(next+1:next+n_nodes,k+1)*E(i,:) <= -current_edge(c,k+1),...
%                 (0-M_edge)*(2-gam(c,k)-sum(xhi{c}(xhi_idx,k)))...
%                 + xhi{c}(next+1:next+n_nodes,k+1)*E(i,:) <= current_edge(c,k+1),...
%                 ];

        end
        
        % take care of the case where it's not at a node
%         con = [con,...%         con = [con,...
%             (0-M_edge)*gam(c,k) <= -current_edge(c,k+1) + current_edge(c,k),...
%             (0-M_edge)*gam(c,k) <=  current_edge(c,k+1) - current_edge(c,k)];

%             (0-M_edge)*gam(c,k) <= -current_edge(c,k+1) + current_edge(c,k),...
%             (0-M_edge)*gam(c,k) <=  current_edge(c,k+1) - current_edge(c,k)];
        
        % time invariant box constraints on the state
        constraints = [constraints,...
            current_edge(c,k+1) >= 1,...
            0 <= eps_var <= soft*M_edge,...
            sum(xhi{c}(:,k+1)) == 1];
        
    end
end

% AdZ the initial conditions and terminal constraints
initial_conditions = [];
terminal_constraints = [];
for c = 1:C
    states = [states,...
        x{c}(:,1) == x0,...
        m_x <= x{c}(:,1) <= M_x,...
        ];
    
    initial_conditions = [initial_conditions,...
        sum(xhi{c}(:,1)) == 1,...
        current_edge(c,1) == start.edge(c),...
        xhi{c}(start.node(c),1) == true,...
        eps_var(c,1) == start.edge(c)
        ];
    % 
    
%     k = N+1;
%     terminal_constraints = [terminal_constraints, ...
%         beta_var{c}(1,k) + beta_var{c}(2,k) + beta_var{c}(3,k) == 2,...
%         ];
    
end

% add all of the constraint objects
constraints = [constraints, con,...
    states, logic, initial_conditions, terminal_constraints];

% figure(4); 
% subplot(121);
% plot(constraints,x{1}(1,:),[],[],sdpsettings('relax',1));
% subplot(122);
% plot(constraints,x{1}(2,:),[],[],sdpsettings('relax',1))

options = sdpsettings('verbose',0,'solver','gurobi');
obj = (100-x{1}(2,N+1))^2 +(100-x{2}(2,N+1))^2;
disp('Starting problem');
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
    
      
    for c = 1:2
        figure(1);
        subplot(221);
        plot(states{c}(1,:),'-x'); title('Energy');
        legend(['Car ' num2str(c)]); hold on
        xlabel('Time'); ylabel('Energy in battery (kWh)');
        subplot(222); 
        plot(states{c}(2,:),'-x'); title('Cumulative distance'); 
        legend(['Car ' num2str(c)]); hold on
        xlabel('Time'); ylabel('Distance (meters)');
        subplot(223);
        plot(value(current_edge(c,:)),'-x'); title('Current edge length');
        legend(['Car ' num2str(c)]); hold on
        xlabel('Time'); ylabel('Length (meters)');
        subplot(224);
        edge_ = [];
        for k = 1:N+1
            edge_ = [edge_, find(value(xhi{c}(:,k)) == 1)];
        end
        plot(1:N+1,edge_,'x'); title('Current edge'); hold on;
        xlabel('Time'); ylabel('Current edge (vector index)');
        
        figure(2);
        subplot(2,1,c);
        for k = 1:N+1
            if value(beta_var{c}(1,k)) == true
                v_k = 3;
            elseif value(beta_var{c}(2,k)) == true
                v_k = 2;
            else
                v_k = 1;
            end
            
            scatter(k,v_k,'kx'); hold on
            title(['State: Car ' num2str(c)]); grid on
            axis([1 N+1 1 3]);
            
            
        end
        legend('3: Charging','2: Waiting','1: Driving');
    end
%% See some values for checking the decision-making logic    
    disp('Decision-making logic');
    gamma__y = [value(gam(1,:))',value(y(1,:))']
    
    disp('Edge value | Epsilon counter value | Gamma')
    edge__eps__gamma = [value(current_edge(1,1:end-1))',value(eps_var(1,1:end-1))', value(gam(1,:))']
end