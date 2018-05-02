% uses the array form of the variables for clarity
%% Define the problem parameters
clear, close all, clc
clear('yalmip');
N = 10; % prediction horizon
Ts = 1; % sample time
n_nu = 5; % number of nodes (>= to the # of stations)
capacity_nu = 1*ones(n_nu,1); % vector of vehicle capacity per station
Ec = 5; % Power threshold (switches power to decrementing)
E_max = 24; %kWh (battery capacity)
E_min = 0;
d_max = 1000;
Pmax = 45; %kW (max charging power)
P_drive_test = 2; % kW
P_charge_test = 1; % kW
D_test = 5;

%% Graph

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
m_x = [E_min; 0];
M_x = [E_max; d_max];

M_beta_var = [1 1 1]';
M_edge = max(max(E));
m_edge = min(min(E(E>0)));


%% Initialization conditions

origin.edge = m_edge*[1 1]; % Leave as is (initialization condition)
origin.node = [2 3]; % The origin node's index
dest.edge = [1 1]; % Will be the radius that we must arrive within
dest.node = [4 4]; % The index of the node to arrive at
C = length(origin.edge); % number of cars
n_nodes = 5;
x0 = [E_max; 0];

%% Set up the MIQP problem (YALMIP) using big-M relaxation

%-----Continous variables-----
% x{car}(state,k)
x = sdpvar(repmat(2,1,C),repmat(N+1,1,C));
% eps_var(c,k)
eps_var = sdpvar(C,N+1);
% current_edge(c,k)
current_edge = sdpvar(C,N+1);

%-----Binary variables-----
% y(c,k)
y = binvar(C,N+1);
% gam(c,k)
gam = binvar(C,N+1);
% xi{c}(node,N+1)
xi = binvar(repmat(n_nodes^2,1,C),repmat(N+1,1,C));
% beta_var{car}(mode(1/2/3),k)
beta_var = binvar(repmat(3,1,C),repmat(N+1,1,C));
% rho{c}(mode(1/2),k)
rho_var = binvar(repmat(2,1,C),repmat(N,1,C));
% state_preserve(c,k)
state_preserve = binvar(C,N+1);

% ----To add still----
% delta_dist(c,k)
% delta_dist = binvar(C,N);

% Slack variable for debugging
% (Get rid of at some point)
xtra = sdpvar(C,N+1);

% define the dynamical constraints
constraints = [];
for k = 1:N
    for c = 1:C
        % disjunctive implies
        constraints = [constraints,...
            implies(2 == gam(c,k+1) + y(c,k+1),...
                    beta_var{c}(1,k+1)),...
            implies(1 == gam(c,k+1) + y(c,k+1),...
                    beta_var{c}(2,k+1)),...
            implies(0 == gam(c,k+1) + y(c,k+1),...
                    beta_var{c}(3,k+1)),...
            sum(beta_var{c}(:,k+1)) == 1,...
            y(c,k+1) <= gam(c,k+1)];
        
        % mode switching
        constraints = [constraints,...
            implies(beta_var{c}(1,k+1),...
                    x{c}(:,k+1) == x{c}(:,k) + [P_charge_test; 0]*Ts),...
            implies(beta_var{c}(2,k+1),...
                    x{c}(:,k+1) == x{c}(:,k))...
            implies(beta_var{c}(3,k+1),...
                    x{c}(:,k+1) == x{c}(:,k) + [-P_drive_test; D_test]*Ts),...
                    ];
                
        % constraint on gamma        
        constraints = [constraints,...
            implies(current_edge(c,k) <= eps_var(c,k+1),...
                    gam(c,k+1)),...
            implies(0.01 <= eps_var(c,k) <= 0.99*current_edge(c,k),...
                    1-gam(c,k+1)),...
                    ];
                
        % constraints on state_preserve for memory
%         constraints = [constraints,...
%             implies(current_edge(c,k) <= eps_var(c,k),...
%                     state_preserve(c,k+1)),...
%             implies(current_edge(c,k) >= 1.01*eps_var(c,k),...
%                     1-state_preserve(c,k+1))];
                
        % epsilon value (both cases)        
        constraints = [constraints,...
            implies(gam(c,k+1),...
                    eps_var(c,k+1) == 0),...
            implies(1-gam(c,k+1),...
                    eps_var(c,k+1) == eps_var(c,k) + D_test*Ts)];
                
        % stay on current edge (gamma = 0)
        constraints = [constraints,...
            implies(1-gam(c,k+1),...
                    [current_edge(c,k+1) == current_edge(c,k),...
                    xi{c}(:,k+1) == xi{c}(:,k)]),...
                    ];  

        
        % create the switching conditions for edges
        for i = 1:n_nodes
            % Consider the cases where a switch must be made due to reaching a node
            % There are n of these cases since this can occur at each node in a
            % network.
            next = (i-1)*n_nodes;
            xi_idx = i:n_nodes:n_nodes^2;
%             disp([num2str(i) ' | ' num2str(xi_idx) ' | ' num2str(next+1:next+n_nodes)]);

            constraints = [constraints,...
                implies(gam(c,k+1) + sum(xi{c}(xi_idx,k))== 2,...
                        current_edge(c,k+1) == xi{c}(next+1:next+n_nodes,k+1)'*E(:,i)),...
                        ];
%             xi{c}(next+1:next+n_nodes,k+1)'*E(:,i)
        end
       
        
        constraints = [constraints,...
                implies(current_edge(c,k) >= 1.01*eps_var(c,k),...
                        current_edge(c,k+1) == current_edge(c,k)),...
                        ];
        
        % time-invariant box constraints on the state
        constraints = [constraints,...
            sum(xi{c}(:,k+1)) == 1
            m_x <= x{c}(:,k+1) <= M_x,...
            m_edge <= current_edge(c,k+1) <= M_edge,...
            0 <= eps_var(c,k+1) <= 2*M_edge,... % looser upper bound to allow for overstep
            ];
        
    end
end

% Add the initial conditions and terminal constraints
initial_conditions = [];
terminal_constraints = [];
for c = 1:C
    
    initial_conditions = [initial_conditions,...
        sum(xi{c}(:,1)) == 1,...
        x{c}(:,1) == x0,...
        m_x <= x{c}(:,1) <= M_x,...
        m_edge <= current_edge(c,1) <= M_edge,...
        beta_var{c}(1,1) == true,...
        ];
    
    initial_conditions = [initial_conditions,...
        current_edge(c,1) == origin.edge(c),...
        xi{c}(origin.node(c),1) == true,...
        eps_var(c,1) == origin.edge(c),...
        ];

    target = [(dest.node(c)-1)*n_nodes+1:dest.node(c)*n_nodes];
    terminal_constraints = [terminal_constraints,...
        
        sum(xi{c}(target,end)) == 1
        ];
%     terminal_constraints = [terminal_constraints,...
%         x{1}(2,:) >= 15];
end

% add all of the constraint objects
constraints = [constraints,...
                initial_conditions, terminal_constraints];

% figure(4); 
% subplot(121);
% plot(constraints,x{1}(1,:),[],[],sdpsettings('relax',1));
% subplot(122);
% plot(constraints,x{1}(2,:),[],[],sdpsettings('relax',1))

options = sdpsettings('verbose',0,'solver','gurobi');
obj = (0-x{1}(2,N+1))^2 +(0-x{2}(2,N+1))^2;
disp('starting problem');
p = optimize(constraints,obj,options);
%% Show the results
if p.problem == 1
    disp(p);
    error('Infeasible');
elseif p.problem ~= 0
    disp(p);
    error('The above went wrong');
else
    disp(p);
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
        plot(round(value(current_edge(c,:)),2),'-x'); title('Current edge length');
        legend(['Car ' num2str(c)]); hold on
        xlabel('Time'); ylabel('Length (meters)');
        subplot(224);
        edge_ = [];
        for k = 1:N+1
            edge_ = [edge_, find(value(xi{c}(:,k)) >= .95)];
        end
        plot(1:length(edge_),edge_,'x'); title('Current edge'); hold on;
        xlabel('Time'); ylabel('Current edge (vector index)');
        axis([0 N+1 0 n_nodes^2])
        
        figure(2);
        subplot(2,1,c);
        for k = 1:N
            if value(beta_var{c}(1,k)) == true
                v_k = 1;
            elseif value(beta_var{c}(2,k)) == true
                v_k = 2;
            elseif value(beta_var{c}(3,k)) == true
                v_k = 3;
            end
            
            scatter(k,v_k,'kx'); hold on
            title(['Mode: Car ' num2str(c)]); grid on
            xlabel('Time');
            ylabel('Mode');
            axis([1 N+1 1 3]);
        end
        
        legend('3: Driving','2: Waiting','1: Charging');
    end
%% See some values for checking the decision-making logic    
    disp('Decision-making logic');
    gamma__y = [value(gam(1,:))',value(y(1,:))']
    
    disp('Edge value | Epsilon counter value | Gamma |  Y')
    car = 2;
    edge__eps__gam__y = round([value(current_edge(car,:))',value(eps_var(car,:))',...
        value(gam(car,:))',value(y(car,:))'])
end