function costs(evarray);

% Cost functions

% Assume all cars driving at a constant velocity throughout
% Program focuses only on routing them and scheduling charges
v = 45; % mi/h

%% Cost to customers (time)
% Cost = C_driving + C_waiting + C_charging;

% feas_traj only computes at the beginning, 
% don't need to redo at each timestep
Cost = 0;

% 1st find feasible trajectories for graph G(V,E,A)
traj = feas_traj(C); 

% Finding number of cars along each edge for graph 
% Still need a 3D matrix or cell array (better)
edgeCount = zeros(V, V, N/Ts); % prediction horizon N and sampling time Ts

for k = 1:N/Ts
    for i = 1:size(traj.sequence)(1,2) 
        for j = 1:size(traj.sequence{i}) - 1
            edgeCount(k,

% Driving time for all cars for whole journey from start to dest
for i = 1:size(traj.distances)(1,2) 
    for j = 1:size(traj.distances{i})
        C_driving = (traj.distances{i,j}/v)*(1 + 0.15*(
        
%% Cost to stations

end