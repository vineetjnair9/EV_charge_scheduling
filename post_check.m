clc
ctr = 0;
n_nodes = 5;
observed = [];
for k = 1:N
    for c = 1
        total = 0;
%         disp('Current location-----');
        for i = 1:n_nodes
            next = (i-1)*n_nodes;
            xi_idx = i:n_nodes:n_nodes^2;
            
            if value(sum(xi{c}(xi_idx,k)))
                disp([num2str(i) ' | ' num2str(xi_idx) ' | ' num2str(next+1:next+n_nodes)]);
                disp('From:');
                disp(find(value(xi{c}(:,k)) == 1));
                disp('To: ');
                disp(find(value(xi{c}(:,k+1)) == 1));
                disp(['Gamma: ' num2str(value(gam(c,k))) ' | ' 'xi: ' num2str(value(sum(xi{c}(xi_idx,k))))]);
            end
            
            if value(gam(c,k)) && value(sum(xi{c}(xi_idx,k)))
                E(:,i)
                observed = [observed, value(xi{c}(next+1:next+n_nodes,k+1)'*E(:,i))];
            end
        end
            
    end
end
disp(observed');
