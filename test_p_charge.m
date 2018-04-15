% Checks to see that p_charge is working as anticipated
close all
E = 0; Ec = 5; Ecap = 15; Pmax = 45;
p_v = [];
E_v = [];
while abs(Ecap-E) > .01
E_v = [E_v, E];
p = p_charge(E,Pmax,Ec,Ecap);
p_v = [p_v, p];
E = E + p*.01;
pause(.1)
end
figure(1);
plot(E_v,p_v)