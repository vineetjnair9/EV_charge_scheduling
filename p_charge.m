function power = p_charge(E,Pmax,Ec,Ecapacity)

if E <= Ec
    power = Pmax;
else
    power = -Pmax*E /(Ecapacity-Ec) + Pmax*Ecapacity/(Ecapacity-Ec);
end

end