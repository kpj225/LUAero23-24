function [L,D] = Lift_calcs_cr001sm(V,Chord,Span,density,alpha)
%find lift and drag of cr001sm at various conditions
%   Detailed explanation goes here
    AR = Span/Chord; %find aspect ratio
    area = Span*Chord; 
    a0 = 7;
    alpha_0 = -2;
    cd0 = .03;
    e = .7;
    ar = (a0)/(1+(a0/(pi*e*AR)));
    ad = ar*(2*pi/360);
    Cl = ad*(alpha-alpha_0);
    q = .5*density*(V^2);
    L=Cl*q*area;
    cd = cd0 + (Cl^2)/(pi*e*AR);
    D= cd*q*area;
end