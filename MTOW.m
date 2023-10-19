
function [MTOW] = MTOW(S_lo,P_watts,span,chord,cl_max)
% using S_lo = 20ft, P_watts = 3000, span = 2.5ft, chord = .5ft, cl_max =
% .89
S = span*chord;
W = 5;
err = 10;
    while err > .1 
        W_guess = W;
        V_lo = 1.2*sqrt(2*W_guess/(.00238*2*S*cl_max));
        P_ft_lb = (550/.7456)*(P_watts/1000);
        T = .9*P_ft_lb/V_lo;
        [L,D] = Lift_calcs_cr001sm(V_lo*.7,chord,span,.00238,10);
        L = L*2;
        D = D*2 + .5*.00238*((V_lo*.7)^2)*.6;
        resistance = T-D;
        W = sqrt((S_lo*(32.2*.00238*S*cl_max*resistance))/1.44);
        err = abs(W_guess-W);
    end
    MTOW = W;
    %S_lo = (1.44*(W^2))/(32.2*.00238*S*cl_max*resistance);
end