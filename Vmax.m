function [outputArg1,outputArg2] = Vmax(Chord,span,W,Cd0,P)
% using chord = .5ft, span = 2.5ft, W = aircraft weight from MTOW, Cd_0 = .03, 
%P = 3000 
%P = WH*(max burst rate(C)) = 60WH * 50C
    rho = .00238;
    S = Chord*span;
    AR = span/Chord;
    e = .7;
    x0 = [150];
    %P = (550/.745)*P
    fun = @(x)root1d(x,rho,Cd0,W,e,AR,S,P);
    fsolve(fun,x0)

end
function F = root1d(Vmax,rho,Cd0,W,e,AR,S,P);
    cl = W/(2*S*(.5*.00238*Vmax^2));
    alpha = (cl/.0746) + (-2);
    [L,D] = Lift_calcs_cr001sm(Vmax,.5,2.5,rho,alpha);
    D = 2*D + .5*.00238*(Vmax^2)*.6;
    T = ((550/.7456)*(P/1000))/Vmax;
    F(1) = D-T;
end