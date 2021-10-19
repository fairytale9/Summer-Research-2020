%---Calculate optimized thickness vs. AMX spectrum---%
%---Considering spectrum

S = readmatrix('step8S.xlsx');
lambda = S(:,1);
AM1_5 = S(:,2);
Zlambda = S(:,3);

Q = readmatrix('step8.xlsx');
alpha = Q(:,1);
nlambda = Q(:,2);

S0 = readmatrix('spectrum.xlsx');
lambda0 = S0(:,1);
AM1_50 = S0(:,3);
Zlambda0 = S0(:,4);

C = 56.804;
G = (1:291)';

AM = 1:0.01:3.9;

I = trapz(lambda0,AM1_50); 

for i=1:291
    IX0 = AM1_50.*exp(Zlambda0.*(AM(i)-1.5));
    P1 = trapz(lambda0,IX0);
    Tx = P1./I;
    bestH = 0;
    default = 0;
    for H=84:0.1:110
        Ilambda = 1./(1+1./(4*nlambda.^2.*alpha*H/10000));
        J0 = 5.1684E-19*H;
        IX = AM1_5.*exp(Zlambda.*(AM(i)-1.5));
        Slambda = IX.*lambda/1239.8;                                            
        SI = Ilambda.*Slambda;                                 
        JL = trapz(lambda,SI);
        
        A = JL*I*0.6/1000./J0;                  
        func = @(x)exp(x)*(x+1)-A;
        x = fzero(func,[-1 41]);
        max = A*J0*x^2/(C*(x+1));    
        efficiency = max/(1000.*Tx);
        
        if efficiency>default
            default = efficiency;
            bestH = H;
        end
    end
    G(i) = bestH;
end
writematrix(G,'H.xlsx')
