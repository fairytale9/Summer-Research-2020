%---Calculate solar cell efficiency vs. AMX spectrum

max = 1:39;
C = 56.804;

%---Use the known AM0 and AM1.5 spectrum to calculate AMX spectrum
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

H = 100;
AM = 1:39;

%---All the calculations are executed by matrix operations
AM1 = AM1_50.*exp(Zlambda0.*(-0.5));
I = trapz(lambda0,AM1); 
IX0 = AM1_50.*exp(Zlambda0.*(AM-1.5));
P1 = trapz(lambda0,IX0)
Tx = P1./I;
Ilambda = 1./(1+1./(4*nlambda.^2.*alpha*H/10000));
J0 = 5.1684E-19*H;
IX = AM1_5.*exp(Zlambda.*(AM-1.5));
Slambda = IX.*lambda/1239.8;                                            
SI = Ilambda.*Slambda;                                 
JL = trapz(lambda,SI);

A = JL./J0;
for j=1:39              
    func = @(x)exp(x)*(x+1)-A(j);
    x = fzero(func,[-1 41]);
    max(j) = A(j)*J0*x^2/(C*(x+1));
end
efficiency = max./(1000.*Tx)    
plot(AM,efficiency)
    