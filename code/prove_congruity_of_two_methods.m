%---Prove the congruity of irradiance calculated from spectrum and 
%   empirical equations


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

AM = 1.5;
H = 100;

%H and AM -> JL (in step10)
Ilambda = 1./(1+1./(4*nlambda.^2.*alpha*H/10000));
IX = AM1_5.*exp(Zlambda.*(AM-1.5));
Slambda = IX.*lambda/1239.8;                                            
SI = Ilambda.*Slambda;                                 
JL10 = trapz(lambda,SI);

%get P0
AM1 = AM1_50.*exp(Zlambda0.*(-0.5));
P010 = trapz(lambda0,IX0);

JL10*P010
%get JL (in step5)
LT = log10(H);
JL5 = (33.11582+8.62434*LT-2.17858*LT^2+0.21387*LT^3)*10;

%get P1
I0 = 1353;
P05 = 1.1.*I0.*0.7*1.074231;
P1 = 1.1.*I0.*0.7.^(AM.^(0.678))*1.074231;

JL5*P1
