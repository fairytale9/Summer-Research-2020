%---Based on the model of solar irradiance I have bulit
%   Calculate best tilt, optimized thickness and maximum energy captured
%   Considering the intensity and the spectrum---%


S0 = readmatrix('AM1.5__full_spectrum.xlsx');
lambda0 = S0(:,1);
AM1_50 = S0(:,3);

I = trapz(lambda0,AM1_50); %power of AM1.5 spectrum

S = readmatrix('AM1.5_spectrum.xlsx');
lambda = S(:,1);
AM1_5 = S(:,2);
Zlambda = S(:,3);

Q = readmatrix('parameters.xlsx');
alpha = Q(:,1);
nlambda = Q(:,2);

Z = readmatrix('tilt.xlsx');
Y = readmatrix('H.xlsx');

C = 56.804;

E = (1:91)';
F = (1:91)';
G = (1:91)';
max = zeros(240,365);

for i=1:90
    latitude = (i-1)*pi/180;
    besttilt = 0;
    bestH = 0;
    default = 0;
    N = 0:364;
    declination = -asin(0.39779*cos(0.98565*(N+10)*pi/180)+1.914*sin(0.98565*(N-2)*pi/180)*pi/180);
    
    JD = 2458850+N;
    n = JD-2451545;
    g = 357.528+0.9856003*n;
    for column=1:365
        while(g(column)>360)
                g(column)=g(column)-360;
        end
    end
    distance = (1.00014-0.01671*cos(g*pi/180)-0.00014*cos(2*g*pi/180))*1.496E11;
    r = 1.496E11;
    I0 = I*(r./distance).^2;
    
    h = (-11.9:0.1:12)';
    hourangle = h*15*pi/180;
    zenithangle = acos(sin(latitude)*sin(declination)+cos(latitude)*cos(declination).*cos(hourangle));
    azimuthangle = acos((sin(declination)-cos(zenithangle)*sin(latitude))./(sin(zenithangle)*cos(latitude))); %at latitude of 90, use other forms                                     
    azimuthangle(h>0,:) = 2*pi-acos((sin(declination)-cos(zenithangle(h>0,:))*sin(latitude))./(sin(zenithangle(h>0,:))*cos(latitude)));                                      
    AM = 1./(cos(zenithangle)+0.50572*(96.07995-zenithangle*180/pi).^(-1.6364));
    for tilt=Z(i)*pi/180-pi/3600:pi/18000:Z(i)*pi/180+pi/3600
        cospsi = cos(tilt).*cos(zenithangle)-sin(tilt).*sin(zenithangle).*cos(azimuthangle);
        cospsi(zenithangle>pi/2) = 0;
        cospsi(cospsi<0) = 0;
        P0 = I0.*cospsi;
        for H=Y(i)-0.5:0.1:Y(i)+0.5            
            Ilambda = 1./(1+1./(4*nlambda.^2.*alpha*H/10000));
            J0 = 5.1684E-19*H;                      
                                                              
            for row=1:240
                for column=1:365              
                    if P0(row,column)>0
                        IX = AM1_5.*exp(Zlambda.*(AM(row,column)-1.5));
                        Slambda = IX.*lambda/1239.8;                                            
                        SI = Ilambda.*Slambda;                                 
                        JL = trapz(lambda,SI);
                        A = real(JL*P0(row,column)/1000/J0);
                        func = @(x)exp(x)*(x+1)-A;
                        x = fzero(func,[-1 41]);
                        max(row,column) = A*J0*x^2/(C*(x+1));
                    else
                        max(row,column) = 0;
                    end
                end
            end           
            totalenergy = sum(max,'all')*0.1/1000;
            
            if totalenergy>default
                besttilt = tilt;
                default = totalenergy;
                bestH = H;
            end                    
        end
    end
    E(i) = besttilt*180/pi;
    F(i) = default;
    G(i) = bestH;
end

writematrix(E,'best_tilt.xlsx')
writematrix(F,'total_energy.xlsx')
writematrix(G,'optimized_H.xlsx') 
