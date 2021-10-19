%---Based on the model of solar irradiance I have bulit
%   Calculate best tilt, optimized thickness and maximum energy captured
%   Only considering the intensity not the spectrum---%

%---load data calculated without considering the property of solar panels
rough_tilt = readmatrix('tilt.xlsx');
rough_thickness = readmatrix('H.xlsx');

P0 = 1000;
C = 56.804;

D = 0:90; %latitude
E = (1:91)'; %store the best tilt at each latitude
F = (1:91)'; %store the maximum energy captured at each latitude
G = (1:91)'; %store the optimized thickness at each latitude
max = zeros(240,365);

for i = 1:90
    latitude = D(i)*pi/180;
    besttilt = 0;
    bestH = 0;
    default = 0;
    for tilt = rough_tilt(i)*pi/180-pi/9000:pi/18000:rough_tilt(i)*pi/180+pi/9000
        for H = rough_thickness(i)-0.2:0.1:rough_thickness(i)+0.2
            LT = log10(H);
            JL = (33.11582+8.62434*LT-2.17858*LT^2+0.21387*LT^3)*10;
            J0 = 5.1684E-19*H;
            
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
            %---adding the influence of the changing of distance between sun and earth---%
            distance = (1.00014-0.01671*cos(g*pi/180)-0.00014*cos(2*g*pi/180))*1.496E11;
            r = 1.496E11;
            I0 = 1.353*(r./distance).^2;
            
            %---matrix calculation to accelerate the program---%
            h = (-11.9:0.1:12)';
            hourangle = h*15*pi/180;
            zenithangle = acos(sin(latitude)*sin(declination)+cos(latitude)*cos(declination).*cos(hourangle));
            azimuthangle = acos((sin(declination)-cos(zenithangle)*sin(latitude))./(sin(zenithangle)*cos(latitude)));                                                 
            azimuthangle(h>0,:) = 2*pi-acos((sin(declination)-cos(zenithangle(h>0,:))*sin(latitude))./(sin(zenithangle(h>0,:))*cos(latitude)));                                      
            AM = 1./(cos(zenithangle)+0.50572*(96.07995-zenithangle*180/pi).^(-1.6364));            
            I = 1.1.*I0.*0.7.^(AM.^(0.678));
            I(zenithangle>pi/2) = 0;
            cospsi = cos(rough_tilt).*cos(zenithangle)-sin(rough_tilt).*sin(zenithangle).*cos(azimuthangle);
            cospsi(cospsi<0) = 0;                                                      
            result = cospsi.*I;
            P1 = result*1000;
            
            A = JL.*P1/P0/J0;            
            for row=1:240
                for column=1:365
                    if P1(row,column)>0
                        AA = real(A(row,column));               
                        func = @(x)exp(x)*(x+1)-AA;
                        x = fzero(func,[-1 40]);
                        max(row,column) = AA*J0*x^2/(C*(x+1));
                    else
                        max(row,column) = 0;
                    end
                end
            end                     
            totalenergy = sum(max,'all')*0.1/1000;   
            
            if  totalenergy>default
                besttilt = rough_tilt;
                default = totalenergy;
                bestH = H;
            end
        end
    end
    E(i) = besttilt*180/pi;
    F(i) = default;
    G(i) = bestH;
end

writematrix('best_tilt.xlsx',E)
writematrix('total_energy.xlsx',F)
writematrix('optimized_H.xlsx',G)
