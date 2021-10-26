clear all
clc
close all

%Constants%%%
G = 6.6742*10^(-20); %gravitational constant [km^3/(kg*s^2)]
m1 = 5.972*10^24; %mass of Earth [kg]
mu = G*m1; %gravitational parameter [km^3/s^2]
%%%%%%%%%%%%

% orbital elements at initial time
a=6897.70030; %semi-major axis [kilometers] 
enorm=0.0017946; %eccentricity. 
inclination=16; %inclination [degrees]
RAAN=162; %  right ascension of the ascending node [degrees]
arg_per=29; % argument of perigee [degrees]
true_anomaly=82; %True anomaly [degrees]

%converting orbital parameters to ECI cartsian coordinates 
r=zeros(3,1); %fill r vector with zeroes
v=zeros(3,1); %fill v vector with zeroes

p=a*(1-enorm^2); %intermediate variable
q=p/(1+enorm*cosd(true_anomaly));%intermediate variable

% Creating r vector in pqw coordinates
R_pqw(1,1) = q*cosd(true_anomaly);
R_pqw(2,1) = q*sind(true_anomaly);
R_pqw(3,1) = 0;
    
% Creating v vector in pqw coordinates    
V_pqw(1,1) =-(mu/p)^.5*sind(true_anomaly);
V_pqw(2,1) =((mu/p)^.5)*(enorm + cosd(true_anomaly));
V_pqw(3,1) =   0;

%Solving for 313 rotation matrices
R1_i = [1 0 0; 0 cosd(inclination) sind(inclination); 0 -sind(inclination) cosd(inclination)];
R3_Om = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
R3_om = [cosd(arg_per) sind(arg_per) 0; -sind(arg_per) cosd(arg_per) 0; 0 0 1];

support_var = R3_om*R1_i*R3_Om; %Intermediate variable
support_var=support_var'; %Transposed

r(:,1)=support_var*R_pqw; %Radius r [km] in ECI Cartesian
v(:,1)=support_var*V_pqw; %Velocity v [km/s] in ECI Cartesian
%%%%%FROM ORBITAL ELEMENTS TO STATE VECTOR (END)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%numerical integration of Cartesian State Vector
Y0 = [r; v];
tf = 10*3600; %s 
TSPAN = [0 tf]; %duration of integration in seconds
OPTIONS = odeset('Maxstep', 10);
%Uncomment for Problem 1
%[TOUT,YOUT] = ode45(@relacc_with_j2_HW3,TSPAN,Y0,OPTIONS); 
%Comment for problem 1
[TOUT,YOUT] = ode45(@relacc_with_j2_and_thrust_HW3,TSPAN,Y0,OPTIONS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%numerical integration of Gauss Variational equations
OrbEl0 = [a; enorm; inclination*pi/180; RAAN*pi/180; arg_per*pi/180; true_anomaly*pi/180];
[TOUTorb,YOUTorb] = ode45(@GVE_HW3,TSPAN,OrbEl0,OPTIONS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot3(YOUT(:,1),YOUT(:,2),YOUT(:,3),'r')
grid on
hold on
plot3(linspace(0,2*a),zeros(100),zeros(100),'k','linewidth',2)
plot3(zeros(100),linspace(0,2*a),zeros(100),'k','linewidth',2)
plot3(zeros(100),zeros(100),linspace(0,2*a),'k','linewidth',2)
axis equal
load('topo.mat', 'topo', 'topomap1');
[X,Y,Z] = sphere;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.CData = topo;
s=surface(6378*X,6378*Y,6378*Z,props);
xlabel('x [km]');ylabel('y [km]');zlabel('z [km]')

%orbital elements, latitude and longitude (to be corrected with Earth's
%rotation)
[size1 size2] = size(YOUT);
for i = 1:size1
    [a_out(i), enorm(i), inclination(i), RAAN(i), arg_per(i), true_anomaly(i)] = SV2OE_HW3(YOUT(i,1:3), YOUT(i,4:6));
end


figure(2)
plot(TOUT/3600, a_out, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,1),'r')
xlabel('time [hr]')
ylabel('semi-major axis a [km]')
legend("GVE's","Converted from Cartesian")

enorm = wrapToPi(enorm);
figure(3)
plot(TOUT/3600, enorm, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,2),'r')
xlabel('time [hr]')
ylabel('eccentricity e')
legend("GVE's","Converted from Cartesian")

inclination = wrapToPi(inclination);
figure(4)
plot(TOUT/3600, inclination*180/pi, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,3)*180/pi,'r')
xlabel('time [hr]')
ylabel('inclination i [deg]')
legend("GVE's","Converted from Cartesian")

RAAN = wrapToPi(RAAN);
figure(5)
plot(TOUT/3600, RAAN*180/pi, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,4)*180/pi,'r')
xlabel('time [hr]')
ylabel('RANN \Omega [deg]')
legend("GVE's","Converted from Cartesian")

arg_per = wrapToPi(arg_per);
figure(6)
plot(TOUT/3600, arg_per*180/pi, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,5)*180/pi,'r')
xlabel('time [hr]')
ylabel('argument of perigee \omega [deg]')
legend("GVE's","Converted from Cartesian")

figure(7)
plot(TOUT/3600, true_anomaly*180/pi, 'b')
grid on
hold on
plot(TOUTorb/3600, YOUTorb(:,6)*180/pi,'r')
xlabel('time [hr]')
ylabel('true anomaly \theta [deg]')
legend("GVE's","Converted from Cartesian")