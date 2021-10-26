%Orbital Paramaters with J2
clear all;
clc; close all;

%% Constants of the Problem
%Start
G = 6.6742*10^(-20); %[km^3/(kg*s^2)] gravitational constant
m1 = 5.9726*10^24; %[kg] mass of earth
mu = G*m1; %[km^3/s^2]
%Convert date to a large number representing days
JulianDate = juliandate([2017, 03, 03, 12, 00,00]);%days Julian Date (input is GMT time) - equivalent to 2017/062/12:00:00.000
%Find GST WRT JulianDate
gst = gstime(JulianDate); %[rad] Greenwich location on julian date (longitude in ECI)
omega_E = 72.9217*10^(-6); %[rad/s] Earth's Angular Velocity
%End
%% ISS Orbital Parameters
a = 6778928.17*10^(-3); %[km]
enorm = 0.0008446;
inclination = 51.46843*pi/180; %[rad]
RAAN = 203.98784*pi/180; %[rad]
arg_per = 141.71208*pi/180; %[rad]
true_anomaly = 149.88842*pi/180; %[rad]

hnorm = sqrt(a*mu*(1-enorm^2)); %[] norm of angular momentum
T = (2*pi/sqrt(mu))*a^(3/2); %[s] Period

%% Convert Orbital Parameters to State Vector
R1_i = [1 0 0; 0 cos(inclination) sin(inclination); 0 -sin(inclination) cos(inclination)];
R3_Om = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
R3_om = [cos(arg_per) sin(arg_per) 0; -sin(arg_per) cos(arg_per) 0; 0 0 1];

%Align Axis with Euler Angle Sequence
support_var = R3_om*R1_i*R3_Om; %Total rotation Matrix
x = support_var(1,1:3);
y = support_var(2,1:3);
z = support_var(3,1:3);
%taking first, second, and third rows of R to tell you x y and z
%x is e direction, z is h directon
e_orb = enorm*x;
h_orb = hnorm*z;
%Now can put everything together

rnorm = ((hnorm^2)/mu)*(1/(1+enorm*cos(true_anomaly)));
r = (rnorm*cos(true_anomaly)*x + rnorm*sin(true_anomaly)*y)';
u_r = r/norm(r);
u_per = (cross(z, u_r)/norm(cross(z, u_r)))';
v = (mu/hnorm)*enorm*sin(true_anomaly)*u_r + (mu/hnorm)*(1+enorm*cos(true_anomaly))*u_per;

Y0 = [r; v];%State Vector
%End Orbital Elements to State Vector

%% Numerical Integration of cartesian state vector
tf = 1*6*3600; %[s]
TSPAN = [0 tf];%duration of integration in seconds
OPTIONS = odeset('Maxstep', 10);
[TOUT, YOUT] = ode45(@relacc_with_j2, TSPAN, Y0, OPTIONS);

%Plot Results
figure(1)
plot3(YOUT(:,1), YOUT(:,2), YOUT(:,3),'r')
grid on
hold on
axis equal
%Plot Axes and Earth
plot3(linspace(0,2*a), zeros(100), zeros(100),'k','linewidth',2)
plot3(zeros(100), linspace(0,2*a), zeros(100),'k','linewidth',2)
plot3(zeros(100), zeros(100), linspace(0,2*a),'k','linewidth',2)
RE = 6378;
[X,Y,Z] = sphere;
s = surface(RE*X,RE*Y,RE*Z);
xlabel('x [km]');ylabel('y [km]');zlabel('z [km]');


%% Orbital Elements, latitude, and longitude
%(To be corrected by Earth's Rotation)
[size1,size2] = size(YOUT);
for i = 1:size1
    %State Vector to Orbital Elements
    [a_out(i), enorm(i), inclination(i), RAAN(i), arg_per(i), true_anomaly(i)] = SV2OE(YOUT(i,1:3), YOUT(i,4:6));
    lat(i) = asin(YOUT(i,3)/norm(YOUT(i,1:3)));
    lon(i) = atan2(YOUT(i,2),YOUT(i,1));
end

figure(2)
plot(TOUT/3600, a_out)
grid on
xlabel('time [hr]')
ylabel('Semi_major Axis, a [km]')

figure(3)
plot(TOUT/3600, enorm)
grid on
xlabel('time [hr]')
ylabel('Eccentricity, e')

figure(4)
plot(TOUT/3600, inclination*180/pi)
grid on
xlabel('time [hr]')
ylabel('Inclination, i [degrees]')

figure(5)
plot(TOUT/3600, RAAN*180/pi)
grid on
xlabel('time [hr]')
ylabel('RAAN, \Omega [degrees]')

figure(6)
plot(TOUT/3600, arg_per*180/pi)
grid on
xlabel('time [hr]')
ylabel('Argument of Perigee, \omega [degrees]')

%Correct Values of Longitude
correct_lon = mod((lon' -(gst+omega_E*TOUT)),2*pi); %correcting values of longitude for earth's rotation
%Plot corrected values of latitude and longitude
figure(7)
plot(correct_lon*180/pi, lat*180/pi, 'LineStyle', 'none', 'Marker', '.')
hold on
plot(correct_lon(1)*180/pi, lat(1)*180/pi, 'r*')
plot(correct_lon(end)*180/pi, lat(end)*180/pi, 'k*')
xlabel('East Longitude [deg]')
ylabel('Latitude [deg]')
