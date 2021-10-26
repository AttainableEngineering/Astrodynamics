function dy = GVE(t, orb_el)
%implements Gauss Variational Equations on osculating orbital elements

%constants of the problem START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 6.6742*10^(-20); %km^3/(kg*s^2) gravitational constant
m1 = 5.9726*10^24; %kg mass of Earth
mu = G*m1; %km^3/s^2 gravitational parameter
J2 = 0.00108263; % [-] second zonal harmonic
Re = 6378;%km radius of Earth
%constants of the problem END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input orbital elements
a = orb_el(1);
enorm = orb_el(2);
inclination = orb_el(3); %rad
RAAN = orb_el(4); %rad
arg_per = orb_el(5); %rad
true_anomaly = orb_el(6); %rad

%%%%%FROM ORBITAL ELEMENTS TO STATE VECTOR (START)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hnorm = sqrt(a*mu*(1-enorm^2));
T = (2*pi/sqrt(mu))*a^(3/2); %s period

R1_i = [1 0 0; 0 cos(inclination) sin(inclination); 0 -sin(inclination) cos(inclination)];
R3_Om = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
R3_om = [cos(arg_per) sin(arg_per) 0; -sin(arg_per) cos(arg_per) 0; 0 0 1];

support_var = R3_om*R1_i*R3_Om;
x = support_var(1,1:3);
y = support_var(2,1:3);
z = support_var(3,1:3);

e_orb = enorm*x;
h_orb = hnorm*z;

rnorm = ((hnorm^2)/mu)*(1/(1+enorm*cos(true_anomaly)));
r = (rnorm*cos(true_anomaly)*x + rnorm*sin(true_anomaly)*y)';
u_r = r/norm(r);
u_per = (cross(z, u_r)/norm(cross(z, u_r)))';
v = (mu/hnorm)*enorm*sin(true_anomaly)*u_r + (mu/hnorm)*(1+enorm*cos(true_anomaly))*u_per;
p = a*(1-enorm^2);
%%%%%FROM ORBITAL ELEMENTS TO STATE VECTOR (END)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%perturbation of interest (pr, ps, pw, in LVLH)

% %%J2
pr = (-3/2)*((J2*mu*Re^2)/(rnorm^4))*(1-3*(sin(inclination)^2*sin(true_anomaly+arg_per)^2));
ps = (-3/2)*((J2*mu*Re^2)/(rnorm^4))*((sin(inclination)^2*sin(2*(true_anomaly+arg_per))));
pw = (-3/2)*((J2*mu*Re^2)/(rnorm^4))*((sin(2*inclination)*sin(true_anomaly+arg_per)));

%continuous thrust from uncollected HW 6
% 1) .0008 km/s^2 always pushing in the positive radial direction
% 2) .000005 km/s^2 always pushing in the positive velocity direction
% 3) .002 km/s^2 always pushing in the positive h (angular momentum) direction

%Comment these out for Problem 1
 pr = pr + .0008;
 ps = ps + .000005;
 pw = pw + .002;

%GVEs
a_dot            = 2*(a^2*enorm*sin(true_anomaly)/hnorm)*pr + 2*(p*a^2/(hnorm*rnorm))*ps;
e_dot            = (hnorm/mu)*sin(true_anomaly)*pr + (1/(hnorm*mu))*((hnorm^2+mu*rnorm)*cos(true_anomaly)+mu*enorm*rnorm)*ps;
inclination_dot  = (rnorm/hnorm)*cos(true_anomaly+arg_per)*pw;
RAAN_dot         = (rnorm/(hnorm*sin(inclination)))*sin(true_anomaly+arg_per)*pw;
arg_per_dot      = -(1/(hnorm*enorm))*((hnorm^2/mu)*cos(true_anomaly)*pr - (rnorm+hnorm^2/mu)*sin(true_anomaly)*ps) - ...
                    ((rnorm*sin(true_anomaly+arg_per))/(hnorm*tan(inclination)))*pw;
true_anomaly_dot = hnorm/rnorm^2 + (1/(hnorm*enorm))*((hnorm^2/mu)*cos(true_anomaly)*pr - (rnorm+hnorm^2/mu)*sin(true_anomaly)*ps);

%output
dy = [a_dot; e_dot; inclination_dot; RAAN_dot; arg_per_dot; true_anomaly_dot];

end

