function [a, enorm, inclination, RAAN, arg_per, true_anomaly] = SV2OE_HW3(r0, v0)
%from State Vector to Orbital Elements

%constants of the problem START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 6.6742*10^(-20); %km^3/(kg*s^2) gravitational constant
m1 = 5.972*10^24; %kg mass of Earth
mu = G*m1; %km^3/s^2 gravitational parameter
%constants of the problem END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = cross(r0, v0); %km^2/s angular momentum vector
hnorm = norm(h);
arg_acos = h(3)/norm(h);
if norm(arg_acos) >= 1
    arg_acos = sign(arg_acos);
end
inclination = acos(arg_acos); %inclination
%inclination_deg = inclination*180/pi
if inclination == 0 || inclination == pi
    node_line = [1;0;0]; %picking x axis
    RAAN = 0;   
else
    node_line = cross([0;0;1],h)/norm(cross([0;0;1],h));
    arg_acos = node_line(1);
    if norm(arg_acos) >= 1
        arg_acos = sign(arg_acos);
    end
    RAAN = acos(arg_acos); %right ascension of the ascending node
end
if node_line(2) < 0
    RAAN = 2*pi - RAAN;
% else
%     RAAN
end
%RAAN_deg = RAAN*180/pi;
e = cross(v0, h)/mu - r0/norm(r0); %[-] eccentricity vector
enorm = norm(e); %printing norm of e vector
if enorm < 10^(-5)
    arg_acos = dot(node_line, r0)/norm(r0);
    if norm(arg_acos) >= 1
        arg_acos = sign(arg_acos);
    end
    arg_per = acos(arg_acos); % argument of perigee
    if r0(3) < 0
        arg_per = 2*pi - arg_per;
%     else
%         arg_per
    end
else
    arg_acos = dot(node_line, e)/enorm;
    if norm(arg_acos) >= 1
        arg_acos = sign(arg_acos);
    end
    arg_per = acos(arg_acos); % argument of perigee
    if e(3) < 0
        arg_per = 2*pi - arg_per;
%     else
%         arg_per
    end
end
%arg_per_deg = arg_per*180/pi
rp = (norm(h)^2/mu)*(1/(1+norm(e))); %km perigee
ra = (norm(h)^2/mu)*(1/(1-norm(e))); %km apogee
if ra < 0
    ra = -ra;
end
a = (rp+ra)/2; %km semi-major axis
T = (2*pi/sqrt(mu))*a^(3/2); %s period
if enorm < 10^(-5)
    true_anomaly = 0;
else
    arg_acos = dot(e, r0)/(enorm*norm(r0));
    if norm(arg_acos) >= 1
        arg_acos = sign(arg_acos);
    end
    true_anomaly = acos(arg_acos); % true anomaly
end
ur = r0/norm(r0);
if dot(v0, ur) < 0
    true_anomaly = 2*pi - true_anomaly;
% else
%     true_anomaly
end
%true_anomaly_deg = true_anomaly*180/pi

end

