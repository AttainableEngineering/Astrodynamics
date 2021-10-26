function ds = HW1_2017_P2_EOM(t, s)
mu_e = 398600; %[km^3/s^2]
r = s(1:3);
v = s(4:6);
a = -(mu_e/(norm(r)^3))*r;
ds = [v;a];