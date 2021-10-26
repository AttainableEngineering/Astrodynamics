function dy_dt = CWeqq_HW3(t, u)
%Clohessy-Wiltshire differential equations of S/C relative motion
global theta_dot_c

om = theta_dot_c;

x = u(1);
y = u(2);
z = u(3);
xdot = u(4);
ydot = u(5);
zdot = u(6);

xddot = 2*om*ydot + 3*(om^2)*x;
yddot = -2*om*xdot;
zddot = -(om^2)*z;

dy_dt = [xdot; ydot; zdot; xddot; yddot; zddot];

end

