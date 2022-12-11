clc; clear;

U11 = sym('U11');
U12 = sym('U12');
U21 = sym('U21');
U22 = sym('U22');
U31 = sym('U31');
U32 = sym('U32');

V11 = sym('V11');
V12 = sym('V12');
V21 = sym('V21');
V22 = sym('V22');
V31 = sym('V31');
V32 = sym('V32');

%% step size
t = sym('t');

%% initial position
U1 = [U11,U12];
U2 = [U21,U22];
U3 = [U31,U32];

%% moving direction
V1 = [V11,V12];
V2 = [V21,V22];
V3 = [V31,V32];

%% edges after moving
A = [(U2+V2*t) - (U1+ V1*t)];
B = [(U3+V3*t) - (U1+ V1*t)];
C = [A; B];

solve(det(C), t);
cf = coeffs(det(C),t); % Now cf(1),cf(2),cf(3) holds the coefficients for the polynom. at order c,b,a

Ds = [U2-U1; U3-U1];
D  = [V2-V1; V3-V1];


