% This is the unit test script to test the curvature 
ERROR_TOLERANCE = 1e-10;

% preconditions
eta = 1:25;
dr = .01; % To ensure that dr * newCPoints < 1
z_k = 0;
for ii = 0:3:15
    assert(length(curvature(eta, z_k, dr, ii)) == 2 * length(eta) - ii + 2, "Assertion failed with newCPoints = %d", ii);
end


%% Test 1.1: No contact points. Linear Eta

eta = 1:5;
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;

C = dr * sqrt(1 + 1/dr^2);
assert(norm(curvature(eta, z_k, dr, 0) - [zeros(5, 1); 4/dr^2; (1/dr)/(C); ...
    1/(2*dr*C); 1/(3*dr*C); ((-10+4)/dr^2)/(1+(2/dr)^2)^(3/2) + (-2/dr)/(dr*4*(1+(-2/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
        "Error with no contact points and linear eta. z_k = %6.6g", z_k);


%% Test 1.2: No contact points. Constant Eta.

eta = ones(1, 5);
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;

assert(norm(curvature(eta, z_k, dr, 0) - [zeros(5, 1); 0; 0; 0; 0; ...
    ((-1)/dr^2)/(1+(1/(2*dr))^2)^(3/2) + (-.5/dr)/(dr*4*(1+(-.5/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
    "Error with no contact points and constant eta. z_k = %6.6g", z_k);

%% Test 1.3: No contact points. Zig-Zag Eta.

eta = [1, 3, 1, 3, 1];
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;

assert(norm(curvature(eta, z_k, dr, 0) - [zeros(5, 1); 8/dr^2; -4/dr^2; 4/dr^2; -4/dr^2; ...
    (1/dr^2)/(1+(-3/(2*dr))^2)^(3/2) + (-1.5/dr)/(dr*4*(1+(-1.5/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
    "Error with zero contact points and zig-zag eta. z_k = %6.6g", z_k);



%% Test 1.4: One contact point. Linear Eta

eta = 1:5;
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;
Eta_before = z_k - sqrt(1 - dr^2 * (1-1)^2);
C = dr * sqrt(1 + 1/dr^2);
assert(norm(curvature(eta, z_k, dr, 1) - [zeros(4, 1); 2;...
    ((-1+Eta_before)/dr^2)/(1+((3-Eta_before)/(2*dr))^2)^(3/2)+...
    ((3-Eta_before)/(2*dr))/(1*dr*(1+((3-Eta_before)/(2*dr))^2)^(1/2)); ...
    1/(2*dr*C); 1/(3*dr*C); ...
    ((-10+4)/dr^2)/(1+(2/dr)^2)^(3/2) + (-2/dr)/(dr*4*(1+(-2/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
    "Error with one contact point and linear eta. z_k = %6.6g", z_k);


%% Test 1.5: One contact point. Zig-Zag Eta

eta = [1, 3, 1, 3, 1];
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;
Eta_before = z_k - sqrt(1 - dr^2 * (1-1)^2);
assert(norm(curvature(eta, z_k, dr, 1) - [zeros(4, 1); 2; ...
    ((-5+Eta_before)/dr^2)/(1+((1-Eta_before)/(2*dr))^2)^(3/2)+...
    ((1-Eta_before)/(2*dr))/(1*dr*(1+((1-Eta_before)/(2*dr))^2)^(1/2)); ...
    4/(dr^2); -4/(dr^2); ...
    (1/dr^2)/(1+(-1.5/dr)^2)^(3/2) + (-1.5/dr)/(dr*4*(1+(-1.5/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
    "Error with one contact point and Zig-Zag eta. z_k = %6.6g", z_k);


% Test 1.6: One contact point. Constant Eta

% eta = ones(1, 5);
% dr = 0.1;
% z_k = 0;
% Eta_before = z_k - sqrt(1 - dr^2 * (1-1)^2);
% C = dr * sqrt(1 + 1/dr^2);
% assert(norm(curvature(eta, z_k, dr, 1) - [zeros(4, 1); 2; ((3+1)/(2*dr))/(1*dr*(1+(2/dr)^2)^(1/2));...
%     1/(2*dr*C); 1/(3*dr*C); ((-10+4)/dr^2)/(1+(2/dr)^2)^(3/2) + (-2/dr)/(dr*4*(1+(-2/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
%     "Error with one contact point and linear eta");

%% Test 1.7: Two contact points. Zig-Zag Eta

eta = [1, 3, 1, 3, 1];
dr = 0.01 + rand()/10;
z_k = eta(1) - rand()/10;
Eta_before = z_k - sqrt(1 - dr^2 * 1^2);
C = dr * sqrt(1 + 1/dr^2);
assert(norm(curvature(eta, z_k, dr, 2) - [zeros(3, 1); 2; 2; ...
    ((1+Eta_before)/dr^2)/(1+((3-Eta_before)/(2*dr))^2)^(3/2)+...
    ((3-Eta_before)/(2*dr))/(2*dr*(1+((3-Eta_before)/(2*dr))^2)^(1/2)); ...
     -4/(dr^2); ...
    (1/dr^2)/(1+(-1.5/dr)^2)^(3/2) + (-1.5/dr)/(dr*4*(1+(-1.5/dr)^2)^(1/2)); 0; 0]) < ERROR_TOLERANCE, ...
    "Error with one contact point and Zig-Zag eta. z_k = %6.6g", z_k);



