% preconditions

% Correct number of entries:
eta = ones(4, 1); Ntot = length(eta);
z_k = 0;
dr = 0.01;
newCPoints = 0;

% Check that output has correct dimensions
assert(norm(size(JacobianCurvature(eta, z_k, dr, newCPoints)) - ...
    [2 * Ntot - newCPoints + 2, 2 * Ntot - newCPoints + 2]) < 1e-12);

%% Test 1: Linear Eta. No contact points.
eta = 1:5;
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 0;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);

%% Test 2: Constant Eta. No contact points.
eta = ones(5, 1);
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 0;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);

%% Test 3: Linear Eta. One contact point.
eta = 1:5;
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 1;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);

%% Test 4: Constant Eta. One contact points.
eta = ones(5, 1);
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 1;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);

%% Test 5: Linear Eta. Two contact points.
eta = 1:5;
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 2;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);

%% Test 6: Constant Eta. Two contact points.
eta = ones(5, 1);
z_k = eta(1) + rand()/5;
dr = 0.01 + rand();
newCPoints = 2;
assert(norm(numericalJacobian(eta, z_k, dr, newCPoints, 1e-7) - ...
    JacobianCurvature(eta, z_k, dr, newCPoints)) < 1e-6, ...
    "Error with dr = %6.6g and zk = 6.6g", dr, z_k);


%% Auxiliary function
function numJacobian = numericalJacobian(eta, z_k, dr, newCPoints, epsilon)
    % This function will caculate the numerical jacobianof the curvature
    % function
    if size(eta, 2) > 1
        eta = eta';
    end
    Ntot = length(eta);
    number_of_variables = length(eta) - newCPoints;
    step_ahead_matrix = zeros(2 * Ntot - newCPoints + 2, number_of_variables);
    step_back_matrix  = zeros(2 * Ntot - newCPoints + 2, number_of_variables);
    sphere_step_ahead = curvature(eta, z_k + epsilon, dr, newCPoints);
    sphere_step_back  = curvature(eta, z_k - epsilon, dr, newCPoints);
    for ii = 1:number_of_variables
        step_vector = zeros(Ntot, 1);
        step_vector(ii + newCPoints) = epsilon;
        step_ahead_matrix(:, ii) = curvature(eta + step_vector, z_k, dr, newCPoints);
        step_back_matrix(:, ii)  = curvature(eta - step_vector, z_k, dr, newCPoints);
    end
    numJacobian = [(step_ahead_matrix - step_back_matrix)/(2*epsilon), ...
        zeros(2 * Ntot - newCPoints + 2, Ntot), ...
        (sphere_step_ahead - sphere_step_back)/(2*epsilon), ...
        zeros(2 * Ntot - newCPoints + 2, 1)];
end