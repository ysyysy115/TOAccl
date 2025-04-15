% Function to compute element stiffness and mass matrices
function [Ke, Me] = element_matrices(E, nu, rho, dx, dy, thickness)
    % Plane stress constitutive matrix D
    D = (E / (1 - nu^2)) * [1, nu, 0;
                            nu, 1, 0;
                            0, 0, (1 - nu)/2];
    % Gauss points and weights for 2x2 Gauss quadrature
    gp = [-1/sqrt(3), 1/sqrt(3)];
    gw = [1, 1];

    Ke = zeros(8,8);
    Me = zeros(8,8);

    % Coordinates of element nodes (local coordinates)
    x_nodes = [0, dx, dx, 0];
    y_nodes = [0, 0, dy, dy];

    % Loop over Gauss points
    for i = 1:2
        xi = gp(i);
        wi = gw(i);
        for j = 1:2
            eta = gp(j);
            wj = gw(j);
            w = wi * wj;

            % Shape functions N
            N = 1/4 * [(1 - xi)*(1 - eta);
                       (1 + xi)*(1 - eta);
                       (1 + xi)*(1 + eta);
                       (1 - xi)*(1 + eta)];

            % Derivatives of shape functions w.r.t xi and eta
            dN_dxi = 1/4 * [-(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta)];
            dN_deta = 1/4 * [-(1 - xi), -(1 + xi),  (1 + xi),   (1 - xi)];

            % Jacobian matrix
            J = [dN_dxi; dN_deta] * [x_nodes', y_nodes'];

            detJ = det(J);

            if detJ <= 0
                error('Jacobian determinant is non-positive.');
            end

            % Inverse Jacobian
            invJ = inv(J);

            % Derivatives of shape functions w.r.t x and y
            dN_dx_dy = invJ * [dN_dxi; dN_deta];

            % B matrix (strain-displacement matrix)
            B = zeros(3,8);
            for node = 1:4
                B(1, 2*node-1) = dN_dx_dy(1, node);
                B(2, 2*node)   = dN_dx_dy(2, node);
                B(3, 2*node-1) = dN_dx_dy(2, node);
                B(3, 2*node)   = dN_dx_dy(1, node);
            end

            % Element stiffness matrix
            Ke = Ke + B' * D * B * detJ * w * thickness;

            % Element mass matrix (consistent mass matrix)
            N_mat = zeros(2,8);
            for node = 1:4
                N_mat(1, 2*node-1) = N(node);
                N_mat(2, 2*node)   = N(node);
            end
            Me = Me + (N_mat' * N_mat) * rho * detJ * w * thickness;
        end
    end
end