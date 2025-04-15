function [frequencies, dfreqfx_full,V] = FEAandSens_start(x,Lx,Ly)
% profile on

    % Material properties
    E = 200000;      % Young's modulus in MPa
    nu = 0.3;        % Poisson's ratio
    rho = 8e-9;      % Density in kg/mm^3
    thickness = 1;   % Thickness in mm

    % Penalization factors for stiffness and mass matrices
    p = 3;  % Penalization for stiffness (SIMP)
    q = 6;  % Penalization for mass

    % Read the size of the design density matrix
    [nely, nelx] = size(x);

    % Element dimensions
    dx = Lx/nelx; % Element size in x (mm)
    dy = Ly/nely; % Element size in y (mm)

    % Number of nodes and elements
    nndx = nelx + 1;
    nndy = nely + 1;
    nnd = nndx * nndy; 
    nel = nelx * nely; 
    nne = 4; 
    nodof = 2; 
    eldof = nne * nodof;

    % Generate Node Map and store node_row and node_col
    NODE_MAP = zeros(nndy, nndx);
    node_row = zeros(nnd, 1);
    node_col = zeros(nnd, 1);
    inode = 0;
    for col = 1:nndx
        for row = 1:nndy
            inode = inode + 1;
            NODE_MAP(row, col) = inode;
            node_row(inode) = row;
            node_col(inode) = col;
        end
    end
    NODE_MAP = flipud(NODE_MAP); % Flip up-down to start numbering from bottom-left
    % Adjust node_row accordingly
    node_row = nndy - node_row + 1;

    % Node coordinates
    geom = zeros(nnd, 2);
    for inode = 1:nnd
        col = node_col(inode);
        row = node_row(inode);
        geom(inode, 1) = (col - 1) * dx; % x coordinate
        geom(inode, 2) = (row - 1) * dy; % y coordinate
    end

    % Generate Element Map and store element_row and element_col
    ELE_MAP = zeros(nely, nelx);
    element_row = zeros(nel, 1);
    element_col = zeros(nel, 1);
    iel = 0;
    for col = 1:nelx
        for row = 1:nely
            iel = iel + 1;
            ELE_MAP(row, col) = iel;
            element_row(iel) = row;
            element_col(iel) = col;
        end
    end
    ELE_MAP = flipud(ELE_MAP); % Flip up-down
    % Adjust element_row accordingly
    element_row = nely - element_row + 1;

    % Generate Element Connectivity
    elements = zeros(nel, 4);
    for iel = 1:nel
        hang = element_row(iel);
        lie = element_col(iel);
        elements(iel, 1) = NODE_MAP(hang + 1, lie);     % Node 1 (lower-left)
        elements(iel, 2) = NODE_MAP(hang + 1, lie + 1); % Node 2 (lower-right)
        elements(iel, 3) = NODE_MAP(hang, lie + 1);     % Node 3 (upper-right)
        elements(iel, 4) = NODE_MAP(hang, lie);         % Node 4 (upper-left)
    end
    ELE_COOR = elements;

    % Total Degrees of Freedom (DOF)
    ndof = nnd * nodof;

    % Compute element stiffness and mass matrices
    [Ke_s, Me_s] = element_matrices(E, nu, rho, dx, dy, thickness);

    % Preallocate arrays for sparse assembly
    eldof_sq = eldof^2;
    indexi = zeros(eldof_sq * nel, 1);
    indexj = zeros(eldof_sq * nel, 1);
    dataK = zeros(eldof_sq * nel, 1);
    dataM = zeros(eldof_sq * nel, 1);

    idx = 0;

    for e = 1:nel
        % Nodes of the element
        nodes = elements(e, :);
        % Global DOF indices
        edof = zeros(eldof, 1);
        for k = 1:4
            edof(2 * k - 1) = 2 * nodes(k) - 1;
            edof(2 * k) = 2 * nodes(k);
        end

        % Mapping element number to design variable indices
        row = element_row(e);
        col = element_col(e);
        x_e = x(row, col);

        % Penalization for stiffness (SIMP method)
        Ke = x_e^p * Ke_s;

        % Penalization for mass
        if x_e > 0.1
            Me = x_e * Me_s;
        else
            Me = x_e^q * Me_s;
        end

        % Assemble K and M
        [row_idx, col_idx] = meshgrid(edof, edof);
        idx_range = idx + 1:idx + eldof_sq;
        indexi(idx_range) = row_idx(:);
        indexj(idx_range) = col_idx(:);
        dataK(idx_range) = Ke(:);
        dataM(idx_range) = Me(:);
        idx = idx + eldof_sq;
    end

    % Assemble global stiffness and mass matrices
    K = sparse(indexi, indexj, dataK, ndof, ndof);
    M = sparse(indexi, indexj, dataM, ndof, ndof);

    % Apply boundary conditions
    left_nodes = NODE_MAP(:, 1);
    right_nodes = NODE_MAP(:, end);

    % Initialize nf matrix
    nf = ones(nnd, nodof);
    % Fix DOFs on left edge (both DOFs)
    nf(left_nodes, :) = 0;
    % Fix DOFs on right edge (only x DOF)
    nf(right_nodes, 1) = 0; % x DOF fixed
    % % Fix DOFs on right edge (only x DOF)
    % nf(right_nodes, :) = 0; % x DOF fixed

    % Counting free DOFs
    n = 0;
    id = zeros(nnd, nodof);
    free_dofs = [];
    fixed_dofs = [];
    for i = 1:nnd
        for j = 1:nodof
            if nf(i, j) == 0
                dof = nodof * (i - 1) + j;
                fixed_dofs = [fixed_dofs; dof];
            else
                n = n + 1;
                nf(i, j) = n;
                dof = nodof * (i - 1) + j;
                free_dofs = [free_dofs; dof];
            end
        end
    end

    % Extract submatrices
    K_ff = K(free_dofs, free_dofs);
    M_ff = M(free_dofs, free_dofs);

    % Number of eigenvalues to compute
    num_eigenvalues = 9;

    % Solve generalized eigenvalue problem
    opts.maxit = 1000;
    opts.tol = 1e-6;
    [V, omega_sq] = eigs(K_ff, M_ff, num_eigenvalues, 'sm', opts);

    % Eigenfrequencies in Hz
    omega = sqrt(diag(omega_sq));
    frequencies = omega / (2 * pi);

    % Sensitivity analysis
    d_omega_dx = zeros(nely, nelx, num_eigenvalues);

    % Reconstruct full eigenvectors
    TRUE_EIGVEC = zeros(ndof, num_eigenvalues);
    TRUE_EIGVEC(free_dofs, :) = -V;

    %% Sensitivity Analysis with Parallel Computing
    % Preallocate memory for results
    dfreqfx_full = zeros(nely, nelx, num_eigenvalues);

    % Precompute element degrees of freedom (dofs)
    element_dofs = [...
        ELE_COOR(:,1)*2-1, ELE_COOR(:,1)*2, ...
        ELE_COOR(:,2)*2-1, ELE_COOR(:,2)*2, ...
        ELE_COOR(:,3)*2-1, ELE_COOR(:,3)*2, ...
        ELE_COOR(:,4)*2-1, ELE_COOR(:,4)*2];

    for mode = 1:num_eigenvalues
        % Extract eigenvector for current mode
        eigvec_mode = TRUE_EIGVEC(:, mode);

        % Assemble u for all elements (vectorized)
        u = eigvec_mode(element_dofs)';  % u is 8 x nel

        % Sensitivity computation
        dfreqfx1 = zeros(nel, 1);

        % Start parallel loop over elements
        for iel = 1:nel
            row = element_row(iel);
            col = element_col(iel);
            x_e = x(row, col);

            u_e = u(:, iel);

            if x_e > 0.1
                dfreqfx1(iel) = u_e' * (p * x_e^(p-1) * Ke_s - frequencies(mode) * Me_s) * u_e;
            else
                dfreqfx1(iel) = u_e' * (p * x_e^(p-1) * Ke_s - frequencies(mode) * q * x_e^(q-1) * Me_s) * u_e;
            end
        end

        % Normalize and reshape dfreqfx1
        dfreqfx1 = dfreqfx1 / max(abs(dfreqfx1));  % Normalize
        dfreqfx2 = flipud(reshape(dfreqfx1, [nely, nelx]));  % Reshape and flip vertically

        % Store sensitivity for this mode
        dfreqfx_full(:, :, mode) = dfreqfx2;
    end
    % profile viewer
end


