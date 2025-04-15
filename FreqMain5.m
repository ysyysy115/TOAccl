clc
clear all


nelx = 500;
nely = 45;

Lx = 500;
Ly = 45;

volfrac = 0.3;
rmin = 1.5;

x(1:nely, 1:nelx) = volfrac; 
change = 1;
loop = 0;
maxloop = 300;

first_eigenvalue_history=[];

while change > 0.01 && loop < maxloop 
    loop = loop + 1;
    xold = x;

    x(1, :) = 1;
    x(end, :) = 1;

    % FE-ANALYSIS
    if loop ==1
        [frequencies, d_omega_dx,V] = FEAandSens_start(x, Lx, Ly);
    else
        startingV = V(:,1);
        [frequencies, d_omega_dx] = FEAandSens(x, Lx, Ly,startingV, loop);
    end

    
    % Save the first eigenvalue
    first_eigenvalue_history = [first_eigenvalue_history, frequencies(1)];

    % Aggregate sensitivities
    weights = [1 0.6 0.6 0.4 0.3 0.3 0.3 0.2 0.2]; % Define appropriate weights
    dc = zeros(nely, nelx);
    c = 0;
    for mode = 1:size(d_omega_dx, 3)
        dc = dc + weights(mode) * d_omega_dx(:, :, mode);
        c = c + weights(mode) * frequencies(mode);
    end
    dc = -dc;

    % FILTERING OF SENSITIVITIES
    [dc] = check(nelx, nely, rmin, x, dc);    

    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [x] = OC(nelx, nely, x, volfrac, dc); 
    vol=sum(sum(x)) / (nelx * nely);

    if vol < volfrac*0.2
        break
    end

    % PRINT RESULTS
    change = max(max(abs(x - xold)));
    disp([' It.: ' sprintf('%4i', loop) ' Obj.: ' sprintf('%10.4f', c) ...
         ' Vol.: ' sprintf('%6.3f', vol) ...
         ' ch.: ' sprintf('%6.3f', change)]);

    % PLOT DENSITIES  
    figure(2) 
    colormap(gray); imagesc(-real(x)); axis equal; axis tight; axis off; pause(1e-6);
end

% Plot history of the first eigenvalue
figure(3)
plot(1:size(first_eigenvalue_history,2), first_eigenvalue_history, '-o', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('First Eigenvalue');
title('History of the First Eigenvalue');
grid on;