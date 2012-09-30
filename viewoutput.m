M = csvread('output.csv');

% Variables in X
figure(1)
subplot(3,3,1);
plot(M(:,1), M(:,4), '.', 'MarkerSize', 5);
title('X: Particle position along x axis');

subplot(3,3,2);
plot(M(:,1), M(:,6), '.', 'MarkerSize', 5);
title('Vx: Particle velocity along x axis');

subplot(3,3,3);
plot(M(:,1), M(:,8), '.', 'MarkerSize', 5);
title('Ax: Particle acceleration along x axis');

% Variables in Y
subplot(3,3,4);
plot(M(:,1), M(:,5), '.', 'MarkerSize', 5);
title('Y: Particle position along y axis');

subplot(3,3,5);
plot(M(:,1), M(:,7), '.', 'MarkerSize', 5);
title('Vy: Particle velocity along y axis');

subplot(3,3,6);
plot(M(:,1), M(:,9), '.', 'MarkerSize', 5);
title('Ay: Particle acceleration along y axis');

% Angular variables
subplot(3,3,7)
plot(M(:,1), M(:,10), '.', 'MarkerSize', 5);
title('\theta: Particle angular position');

subplot(3,3,8)
plot(M(:,1), M(:,11), '.', 'MarkerSize', 5);
title('\omega: Particle angular velocity');

subplot(3,3,9)
plot(M(:,1), M(:,12), '.', 'MarkerSize', 5);
title('\alpha: Particle angular acceleration');

% Live replay of particle movement
for i = 1:size(M,1)
    figure(2)
    plot(M(i,4),M(i,5), 'o');
    axis([0 10 0 10]);
    title('Particle movement replay');
    pause(0.005);
end
