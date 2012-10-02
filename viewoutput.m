M = csvread('output.csv');

% Number of particles followed
N = (size(M,2)-1)/11;

% Loop over followed particles to show output
for i = 0:N-1
	% Variables in X
	figure('Name', ['Simulation output for Particle ' num2str(M(1,2+i*11))], 'NumberTitle', 'off');
	subplot(3,3,1);
	plot(M(:,1), M(:,4+i*11), '.', 'MarkerSize', 5);
	title('X: Particle position along x axis');

	subplot(3,3,2);
	plot(M(:,1), M(:,6+i*11), '.', 'MarkerSize', 5);
	title('Vx: Particle velocity along x axis');

	subplot(3,3,3);
	plot(M(:,1), M(:,8+i*11), '.', 'MarkerSize', 5);
	title('Ax: Particle acceleration along x axis');

	% Variables in Y
	subplot(3,3,4);
	plot(M(:,1), M(:,5+i*11), '.', 'MarkerSize', 5);
	title('Y: Particle position along y axis');

	subplot(3,3,5);
	plot(M(:,1), M(:,7+i*11), '.', 'MarkerSize', 5);
	title('Vy: Particle velocity along y axis');

	subplot(3,3,6);
	plot(M(:,1), M(:,9+i*11), '.', 'MarkerSize', 5);
	title('Ay: Particle acceleration along y axis');

	% Angular variables
	subplot(3,3,7)
	plot(M(:,1), M(:,10+i*11), '.', 'MarkerSize', 5);
	title('\theta: Particle angular position');

	subplot(3,3,8)
	plot(M(:,1), M(:,11+i*11), '.', 'MarkerSize', 5);
	title('\omega: Particle angular velocity');

	subplot(3,3,9)
	plot(M(:,1), M(:,12+i*11), '.', 'MarkerSize', 5);
	title('\alpha: Particle angular acceleration');
end

% Live replay of particle movement
h = figure('Name', 'Particle movement replay', 'NumberTitle', 'off');
for i = 1:size(M,1)
    set(0, 'CurrentFigure', h)
    for j = 0:N-1
        plot(M(i,4+j*11), M(i,5+j*11), 'o');
		text(M(i,4+j*11), M(i,5+j*11), num2str(M(i,2+j*11)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
		hold on
    end
    axis([0 10 0 10]);
    title('Particle movement replay');
    drawnow;
    hold off
end
