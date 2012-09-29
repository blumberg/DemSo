figure

subplot(1,2,2);
plot(M(:,1), M(:,7), '.', 'MarkerSize', 5);
%axis([0 7 -6 4]);
title('Vy: Particle velocity along y axis');

for i = 1:size(M,1)
    subplot(1,2,1);
    plot(M(i,4),M(i,5), 'o');
    axis([0 10 0 10]);
    title('Particle movement replay');
    % live Vy plot
%     subplot(1,2,2);
%     hold on;
%     plot(M(i,1), M(i,7));
%     axis([0 7 -6 4]);
    pause(0.005);
end