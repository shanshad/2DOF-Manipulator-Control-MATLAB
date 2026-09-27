Xd = Desired.Data(:,1);
Yd = Desired.Data(:,2);
Xa = Actual.Data(:,1);
Ya = Actual.Data(:,2);

% v = VideoWriter('introducing_disturbance.avi');
% v.FrameRate = 30; % You may lower this if desired
% open(v);

figure;
subplot(1,2,1);
h1 = animatedline('Color','r','LineWidth',2);
title('Desired Trajectory');
xlabel('X'); ylabel('Y'); axis equal; grid on;
xlim([min([Xd;Xa])-0.5, max([Xd;Xa])+0.5]);
ylim([min([Yd;Ya])-0.5, max([Yd;Ya])+0.5]);

subplot(1,2,2);
h2 = animatedline('Color','b','LineWidth',2);
title('Actual Trajectory');
xlabel('X'); ylabel('Y'); axis equal; grid on;
xlim([min([Xd;Xa])-0.5, max([Xd;Xa])+0.5]);
ylim([min([Yd;Ya])-0.5, max([Yd;Ya])+0.5]);

N = min(length(Xd), length(Xa));
skip = 10; % Record every 10th sample for faster video

for k = 1:skip:N
    addpoints(h1, Xd(k), Yd(k));
    addpoints(h2, Xa(k), Ya(k));
    drawnow;
    % frame = getframe(gcf);
    % writeVideo(v, frame);
end

% close(v);
% disp('Animation saved');
