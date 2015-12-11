
% Load in positions and paramters.
pos   = load('planetposition_verlet.txt');
%pos   = load('Rungekuttaposition.txt');
%nfile = load('n.dat');

% Parameters from the n-file.
n    = 400;%number of planets
time = 4;
R0   = 20;
timeSteps =4000; 

% Test cases.
% x = pos(:,4);
% y = pos(:,5);
% z = pos(:,6);
% x2 = pos(:,3);
% y2 = pos(:,4);
% x3 = pos(:,5);
% y3 = pos(:,6);

% Actual x, y, and z-data.
startx = zeros(n,timeSteps);
starty = zeros(n,timeSteps);
startz = zeros(n,timeSteps);

% Actual x, y, and z-data.
startx2 = zeros(1,timeSteps);
starty2 = zeros(1,timeSteps);
startz2 = zeros(1,timeSteps);

for j = 1:timeSteps
    for i = 0:n-1
        startx(i+1,j) = pos(j,1+i*3);
        starty(i+1,j) = pos(j,2+i*3);
        startz(i+1,j) = pos(j,3+i*3);
    end
    startx2(1,j) = pos(j,1+n*3);
    starty2(1,j) = pos(j,2+n*3);
    startz2(1,j) = pos(j,3+n*3);
end


% Set speed(ish) of the animation.
pauseL = 0.01;

% Set up and draw initial positions in the 3d plot.
figure(1);
h = plot3(startx(:,1),starty(:,1),startz(:,1), ...
    '.', 'LineWidth', 1);
hold on
h2=plot3(startx2(1),starty2(1),startz2(1), ...
    '.red', 'LineWidth', 1);
hh = title('Time, $t = 0.000 [\tau]$', 'interpreter', 'latex', 'FontSize', 14);
xlabel('x, [ly]', 'interpreter', 'latex', 'FontSize', 14);
ylabel('y, [ly]', 'interpreter', 'latex', 'FontSize', 14);
zlabel('z, [ly]', 'interpreter', 'latex', 'FontSize', 14);

gs = 10;
axis([-gs gs -gs gs -gs gs]);
grid on;

hButton = uicontrol('style', 'pushbutton', 'Callback', 'uiresume(gcbf)');
hSlide  = uicontrol('style','slider',...
                    'min',1,...
                    'max',timeSteps,...
                    'value',1,...
                    'sliderstep',[.001 0.01],...
                    'units','normalized',...
                    'position',[.02 .02 .98 .02]);
figure(1);

while true
    figure = h;

    sliderVal = round(get(hSlide,'value')); % Current value of the slider.
    j = sliderVal;
    
    set(h, ...
        'XData', startx(:,j), ...
        'YData', starty(:,j), ...
        'ZData', startz(:,j));
    set(h2, ...
        'XData', startx2(:,j), ...
        'YData', starty2(:,j), ...
        'ZData', startz2(:,j));
    
    % Update time counter in title of the plot.
    timeString = sprintf('Time, $t = %.3f \\tau$', j * (time / timeSteps));
    set(hh, 'String', timeString);
    pause(1e-10);
end


% filename = 'testnew51.gif';
%
%
% for j = 1:timeSteps
%     % Update positions in the 3d plot.
%     set(h, ...
%         'XData', startx(:,j), ...
%         'YData', starty(:,j), ...
%         'ZData', startz(:,j));
%
%     % Update time counter in title of the plot.
%     timeString = sprintf('Time, $t = %.3f \\tau$', j * (time / timeSteps));
%     set(hh, 'String', timeString);
%     pause(pauseL);
%
% %     % Create .gif from animation.... Maybe.
% %     frame = getframe(1);
% %     im = frame2im(frame);
% %     [imind,cm] = rgb2ind(im,256);
% %
% %     if j == 1;
% %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
% %     else
% %         imwrite(imind,cm,filename,'gif','WriteMode','append');
% %     end
%
% end

%%

% n = 100;
% plot(x(1:n),y(1:n), 'r-');
% hold on;
% plot(x2(1:n),y2(1:n), 'b-');
%
% plot(x3(1:n),y3(1:n), 'g-');
% axis equal


% Plot radial distribution of objects
Compute r.
r = sqrt(startx(:,end-1000).^2 + starty(:,end-1000).^2 + startz(:,end-1000).^2);

% Plot histogram of r.
figure(5);
[bins, R] = hist(r, 500);
xlabel('r, ly]');
ylabel('# of objects');




