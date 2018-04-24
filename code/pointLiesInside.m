function [ isInside ] = pointLiesInside( point, contour )
%POINTLIESINSIDE Determine if a point lies inside a contour
%
%   Author: Andrew Cox
%   Version: Aug 25, 2015

if size(contour,1) == 2
    contour = contour';
end

quad1 = 0;
quad2 = 0;
quad3 = 0;
quad4 = 0;

% For Debugging
% h_fig = figure(); hold on;
% plot(contour(:,1), contour(:,2), 's');
% plot(point(1,1), point(1,2), 'r*');
% hold off; grid on; axis equal;

isInside = false;
for i = 1:size(contour,1)
    dx = contour(i,1) - point(1,1);
    dy = contour(i,2) - point(1,2);
    
    if (dx > 0) && (dy > 0)
        quad1 = quad1 + 1;
    end
    if (dx < 0) && (dy > 0)
        quad2 = quad2 + 1;
    end
    if (dx < 0) && (dy < 0)
        quad3 = quad3 + 1;
    end
    if (dx > 0) && (dy < 0)
        quad4 = quad4 + 1;
    end
    
    if(quad1 > 0 && quad2 > 0 && quad3 > 0 && quad4 > 0)
        isInside = true;
        break;
    end
end

% For Debugging
% if(isInside)
%     figure(h_fig); hold on;
%     plot(point(1,1), point(1,2), 'g*', 'MarkerSize', 15);
% end
end

