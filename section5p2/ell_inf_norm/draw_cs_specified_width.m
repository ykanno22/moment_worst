function [dummy] = draw_cs_specified_width(coord_x,irr,af)
%
dummy = 1;
%
figure;
clf reset;
hold on;
axis equal
axis off;
%
nm = size(irr,1);
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    mem0(1,:,i) = coord_x(j1,:);
    mem0(2,:,i) = coord_x(j2,:);
end
%
% af = af / (max(af) / 10);
af = af * 1.0;

for i=1:nm
    if af(i) > 1.0*10^(-4)
        plot(mem0(:,1,i), mem0(:,2,i), 'b-', 'LineWidth',af(i));
    elseif af(i) < -10^(-4)
        plot(mem0(:,1,i), mem0(:,2,i), 'r-', 'LineWidth',-af(i));
    end
    hold on;
end

for j=1:size(coord_x,1)
    plot(coord_x(j,1), coord_x(j,2), 'ok',...
        'MarkerFaceColor', 'w', 'MarkerSize',10);
end

axis equal
set(gcf,'renderer','painters');

