function [dll,matH,coord_x,ir,irr,ird] = member(dummy)
%
dx = 1.0;
dy = 1.0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coord_x (nk, 2)
% --->
coord_x = [...
    0,  dy;
    dx, dy;
    0,  0];
nk = size(coord_x,1);

%%%% irr(nm, 2)
irr = [...
    1, 2;
    3, 2];
nm = size(irr,1);

%%%% ird(nk, 2)
ird = [...
    0, 0;
    1, 1;
    0, 0];
nd = sum(sum(ird));
%
ii = 0;
for j=1:nk
    for k=1:2
        if ird(j,k) == 1
            ii = ii + 1;
            ird(j,k) = ii;
        else
            ird(j,k) = nd +1;
        end
    end
end
%
%%%% ir(nm, 4)
ir = zeros(nm,4);
for i=1:nm
    for j=1:2
        ir(i,j)   = ird(irr(i,1), j);
        ir(i,j+2) = ird(irr(i,2), j);
    end
end
%
%%%% matH(nd, nm)
dll = zeros(nm,1);
matH = sparse(zeros(nd+1,nm));
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dll(i) = norm([dx; dy], 2);
    dir_cos(1) =-dx/dll(i);
    dir_cos(2) =-dy/dll(i);
    dir_cos(3) = dx/dll(i);
    dir_cos(4) = dy/dll(i);
    for j=1:4
        if abs(dir_cos(j)) < 10^(-16)
            dir_cos(j) = 0;
        end
        matH(ir(i,j),i) = dir_cos(j);
    end
end
%
matH = matH(1:nd,:);
