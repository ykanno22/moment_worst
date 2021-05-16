function [dll,matH,coord_x,ir,irr,ird,Idx] = member(nx, ny, max_length)
%
dx = 1.0;
dy = 1.0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coord_x (nk, 2)
% --->
for i=1:(ny+1)
    for j=1:(nx++1)
        pp = ((i-1) * (nx+1)) + j;
        coord_x(pp,1) = (j-1) * dx;
        coord_x(pp,2) = (i-1) * dy;
    end
end
coord_x =  1.0 * coord_x;
nk = size(coord_x,1);

%%%% irr(nm, 2)

filename_ird = strcat(num2str(nx),'_',num2str(ny),'_irr.dat');
fid = fopen(filename_ird,'r');
if fid ~= -1
    irr = fscanf(fid,'%i %i',[2 inf]);
    irr = irr';
    fclose(fid);
    iFlag_create_irr = 0;
else
    iFlag_create_irr = 1;
end

if iFlag_create_irr == 1
    irr = [];
    for i=1:(nk-1)
        for j=(i+1):nk
            irr = [irr; i, j];
        end
    end
end
nm = size(irr,1);
%
%%%% ird(nk, 2)
ird = ones(nk,2);
ird(1,:) = zeros(1,2);
ird(((nx+1) * ny) + 1,:) = zeros(1,2);
%
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

if iFlag_create_irr == 1
    Idx.long_member = find(dll > max_length);
    %
    Idx.not_long_member = setdiff([1:nm]', Idx.long_member);
    %
    irr  = irr(Idx.not_long_member,:);
    ir   = ir(Idx.not_long_member,:);
    dll  = dll(Idx.not_long_member,:);
    matH = matH(:,Idx.not_long_member);
    nm = size(irr,1);
    
    Ind.overlap = [];
    vec_xi = cell(2,1);
    vec_xj = cell(2,1);
    for i=1:(nm-1)
        for j=(i+1):nm
            iFlag_overlap = 0;
            if dll(i) > dll(j)
                i1 = i;
                i2 = j;
            else
                i1 = j;
                i2 = i;
            end
            if isempty(intersect(i1,Ind.overlap))
                for jj=1:2
                    vec_xi1{jj} = coord_x(irr(i1,jj),:)';
                    vec_xi2{jj} = coord_x(irr(i2,jj),:)';
                end
                vec_ai1 = vec_xi1{2} - vec_xi1{1};
                vec_coef = [-vec_ai1(2); vec_ai1(1)];
                for jj=1:2
                    res_line(jj) = abs(vec_coef' * (vec_xi2{jj} - vec_xi1{1}));
                end
                if sum(res_line) < 10^(-7)
                    for jj=1:2
                        if abs(vec_ai1(1)) > 10^(-4)
                            var_t(jj) = (vec_xi2{jj}(1) - vec_xi1{1}(1)) / vec_ai1(1);
                        else
                            var_t(jj) = (vec_xi2{jj}(2) - vec_xi1{1}(2)) / vec_ai1(2);
                        end
                    end
                    if (var_t(1) >= -10^(-8)) && (var_t(1) <= 1+10^(-8))
                        if (var_t(2) >= -10^(-8)) && (var_t(2) <= 1+10^(-8))
                            iFlag_overlap = 1;
                        end
                    end
                end
            end
            %
            if iFlag_overlap == 1
                Ind.overlap = [Ind.overlap; i1];
                iFlag_overlap = 0;
            end
        end
    end
    
    Ind.notoverlap = setdiff([1:nm]', Ind.overlap);
    
    irr  = irr(Ind.notoverlap,:);
    ir   = ir(Ind.notoverlap,:);
    dll  = dll(Ind.notoverlap,:);
    matH = matH(:,Ind.notoverlap);
    nm = size(irr,1);
    
    fid = fopen(filename_ird, 'w');
    for i=1:nm
        fprintf(fid, ' %g %g \n', irr(i,1), irr(i,2));
    end
    fclose(fid);
end




Idx.memb_node = cell(nk,1);
for j=1:nk
    [jj,jjj] = find(irr==j);
    Idx.memb_node{j} = sort(jj)';
end

Idx.free_node = [];
for j=1:nk
    if max(ird(j,:)) <= nd
        Idx.free_node = [Idx.free_node, j];
    end
end

