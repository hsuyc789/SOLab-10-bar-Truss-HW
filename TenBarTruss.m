function [Q,stress,R,K] = TenBarTruss(r)
nod_coor=[18.28 9.14; 18.28 0; 9.14 9.14; 9.14 0; 0 9.14; 0 0]; % node1~6 x、y座標
ele_con=[3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4]; % elements' node connection

L(1:6) = 9.14;
L(7:10) = 12.9259;
A(1:6) = pi*r(1)^2;
A(7:10) = pi*r(2)^2;

ele_dof=[5 6 9 10; 1 2 5 6; 7 8 11 12; 3 4 7 8; 5 6 7 8; 1 2 3 4; ...
        7 8 9 10; 5 6 11 12; 3 4 5 6;  1 2 7 8]; % elements' node dof(from node1 to node10)組成10x4矩陣
E = 200*10^9;
K = zeros(12);
for e=1:10 %針對10個元素進行運算
    C=(nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))/L(e); %計算各元素cos值
    S=(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))/L(e); %計算各元素sin值
    k=(A(e)*E/L(e)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;-C*C -C*S C*C C*S; ...
        -C*S -S*S C*S S*S]); %建立各元素的子矩陣

    ele_dof_vec=ele_dof(e,:); % 賦予第e列的ele_dof矩陣給列矩陣ele_dof_vec(:全部元素)
    for i=1:4
       for j=1:4
            K( ele_dof_vec(1,i), ele_dof_vec(1,j) ) = ...
            K( ele_dof_vec(1,i), ele_dof_vec(1,j) ) + k(i,j);
       end
    end
end

F(4)=-10^7;
F(8)=-10^7;
Q = inv(K(1:8,1:8))*F';
Q(12)=0;

for e=1:10
    C=(nod_coor(ele_con(e,2),1)-nod_coor(ele_con(e,1),1))/L(e);
    S=(nod_coor(ele_con(e,2),2)-nod_coor(ele_con(e,1),2))/L(e);
    stress(e)=(E/L(e))*[-C -S C S]*Q( (ele_dof(e,:)) );
end

R = K(9:12,1:12)*Q;

