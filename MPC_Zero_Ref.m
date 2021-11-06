function [M,C,Q_bar,R_bar,G,E,H,U_k]= MPC_Zero_Ref(A,B,N,x_k,Q,R,F)

n=size(A,1);  % A是 n x n 矩阵，得到n
p=size(B,2);  % B是 n x p 矩阵，得到p
M=[eye(n);zeros(N*n,n)];   %初始化M矩阵,M矩阵是(N+1)n x n的，它上面是 n x n 个 1，这里先把下半部分写成0
C=zeros((N+1)*n, N*p);
%定义M和C
tmp=eye(n);
for i=1:N
    rows = i*n+(1:n); %定义当前行数，从i x n 开始，共n行
    C(rows,:)=[tmp*B,C(rows-n, 1:end-p)];  %将用上一部分元素将C矩阵填满
    tmp = A*tmp;  %每一次将tmp左乘一次A
    M(rows,:) = tmp;  %将M矩阵写满
end
%定义Q_bar和
S_q=size(Q,1); %找到Q的维度
S_r=size(R,1);  %找到R的维度
Q_bar= zeros((N+1)*S_q,(N+1)*S_q);   %初始化Q_bar为全0矩阵
for i = 0:N
    Q_bar(i*S_q+1:(i+1)*S_q,i*S_q+1:(i+1)*S_q)=Q;  %将Q写到Q_bar的对角线上
end

%定义R_bar
R_bar = zeros(N*S_r,N*S_r);
for i = 0:N-1
    R_bar(i*S_r+1:(i+1)*S_r,i*S_r+1:(i+1)*S_r)=R;
end

G=M'*Q_bar*M;
E=C'*Q_bar*M;
H=C'*Q_bar*C+R_bar;

f=(x_k'*E')';
U_k=quadprog(H,f);

end
    