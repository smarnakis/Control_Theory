A = [1 3;3 1];
B = [1 0]';
C = [0 1];
D = [0];
%syms s;
%temp = C*inv(eye(2)*s-A)*B
%[ N , D ] = numden(temp);
%sys = tf(sym2poly(N),sym2poly(D))
[n,d] = ss2tf(A,B,C,0);
sys = tf(n,d)
sys = ss(A,B,C,D);
Co = ctrb(sys.a,sys.b)
rank(Co)
temp_sys = canon(sys,'companion')
csys = ss(temp_sys.a',temp_sys.c',temp_sys.b',temp_sys.d)

%%
A = [1 3;3 1];
B = [1 0]';
C = [0 1];
D = [0];
sys = ss(A,B,C,D);
Co = ctrb(sys.a,sys.b)
n = length(A(1,:));
rank(Co);
if n==rank(Co)
    Inv_Co = inv(Co);
    q = Inv_Co(end,:);
    T = zeros(n,n);
    for i=1:n
        T(i,:) = q*A^(i-1);
    end 
    Bco = eye(n,n);
    Bco = Bco(:,end);
end
Aco = T*A*inv(T)
a_coed = Aco(end,:)
lamda_desired = [-1 -2];
syms s;
a_d = double(coeffs(prod(s-lamda_desired)))
Kco = -a_coed - a_d(1:n)
%%
clc
A = [1 0 2;0 -1 1;1 1 -2];
B = [1 0 0]';
C = [1 1 1];
D = [0];
sys = ss(A,B,C,D);
desired_poles = [-4 -4 -3];
[K,sys_co,Kco] = pole_placement_ss(sys,desired_poles);
eig(A+B*K)
%% Ackermann (only for rank(Co)==n)
en = eye(size(A));
en = en(end,:)
Co = ctrb(sys.a,sys.b);
Inv_Co = inv(Co);
syms s;
a_d_coeff = (double(coeffs(prod(s-desired_poles))))
x_d_A = zeros(size(A));
for i=1:length(a_d_coeff)
    x_d_A = x_d_A + a_d_coeff(i)*A^(i-1)
end
K = -en*Inv_Co*x_d_A
eig(A+B*K)