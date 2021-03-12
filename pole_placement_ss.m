function [K,sys_co,Kco] = pole_placement_ss(sys,desired_poles)
A = sys.a;B = sys.b;C = sys.c;D = sys.d;
Co = ctrb(sys.a,sys.b);
n = length(A(1,:));
if length(B(1,:))~=1
    display('For now only SISO systems are supported.');
    display('Returning the input system. Exiting...');
    K = zeros(1,n);
    sys_co = ss(A,B,C,D);
    Kco = K;
else
%rank(Co);
syms s;
if n==rank(Co)
    Inv_Co = inv(Co);
    q = Inv_Co(end,:);
    T = zeros(n,n);
    for i=1:n
        T(i,:) = q*A^(i-1);
    end 

Aco = T*A/T;
a_coeff = -Aco(end,:);

a_d_coeff = double(coeffs(prod(s-desired_poles)));
Kco = a_coeff - a_d_coeff(1:n);
K = Kco*T;
Aco = T*A/T;
Bco = T*B;
Cco = C/T;
sys_co = ss(Aco,Bco,Cco,D);
else
    [R,basiccol] = rref(Co);
    r = rank(Co);
    Q = Co(:,basiccol);
    %extras = n - r;
    I = eye(n);
    r_temp = r;
    for i=1:n
        temp = [Q I(:,i)];
        if rank(temp)>r_temp
            Q = temp;
            r_temp = r_temp + 1;
            if rank(temp)==n
                break;
            end
        end
    end
    T = inv(Q);
    Aco = Q\A/T;
    Bco = Q\B;
    Cco = C/T;
    Ac = Aco(1:r,1:r);
    Anotc = Aco(r+1:n,r+1:n);
    Bc = Bco(1:r,:);
    Bnotc = Bco(r+1:n,:);
    Cc = Cco(:,1:r);
    Cnotc = Cco(:,r+1:n);
    %transfer = Cc*inv(s*eye(r)-Ac)*Bc + D;
    %feval(symengine, 'normal', transfer)
    sys_co = ss(Aco,Bco,Cco,D);
    if not(any(eig(Anotc)> 0))
        display('The system is not controlable but it is stabilisable.');
        display('Only the controlable eigenvalues are pushed to the desired values');
        K_syms = sym('K',[1 r]);
        a_c = coeffs(det(s*eye(r)-(Ac+Bc*K_syms)),s);
        a_d_c = coeffs(prod(s-desired_poles(1:r)));
        temp = struct2cell(solve(a_c == a_d_c));
        Kco = [cat(r,temp{:}) zeros(1,n-r)];
        K = double(Kco)/Q;
    else
        display('The system is not stabilisable.');
        display('One or more eigenvalues of the non-controlable part have positive value.');
        display('Returning the same system')
        K = zeros(1,n);
        sys_co = ss(A,B,C,D);
        Kco = K;
    end
    
end
end
end