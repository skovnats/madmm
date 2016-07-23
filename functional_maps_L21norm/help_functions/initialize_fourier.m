function [A, B] = initialize_fourier( Z, P, Phi, Ax, Q, Psi, Ay)
% Function get's matrices P, Phi, Ax, Q, Psi, Ay and find diagonal matrices A, B (K x Z, Z <= K)
% s.t. A,B = argmin | P'*Ax*Phi*A - Q'*Ay*Psi*B |

%
K = size( Phi, 2 );

A = zeros(K, Z);
B = A;

mtr_ptions = [1 1;
    1 -1;
    -1 1;
    -1 -1];
for i = 1:Z
    PPhi_i = P'*Ax*Phi( :, i );
    QPsi_i = Q'*Ay*Psi( :, i );
    
    for j = 1:4
        opt(j) = norm(mtr_ptions(j,1)*PPhi_i - mtr_ptions(j,2)*QPsi_i, 'fro');
    end
    [mn, ind] = min(opt); ind = ind(1);
    
    A( i, i ) = mtr_ptions( ind, 1 );
    B( i, i ) = mtr_ptions( ind, 2 );
end