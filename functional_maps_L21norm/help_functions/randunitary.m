function [U] = randunitary(k,n)
% Generates random unitary matrix kxk

a=rand(k);
a=a+a';
[U,~] = eig(a);

if exist('n','var')
   U=[U;zeros(n-k,k)]; 
end
end