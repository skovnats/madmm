function [ C ] = calc_CMaks( FX0, FY0, basis1, basis2, AX, AY, WX, WY )
% Implementation of the Maks's Ovsjannikov's et al. approach for
% calculating the functional map C. 
%% NOTE: I assume that the matrix C is orthonormal
%% |AA C - BB|
FX=[];FY=[];
%% if FX and FY are cells, make them vectors/matrices
if iscell(FX0)
   for i=1:numel(FX0)
      FX=[FX,FX0{i}]; 
      FY=[FY,FY0{i}]; 
   end
else
    FX=FX0;
    FY=FY0;
    clear FX0 FY0;
end

%% Calculate the Fourier coefficients
BX=basis1.' * AX * FX;
BY=basis2.' * AY * FY;

%% Place the variables in the equation
AA=BX.';
BB=BY.';

%% Terms for commutativity
if nargin > 6
    D1 = basis1.' * AX * WX * basis1;
    D2 = basis2.' * AY * WY * basis2;
    %% Combining everything 
    AA=[AA;(BX.')*D1];
    BB=[BB;(BY.')*D2];
end

%% Solve C
C = solveICP(AA, BB);
end