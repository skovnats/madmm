function [  evecs, evals ] = giveEigVV( W, A, max_num_evecs )
% Receives Laplacian matrixes ( i.e. L = inv(A)W ), and computes general
% eigen problem Wf = lAf

num_vert = size(  W, 1 );
num_evecs = min(num_vert - 1, max_num_evecs);

tic
warning off
if isscalar(A)
    [evecs, evals] = eigs(W/A, num_evecs, 'SM', struct('disp', 0));
%      [evecs, evals] = eigs(W, spdiags(A*ones( num_vert, 1 ), 0, num_vert, num_vert), num_evecs, 'SM', struct('disp', 0));
else
    try
    [evecs, evals] = eigs( W + A, A, num_evecs, 'SM', struct('disp', 0));
    catch
    [evecs, evals] = eigs( W + A, A, num_evecs, -0.5);   
    end
end
warning on

tc=toc; fprintf('done eigs(%6.2f sec)\n', tc);
if isscalar(A)
    evals = diag(evals);
else
    evals = diag(evals) - 1;
end

if ~all(evals > -5e-10)
    %%%
    % UPD: 22.2011 - replaced option of fixing the problem with negative values. 
    %%%
    %evals = max( abs(evals),  -5e-5 );
    mn = min( evals );
    evals = min( evals+abs(mn), 1 );
end
%assert(all(evals >= -5e-5));

% [evecs_, evals_] = eigs(W, A, num_evecs, 0, struct('disp', 0));
% evals_ = diag(evals_);
% max(abs(evals - evals_)) 

[evals, idx] = sort(evals);
evecs = evecs(:,idx);

fprintf('min eig is %d, max is %d \n \n', min(evals), max( evals ));
