I = [1 0;0 1]
ket0 = [1 0;0 0]
ket1 = [0 0;0 1]
ketp = [1 1;1 1]
ketm = [1 -1;-1 1]

outer0 = kron(kron(ket0,ket0),ket0)
outer1 = kron(kron(ket1,ket1),ket1)
outerp = (1/8)*(kron(kron(ketp,ketp),ketp))
outerm = (1/8)*(kron(kron(ketm,ketm),ketm))

Q = (1/4)*((outer0 * outer0') + (outer1 * outer1') + (outerp * outerp') + (outerm * outerm'))
n = 8

%for solving just the primal problem
% cvx_begin sdp
%     variable X(n,n) hermitian semidefinite
%     maximize (trace(Q*X))
%     X >= 0
%     m = TrX(X,1,[4,2]) 
%     m == I
%     
% cvx_end

cvx_begin sdp
    variable X(n,n) hermitian semidefinite
    dual variable y
    maximize (trace(Q*X))
    y: TrX(X,1,[4,2]) == I
    
cvx_end

fprintf('Primal variable(X) :\n')
disp(X)
fprintf('Dual variable(Y) :\n')
disp(y)


function x = TrX(p,sys,dim)

% check arguments
if any(sys > length(dim)) || any(sys < 0)
  error('Invalid subsystem in SYS')
end
disp(prod(dim));
if (length(dim) == 1 && mod(length(p)/dim,1) ~= 0)...
  || length(p) ~= prod(dim)
  error('Size of state PSI inconsistent with DIM');
end


% remove singleton dimensions
if exist('setdiff')
  % matlab
  sys = setdiff(sys,find(dim == 1));
else
  % octave
  sys = complement(find(dim == 1),sys);
end
dim = dim(find(dim ~= 1));


% calculate systems, dimensions, etc.
n = length(dim);
rdim = dim(end:-1:1);
keep = [1:n];
keep(sys) = [];
dimtrace = prod(dim(sys));
dimkeep = length(p)/dimtrace;


if any(size(p) == 1)
  % state vector
  if size(p,1) == 1
    p = p';
  end
  % reshape state vector to "reverse" ket on traced subsystems into a bra,
  % then take outer product
  perm = n+1-[keep(end:-1:1),sys];
  x = reshape(permute(reshape(p,rdim),perm),[dimkeep,dimtrace]);
  x = x*x';


else

  perm = n+1-[keep(end:-1:1),keep(end:-1:1)-n,sys,sys-n];
  x = reshape(permute(reshape(p,[rdim,rdim]),perm),...
              [dimkeep,dimkeep,dimtrace^2]);
  x = sum(x(:,:,[1:dimtrace+1:dimtrace^2]),3);

    end

end
    
    
    
