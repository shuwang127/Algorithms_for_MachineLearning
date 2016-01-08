function  [v,d,flag] = eigs(varargin)
%EIGS   Find a few eigenvalues and eigenvectors.
%   EIGS solves the eigenvalue problem A*v = lambda*v or the generalized
%   eigenvalue problem A*v = lambda*B*v.  Only a few selected eigenvalues,
%   or eigenvalues and eigenvectors, are computed.
%
%   [V,D,FLAG] = EIGS(A)
%   [V,D,FLAG] = EIGS('Afun',N)
%
%   The first input argument is either a square matrix (which can be
%   full or sparse, symmetric or nonsymmetric, real or complex), or a
%   string containing the name of an M-file which applies a linear
%   operator to the columns of a given matrix.  In the latter case,
%   the second input argument must be N, the order of the problem.
%   For example, EIGS('fft',...) is much faster than EIGS(F,...)
%   where F is the explicit FFT matrix.
%
%   The remaining input arguments are optional and can be given in
%   practically any order:
%
%   [V,D,FLAG] = EIGS(A,B,K,SIGMA,OPTIONS)
%   [V,D,FLAG] = EIGS('Afun',N,B,K,SIGMA,OPTIONS)
%
%   where
%
%       B         A symmetric positive definite matrix the same size as A.
%       K         An integer, the number of eigenvalues desired.
%       SIGMA     A scalar shift or a two letter string.
%       OPTIONS   A structure containing additional parameters.
%
%   With one output argument, D is a vector containing K eigenvalues.
%   With two output arguments, D is a K-by-K diagonal matrix and V is
%   a matrix with K columns so that A*V = V*D or A*V = B*V*D.
%   With three output arguments, FLAG indicates whether or not the
%   eigenvalues converged to the desired tolerance.  FLAG = 0 indicated
%   convergence, FLAG = 1 not.  FLAG = 2 indicates that EIGS stagnated
%   i.e. two consecutive iterates were the same.
%
%   If B is not specified, B = eye(size(A)) is used.  Note that if B is
%   specified, then its Cholesky factorization is computed.
%
%   If K is not specified, K = MIN(N,6) eigenvalues are computed.
%
%   If SIGMA is not specified, the K-th eigenvalues largest in magnitude
%   are computed.  If SIGMA is zero, the K-th eigenvalues smallest in
%   magnitude are computed.  If SIGMA is a real or complex scalar, the
%   "shift", the K-th eigenvalues nearest SIGMA are computed using
%   shift-invert.  Note that this requires computation of the LU
%   factorization of A-SIGMA*B.  If SIGMA is one of the following strings,
%   it specifies the desired eigenvalues.
%
%       SIGMA            Specified eigenvalues
%
%       'LM'             Largest Magnitude  (the default)
%       'SM'             Smallest Magnitude (same as sigma = 0)
%       'LR'             Largest Real part
%       'SR'             Smallest Real part
%       'BE'             Both Ends.  Computes k/2 eigenvalues
%                        from each end of the spectrum (one more
%                        from the high end if k is odd.)
%
%   The OPTIONS structure specifies certain parameters in the algorithm.
%
%     Field name       Parameter                             Default
%
%     OPTIONS.tol      Convergence tolerance:                1e-10 (symmetric)
%                      norm(A*V-V*D,1) <= tol * norm(A,1)    1e-6 (nonsymm.)
%     OPTIONS.stagtol  Stagnation tolerance: quit when       1e-6
%                      consecutive iterates are the same.
%     OPTIONS.p        Dimension of the Arnoldi basis.       2*k
%     OPTIONS.maxit    Maximum number of iterations.         300
%     OPTIONS.disp     Number of eigenvalues displayed       20
%                      at each iteration.  Set to 0 for
%                      no intermediate output.
%     OPTIONS.issym    Positive if Afun is symmetric.        0
%     OPTIONS.cheb     Positive if A is a string,            0
%                      sigma is 'LR','SR', or a shift, and
%                      polynomial acceleration should be
%                      applied.
%     OPTIONS.v0       Starting vector for the Arnoldi       rand(n,1)-.5
%                      factorization
%
%   See also EIG, SVDS.

%   Richard J. Radke and Dan Sorensen.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.19 $  $Date: 1998/08/14 18:50:27 $

global ChebyStruct   % Variables used by accpoly.m 

%  ======   Check input values and set defaults where needed

A = varargin{1};
if isstr(A)
   n = varargin{2};
else
   [n,n] = size(A);
   if any(size(A) ~= n)
      error('Matrix A must be a square matrix or a string.')
   end
end

if n == 0
   v = [];
   d = [];
   flag = 0;
   return
end

% No need to iterate - also elimiates confusion with
% remaining input arguments (they are all n-by-n).

if n == 1
   B = 1;
   for j = (2+isstr(A)):nargin
      if isequal(size(varargin{j}),[n,n]) & ~isstruct(varargin{j})
         B = varargin{j};
         if ~isreal(B)
            error('B must be real and symmetric positive definite.')
      	end
         break
      end
   end
   if isstr(A)
      if nargout < 2
         v = feval(A,1) / B;
      else
         v = 1;
         d = feval(A,1) / B;
         flag = 0;
      end
      return
   else
      if nargout < 2
         v = A / B;
      else
         v = 1;
         d = A / B;
         flag = 0;
      end
      return
   end
end

B = [];
b2 = 1;
k = [];
sigma = [];
options = [];
for j = (2+isstr(A)):nargin
   if isequal(size(varargin{j}),[n,n])
      B = varargin{j};
      if ~isreal(B)
         error('B must be real and symmetric positive definite.')
      end
   elseif isstruct(varargin{j})
      options = varargin{j};
   elseif isstr(varargin{j}) & length(varargin{j}) == 2
      sigma = varargin{j};
   elseif max(size(varargin{j})) == 1
      s = varargin{j};
      if isempty(k) & isreal(s) & (s == fix(s)) & (s > 0)
         k = min(n,s);
      else
         sigma = s;
      end
   else
      error(['Input argument number ' int2str(j) ' not recognized.'])
   end
end

% Set defaults

if isempty(k)
   k = min(n,6);
end
if isempty(sigma)
   sigma = 'LM';
end
if isempty(options)
   fopts = [];
else
   fopts = fieldnames(options);
end

% A is the matrix of all zeros (not detectable if A is a string)

if ~isstr(A)
   if nnz(A) == 0
      if nargout < 2
         v = zeros(k,1);
      else
         v = eye(n,k);
         d = zeros(k,k);
         flag = 0;
      end
      return
   end
end

%  Trick to get faster convergence.  The tail end (closest to cut off
%  of the sort) will not converge as fast as the leading end of the
%  "wanted" Ritz values.

ksave = k;
k = min(n-1,k+3);

% Set issym, specifying a symmetric matrix or operator.

if ~isstr(A)
   issym = isequal(A,A');
elseif strmatch('issym',fopts)
   issym = options.issym;
else
   issym = 0;
end

% Set tol, the convergence tolerance.

if strmatch('tol',fopts);
   tol = options.tol;
else
   if issym
      tol = 1e-10;         % Default tol for symmetric problems.
   else
      tol = 1e-6;          % Default tol for nonsymmetric problems.
   end      
end

% Set stagtol, the stagnation tolerance.

if strmatch('stagtol',fopts);
   stagtol = options.stagtol;
else
   stagtol = 1e-6;
end

% Set p, the dimension of the Arnoldi basis.

if strmatch('p',fopts);
   p = options.p;
else
   p = min(max(2*k,20),n);
end

% Set maxit, the maximum number of iterations.

if strmatch('maxit',fopts);
   maxit = options.maxit;
else
   maxit = max(300,2*n/p);
end

% Set display option.

if strmatch('disp',fopts);
   dispn = options.disp;
else
   dispn = 20;
end

if strmatch('v0',fopts);
   v0 = options.v0;
else
   v0 = rand(n,1)-.5;
end

% Check if Chebyshev iteration is requested.

if strmatch('cheb',fopts);
   cheb = options.cheb;
else
   cheb = 0;
end

if cheb & issym
   if isstr(A) & strcmp(sigma,'LR')
      sigma = 'LO';
      ChebyStruct.sig = 'LR';
      ChebyStruct.filtpoly = [];
   elseif isstr(A) & strcmp(sigma,'SR')
      sigma = 'SO';
      ChebyStruct.sig = 'SR';
      ChebyStruct.filtpoly = [];
   elseif isstr(A) & ~isstr(sigma)
      ChebyStruct.sig = sigma;
      sigma = 'LO';
      ChebyStruct.filtpoly = [];
   elseif isstr(A)
      warning('Cheb option only available for symmetric operators and when sigma is ''LR'', ''SR'', or numeric.  Proceeding without cheb.')
      cheb = 0;
   end
elseif cheb
   warning('Cheb option only available for symmetric operators and when sigma is ''LR'', ''SR'', or numeric.  Proceeding without cheb.')
   cheb = 0;
end

if isstr(sigma)
   sigma = upper(sigma);
end

if strcmp(sigma,'SM') & ~isstr(A)
   sigma = 0;
end

if k > n | k < 0
   error('k must be an integer between 0 and n.')
elseif isstr(A) & ~isstr(sigma) & ~cheb
   error('Operator A and numeric sigma not compatible (unless options.cheb > 0).')
elseif fix(k) ~= k
   error('k must be an integer between 0 and n.')
elseif p > n
   error('p must satisfy p <= n.')
elseif issym & ~isreal(sigma)
   error('Use real sigma for a symmetric matrix.')
end

if ~isstr(A) & ~isempty(B);
   pmmd = symmmd(spones(A) | spones(B));
   A = A(pmmd,pmmd);
   B = B(pmmd,pmmd);
elseif ~isstr(A)
   pmmd = symmmd(A);
   A = A(pmmd,pmmd);
end   

B2 = B;

if ~isequal(B,B')
   error('B must be real and symmetric positive definite.')
end
if ~isempty(B)
   [BL,b2] = chol(B);
   if b2
      error('B must be real and symmetric positive definite.')
   end
end

info = 0;
v = v0;
p = p - k;
psave = p;
op = A;
ChebyStruct.op = op;

if ~isempty(B)
   beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
else
   beta = 1.0/norm(v(:,1));
end

v(:,1) = v(:,1)*beta;

%  Compute w = Av;

if ~isempty(B)
   if ~isstr(A) & ~isstr(sigma),
      [L,U] = lu(A - sigma*B);
      condU = condest(U);
      dsigma = n * full(max(max(abs(A)))) * eps;
      if sigma < 0
         sgnsig = -1;
      else
         sgnsig = 1;
      end
      sigitr = 1;
      while condU > 1/eps & ((dsigma <= 1 & sigitr <= 10) | ~isfinite(condU))
         disps1 = sprintf(['sigma = %10e is near an exact solution ' ...
               'to the generalized eigenvalue problem of A and B,\n' ...
               'so we cannot use the LU factorization of (A-sigma*B): ' ...
               'condest(U) = %10e.\n'], ...
            sigma,condU);
         if abs(sigma) < 1
            sigma = sigma + sgnsig * dsigma;
            disps2 = sprintf('We are trying sigma + %10e = %10e instead.\n', ...
               sgnsig*dsigma,sigma);
         else
            sigma = sigma * (1 + dsigma);
            disps2 = sprintf('We are trying sigma * (1 + %10e) = %10e instead.\n', ...
               dsigma,sigma);
         end
         if nargout < 3 & dispn ~= 0             
            disp([disps1 disps2])
         end   
         [L,U] = lu(A - sigma*B);
         condU = condest(U);
         dsigma = 10 * dsigma;
         sigitr = sigitr + 1;
      end
      if sigitr > 1 & nargout < 3 & dispn ~= 0
         disps = sprintf(['LU factorization of (A-sigma*I) for ' ...
               'sigma = %10e: condest(U) = %10e\n'],sigma,condU);
         disp(disps)             
      end
      v(:,1) = U \ (L \ (B*v(:,1)));
      if sum(~isfinite(v(:,1)))
         d = Inf;
         if nargout <= 1
            v = d;
         end
         flag = 1;
         return
      end
      beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
      v(:,1) = v(:,1)*beta;
      v(:,2) = U \ (L \ (B*v(:,1)));
      if sum(~isfinite(v(:,2)))
         d = Inf;
         if nargout <= 1
            v = d;
         end
         flag = 1;
         return
      end
   elseif isstr(A) & b2
      v(:,2) = B \ feval(A,v(:,1));
   elseif isstr(A) & ~b2
      v(:,2) = BL \ feval(A,(BL'\v(:,1)));
   elseif ~b2
      v(:,2) = BL \ (A*(BL'\v(:,1)));
   else
      v(:,2) = B \ (A*v(:,1));
   end
else
   if ~isstr(A) & ~isstr(sigma),
      [L,U] = lu(A - sigma*speye(n));
      condU = condest(U);
      dsigma = n * full(max(max(abs(A)))) * eps;
      if sigma < 0
         sgnsig = -1;
      else
         sgnsig = 1;
      end
      sigitr = 1;
      while condU > 1/eps & ((dsigma <= 1 & sigitr <= 10) | ~isfinite(condU))
         disps1 = sprintf(['sigma = %10e is near an exact eigenvalue of A,\n' ...
               'so we cannot use the LU factorization of (A-sigma*I): ' ...
               ' condest(U) = %10e.\n'],sigma,condU);
         if abs(sigma) < 1
            sigma = sigma + sgnsig * dsigma;
            disps2 = sprintf('We are trying sigma + %10e = %10e instead.\n', ...
               sgnsig*dsigma,sigma);
         else
            sigma = sigma * (1 + dsigma);
            disps2 = sprintf('We are trying sigma * (1 + %10e) = %10e instead.\n', ...
               dsigma,sigma);
         end
         if nargout < 3 & dispn ~= 0
            disp([disps1 disps2])
         end
         [L,U] = lu(A - sigma*speye(n));
         condU = condest(U);
         dsigma = 10 * dsigma;
         sigitr = sigitr + 1;
      end
      if sigitr > 1 & nargout < 3 & dispn ~= 0
         disps = sprintf(['LU factorization of (A-sigma*I) for ' ...
               'sigma = %10e: condest(U) = %10e\n'],sigma,condU);
         disp(disps)             
      end         
      v(:,1) = U \ (L \ v(:,1));
      if sum(~isfinite(v(:,1)))
         d = Inf;
         if nargout <= 1
            v = d;
         end
         flag = 1;
         return
      end
      beta = 1.0/norm(v(:,1));
      v(:,1) = v(:,1)*beta;
      v(:,2) = U \ (L \ v(:,1));
      if sum(~isfinite(v(:,2)))
         d = Inf;
         if nargout <= 1
            v = d;
         else
            v = v(:,1);
         end
         flag = 1;
         return
      end
   elseif isstr(A)
      v(:,2) = feval(A,v(:,1));
   else
      v(:,2) = A*v(:,1);
   end
end

if ~isempty(B)
   beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
else
   beta = 1.0/norm(v(:,1));
end

v(:,1) = v(:,1)*beta;

if ~isempty(B)
   alpha = v(:,1)'*B*v(:,2);
else
   alpha = v(:,1)'*v(:,2);
end

h(1,1) = alpha;

%  Compute the residual in the second column of V

v(:,2) = v(:,2) - alpha*v(:,1);

%  Perform one step of iterative refinement to correct any
%  orthogonality problems

if ~isempty(B)
   alpha = v(:,1)'*B*v(:,2);
else
   alpha = v(:,1)'*v(:,2);
end

v(:,2) = v(:,2) - alpha*v(:,1);

h(1,1) = h(1,1) + alpha;

%  Compute k steps of the Arnoldi sequence

kstart = 1;
ritz = 1.0;
kp1 = k + 1;
kend = k + p;
k1 = 1;

if ~isempty(B)
   if ~isstr(A) & ~isstr(sigma)
      [v,h,info] =  arnold2(k1,k,L,U,B,v,h,tol);
      if info == -1
         d = h(1,1);
         if nargout <= 1
            v = d;
         else
            v = v(:,1);
         end
         flag = 1;
         return
      end
   elseif ~b2
      [v,h] =  arnold(k1,k,A,BL,v,h,1);
   else
      [v,h] =  arnold(k1,k,A,B,v,h,0);
   end
else
   if ~isstr(A) & ~isstr(sigma)
      [v,h,info] =  arnold2nob(k1,k,L,U,v,h,tol);
      if info == -1
         d = h(1,1);
         if nargout <= 1
            v = d;
         else
            v = v(:,1);
         end
         flag = 1;
         return
      end
   else
      [v,h] =  arnoldnob(k1,k,A,v,h);
   end
end

%  Now update the Arnoldi sequence in place

iter = 0;
knew = k;
ritzes = zeros(ksave,1);
ritzests = ones(ksave,1);
stopcrit = 1;
beta = 1;
betanew = 1;
residest = 1;
stagnation = 0;

%  MAIN LOOP

while (((stopcrit > tol) & (iter < maxit) & ~stagnation) | iter < 2)
   
   iter = iter + 1;
   if dispn
      iter
   end
   
   ChebyStruct.iter = iter;
   
   %     Compute p additional steps of the Arnoldi sequence
   
   kold = k;
   k = knew;
   
   if ~isempty(B)
      if ~isstr(A) & ~isstr(sigma)
         [v,h,info] = arnold2(k,kend,L,U,B,v,h,tol);
         if info == -1
            iter = iter - 1;
            flag = 1;
            break
         end
      elseif ~b2
         [v,h,info] = arnold(k,kend,A,BL,v,h,1);
      else
         [v,h,info] =  arnold(k,kend,A,B,v,h,0);
      end
   else
      if ~isstr(A) & ~isstr(sigma)
         [v,h,info] = arnold2nob(k,kend,L,U,v,h,tol);
         if info == -1
            iter = iter - 1;
            flag = 1;
            break
         end
      else
         [v,h,info] = arnoldnob(k,kend,A,v,h);
  		end
   end
   
   %     If we're doing Chebyshev iteration, when the Ritz estimates
   %     on the extreme values of the spectrum are good enough, then
   %     apply the interpolant polynomial instead.
   
   if cheb == 1 & (ritzests(1) < .1) & (ritzests(2) < .1)
      
      if dispn, disp('Starting polynomial acceleration'), end
      
      A = 'accpoly';
      eigh = eig(h);       
      sig = sigma;
      iter = 1;
      ChebyStruct.iter = iter;
      
      %  Restart the basis with the largest/smallest
      %  ritz vector, as appropriate.
      
      v = v(:,2);
      if ~isempty(B)
         beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
      else
         beta = 1.0/norm(v(:,1));
      end
      
      v(:,1) = v(:,1)*beta;
      
      if ~isempty(B)
         if isstr(A) & b2
            v(:,2) = B \ feval(A,v(:,1));
         elseif isstr(A) & ~b2
            v(:,2) = BL \ feval(A,(BL'\v(:,1)));
         elseif ~b2
            v(:,2) = BL \ (A*(BL'\v(:,1)));
         else
            v(:,2) = B \ (A*v(:,1));
         end
         
         beta = 1.0/sqrt(v(:,1)'*B*v(:,1));
         v(:,1) = v(:,1)*beta;
         alpha = v(:,1)'*B*v(:,2);
         h = alpha;
         v(:,2) = v(:,2) - alpha*v(:,1);
         alpha = v(:,1)'*B*v(:,2);
         v(:,2) = v(:,2) - alpha*v(:,1);
         h = h + alpha;
      else
         if isstr(A)
            v(:,2) = feval(A,v(:,1));
         else
            v(:,2) = A*v(:,1);
         end
         
         beta = 1.0/norm(v(:,1));
         v(:,1) = v(:,1)*beta;
         alpha = v(:,1)'*v(:,2);
         h = alpha;
         v(:,2) = v(:,2) - alpha*v(:,1);
         alpha = v(:,1)'*v(:,2);
         v(:,2) = v(:,2) - alpha*v(:,1);
         h = h + alpha;
      end
      
      if ~b2
         [v,h,info] =  arnold(1,kend,A,BL,v,h,1);
      elseif ~isempty(B)
         [v,h,info] =  arnold(1,kend,A,B,v,h,0);
      else
         [v,h,info] =  arnoldnob(1,kend,A,v,h);
      end
      
      k1 = 1;
      knew = k;
      ritzes = zeros(ksave,1);
      cheb = -1;
      residest = 1;
      sigma = 'LR';
      
   end
   
   k = kold;
   
   %     Compute p shifts based on sigma
   
   %     If A is symmetric, keep H tridiagonal (to avoid spurious
   %     imaginary parts).
   
   if issym
      for i=1:kend-2
         h(i,i+2:kend) = zeros(1,kend-i-1);
         hav = mean([h(i,i+1),h(i+1,i)]);
         h(i,i+1) = hav;
         h(i+1,i) = hav;
      end
      hav = mean([h(kend,kend-1),h(kend-1,kend)]);
      h(kend,kend-1) = hav;
      h(kend-1,kend) = hav;
   end   
   
   [w,q1] = shftit(h,kstart,kend,sigma);
   
   %     Update the command window with current eigenvalue estimates
   
   ritzesold = ritzes;
   
   if ~isstr(sigma)
      ritzes = sigma + ...
         1./w((kend-kstart+1):-1:(kend-ksave+1));
      ritzes = [(sigma+1./eig(h(1:kstart-1,1:kstart-1)));...
            ritzes];
   else
      ritzes = w((kend-kstart+1):-1:(kend-ksave+1));
      ritzes = [eig(h(1:kstart-1,1:kstart-1));ritzes];
   end
   
   if dispn
      eigs = ritzes(1:min(dispn,length(ritzes)))
   end
   
   if iter == 1 & cheb > 0
      ChebyStruct.lbd = min(ritzes);
      ChebyStruct.ubd = max(ritzes);
   elseif cheb > 0
      ChebyStruct.lbd = min(ChebyStruct.lbd, min(ritzes));
      ChebyStruct.ubd = max(ChebyStruct.ubd, max(ritzes));
   end
   
   [m1,m2]=size(q1);
   ritz = norm(q1(m1,p+2:m1));
   
   if ~isempty(B)
      betanew = sqrt(v(:,kend+1)'*B*v(:,kend+1));
   else
      betanew = norm(v(:,kend+1));
   end
   
   ritznew = betanew*q1(m1,1:m1);
   jj = m1;
   kconv = 0;
   while(jj > 0),
      if (abs(ritznew(jj)) <= tol),
         jj = jj - 1;
         kconv = kconv+1;
      else
         jj = -1;
      end
   end
   kkconv = kconv;
   
   %     The while loop counts the number of converged ritz values.
   
   ritzests = [w abs(ritznew)'];
   
   ritzests = ritzests(size(ritzests,1) : -1 : ...
      size(ritzests,1)-ksave+1,2);
   
   stopcritold = stopcrit;
   stopcrit = max(ritzests);
   residest = norm(ritzests);    
   
   %  At first, we use the Ritz estimates to estimate convergence.
   %  However, as the algorithm converges, the Ritz estimates become
   %  poor estimates of the actual error in each Ritz pair.  So when the
   %  Ritz estimates become too small, we actually form the matrix of
   %  errors || AV - VD || where V are the estimates for the eigenvectors
   %  and eigenvalues.  This is expensive computationally but ensures
   %  that the user gets the desired eigenpairs to the requested
   %  tolerance.
   
   if max(ritzests) <= tol*max(norm(h),1)
      
      if ~b2 & isstr(sigma)
         vee = BL \ v(:,1:kend);
      else
         vee = v(:,1:kend);
      end
      
      if ~isempty(B2)
         if isstr(A)
            vee = vee * q1(1:kend,kend:-1:kend-ksave+1);
            for i = 1 : ksave
               nvi = norm(v(:,i));
               if isfinite(nvi) & nvi ~= 0
                  vee(:,i) = vee(:,i) / nvi;
               end
            end
            Avee = feval(A,vee);
            errmat = Avee - B2 * vee * diag(ritzes);
            residest = norm(errmat,1) / norm(Avee,1);
         else
            vee = vee * q1(1:kend, kend:-1:kend-ksave+1);
            for i = 1 : ksave
               nvi = norm(vee(:,i));
               if isfinite(nvi) & nvi ~= 0
                  vee(:,i) = vee(:,i) / nvi;
               end
            end
            errmat = A * vee - B2 * vee * diag(ritzes);
            residest = norm(errmat,1) / norm(A,1);
         end
      else
         if isstr(A)
            vee = vee * q1(1:kend,kend:-1:kend-ksave+1);
            for i = 1 : ksave
               nvi = norm(v(:,i));
               if isfinite(nvi) & nvi ~= 0
                  vee(:,i) = vee(:,i) / nvi;
               end
            end
            Avee = feval(A,vee);
            errmat = Avee - vee * diag(ritzes);
            residest = norm(errmat,1) / norm(Avee,1);
         else
            vee = vee * q1(1:kend, kend:-1:kend-ksave+1);
            for i = 1 : ksave
               nvi = norm(vee(:,i));
               if isfinite(nvi) & nvi ~= 0
                  vee(:,i) = vee(:,i) / nvi;
               end
            end
            errmat = A * vee - vee * diag(ritzes);
            residest = norm(errmat,1) / norm(A,1);
         end
      end
      
      for ii = 1:length(ritzes);
         ritzests(ii) = norm(errmat(:,ii));
      end
      
      stopcrit = residest;
      
   end
   
   if (abs(stopcritold-stopcrit) < stagtol*abs(stopcrit))
      stagind = (ritzesold~=0) & (ritzes~=0);
      stagn = abs(ritzesold - ritzes);
      stagn(stagind) = stagn(stagind) ./ ritzes(stagind);
      stagnation = max(stagn) < stagtol;
	end

   if (((stopcrit > tol) & ~stagnation) | iter < 2)
      
      %     Apply the p implicit shifts if convergence has not yet 
      %     happened.  Otherwise don't apply them and get out of the
      %     loop on next loop test. We need to keep same test here as
      %     in the main loop test to   avoid applying shifts and then
      %     quitting, which would lead to a wrong size factorization
      %     on return.
      
      %        If some ritz values have converged then
      %        adjust k and p to move the "boundary"
      %        of the filter cutoff.
      
      if kconv > 0
         kk = ksave + 3 + kconv;
         p = max(ceil(psave/3),kend-kk);
         k = kend - p;
      end
      
      if any(any(imag(v))) | any(any(imag(h)))
         [v,h,knew] = apshft1(v,h,w,k,p);
      else
         [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
      end
      
      if ~isempty(B)
         betanew = sqrt(v(:,kp1)'*B*v(:,kp1));
      else
         betanew = norm(v(:,kp1));
      end   
      
   end
   
   if dispn
      stopcrit
      disp('==========================')
   end
   
end         %  End of Arnoldi iteration  (MAIN LOOP)

if stopcrit <= tol
   flag = 0;
else
   if nargout < 3 & dispn ~= 0
      if iter >= maxit
         disp('Exiting: Maximum number of iterations exceeded.')
      elseif stagnation
         disp('Exiting: Two consecutive iterates were the same.')
      else
         disp('Exiting: Eigenvalues did not converge.')        
      end
   else
      if stagnation
         flag = 2;
      else
         flag = 1;
      end
   end
end

%  Compute the eigenvalues and eigenvectors of h
kend = min(kend,size(h,1));
[w,q1] = shftit(h,1,kend,sigma);

k = ksave;
p = psave;

%  Transform the converged eigenvalues back to 
%  the original problem and return them in the diagonal
%  k by k matrix d.

%  Set v equal to the wanted eigenvectors

v = v(:,1:kend) * q1(1:kend, kend:-1:kend-k+1);

if ~b2 & isstr(sigma)
   v = BL \ v;
end

for i = 1:k
   nvi = norm(v(:,i));
   if isfinite(nvi) & nvi ~= 0
      v(:,i) = v(:,i) / nvi;
   end   
end

%  In Chebyshev iteration, we recover the eigenvalues by Rayleigh
%  quotients.  Otherwise, the eigenvalues are recovered from the
%  set of shifts w.

if cheb == -1
   d = [];
   if ~isempty(B2)
      for i=1:ksave
         d = [d; ( (v(:,i)' * feval(op,v(:,i)) ) ./ ...
               (v(:,i)' * B2 * v(:,i)))];
      end
   else
      for i=1:ksave
         d = [d; ( (v(:,i)' * feval(op,v(:,i)) ) ./ ...
               norm(v(:,i)))];
      end
   end
   d = diag(d);
elseif ~isstr(sigma)
   t = sigma + 1./w;
   d = diag(t(kend:-1:kend-k+1));
else
   d = diag(w(kend:-1:kend-k+1));
end

if ~isstr(A)
   v(pmmd,:) = v;
end

if nargout <= 1
   v = diag(d);
end

% ====================================================

function [v,h,k] =  apshft1(v,h,w,k,p)

%  APSHFT1  Apply shifts to update an Arnoldi factorization
%
%  APSHFT1 was designed to be called by EIGS when V or H is complex.
%
%  [V,H,k] = APSHFT1(V,H,w,k,p) implicitly applies the
%     p real shifts held in w to update the existing Arnoldi
%     factorization AV - VH = re'  .
%                               k+p
%  The routine results in
%
%     A(VQ) - (VQ)(Q'HQ) = re'  Q
%                            k+p
%
%  where the orthogonal matrix Q is the product of the Givens
%  rotations resulting from p bulge chase sweeps.
%
%  The updated residual is placed in V(:,k+1) and the updated
%  Arnoldi factorization V <- VQ, H <- Q'HQ is returned.

%  Dan Sorensen and Richard J. Radke, 11/95.

k1 = 1;
kend = k+p;
kp1 = k + 1;
q = zeros(1,kend);
q(kend) = 1.0;

dh = diag(h,-1);
ix = find(dh(1:end-1)==0);     % Find the column indices of 0
ix = [0 ; ix ; kend-1];       % subdiagonals in H.
nx = size(ix,1);           

for jj = 1:p               % Loop over shifts
   for ii = 1:nx-1         % Loop over blocks in H
      k1 = ix(ii)+1; k2 = ix(ii+1);
      c = h(k1,k1) - w(jj);
      s = h(k1+1,k1);
      [G,R] = qr([c;s]);
      for i = k1:k2        % Loop over rows in the block
         if i > k1 
            [G,R] = qr(h(i:i+1,i-1));
            h(i:i+1,i-1) = R(:,1);
         end
         
         %           apply rotation from left to rows of H
         
         h(i:i+1,i:kend) = G'* h(i:i+1,i:kend);
         
         
         %           apply rotation from right to columns of H
         
         ip2 = i+2;
         if ip2 > kend
            ip2 = kend;
         end
         
         h(1:ip2,i:i+1) = h(1:ip2,i:i+1)*G;
         
         %           apply rotation from right to columns of V
         
         v(:,i:i+1) = v(:,i:i+1)*G;
         
         %           accumulate e'  Q so residual can be updated
         %                       k+p
         
         q(i:i+1) = q(i:i+1)*G;
      end
   end
end

%  Update the residual and store in the k+1 -st column of v

v(:,kend+1) = v(:,kend+1)*q(k);
v(:,kp1) = v(:,kend+1) + v(:,kp1)*h(kp1,k);

% ==========================================================

function [v,h,k] =  apshft2(v,h,wr,wi,k,p)

%  APSHFT2  Apply shifts to update an Arnoldi factorization
%
%  APSHFT2 was designed to be called by EIGS when V and H are real.
%
%  [V,H,k] = APSHFT2(V,H,wr,wi,k,p) implicitly applies
%     the p complex shifts given by w = wr + i*wi, to update the
%     existing Arnoldi factorization AV - VH = re'  .
%                                                k+p
%
%  The routine results in
%
%     A(VQ) - (VQ)(Q'HQ) = re'  Q
%                            k+p
%
%  where the orthogonal matrix Q is the product of the Givens
%  rotations resulting from p bulge chase sweeps.
%
%  The updated residual is placed in V(:,k+1) and the updated
%  Arnoldi factorization V <- VQ, H <- Q'HQ is returned.

%  Dan Sorensen and Richard J. Radke, 11/95.

k1 = 1;
kend = k+p;
kp1 = k + 1;
q = zeros(kend,1);
q(kend) = 1.0;
mark = 0;
num = 0;

dh = diag(h,-1);
ix1 = find(dh(1:end-1)==0);    % Find the column indices of 0
ix = [0 ; ix1 ; kend-1];         % subdiagonals in H.
jx = [0 ; ix1 ; kend-2];
nx = size(ix,1);           

for jj = 1:p
   
   %     compute and apply a bulge chase sweep initiated by the
   %     implicit shift held in w(jj)
   
   if (abs(wi(jj)) == 0.0) 
      
      %        apply a real shift using 2 by 2 Givens rotations
      
      for ii = 1:nx-1      % Loop over blocks in H
         k1 = ix(ii)+1;
         k2 = ix(ii+1);
         c = h(k1,k1) - wr(jj);
         s = h(k1+1,k1) ;
         t = norm([c s]);
         if t == 0.0
            c = 1.0;
         else
            c = c/t;
            s = s/t;
         end
         for i = k1:k2
            if i > k1
               t = norm(h(i:i+1,i-1));
               if t == 0.0
                  c = 1.0;
                  s = 0.0;
               else
                  c = h(i,i-1)/t;
                  s = h(i+1,i-1)/t;
                  h(i,i-1) = t;
                  h(i+1,i-1) = 0.0;
               end
            end
            
            %              apply rotation from left to rows of H
            
            G = [c s ; -s c];
            h(i:i+1,i:kend) = G* h(i:i+1,i:kend);
            
            %              apply rotation from right to columns of H
            
            ip2 = i+2;
            if ip2 > kend
               ip2 = kend;
            end
            
            h(1:ip2,i:i+1) =  h(1:ip2,i:i+1)*G';
            
            %              apply rotation from right to columns of V
            
            v(:,i:i+1) =  v(:,i:i+1)*G';
            
            %              accumulate e'  Q so residual can be updated
            %                          k+p
            
            q(i:i+1) =  G*q(i:i+1);
         end
      end
      num = num + 1;
   else
      
      %        Apply a double complex shift using 3 by 3 Householder 
      %        transformations
      
      if (jj == p | mark == 1) ,
         mark = 0;       % skip application of conjugate shift
      else
         num = num + 2;    % mark that a complex conjugate
         mark = 1;         % pair has been applied
         
         for ii = 1:nx-1   % Loop over blocks in H
            k1 = jx(ii)+1;
            k2 = k1+1;
            k3 = jx(ii+1);
            c = h(k1,k1)*h(k1,k1) + h(k1,k2)*h(k2,k1) ...
               - 2.0*wr(jj)*h(k1,k1);
            c = c + wr(jj)^2 + wi(jj)^2;
            s = h(k2,k1)*(h(k1,k1) + h(k2,k2) - 2.0*wr(jj));
            g = h(k2+1,k2)*h(k2,k1);
            t = norm([c s g]);
            sig = -sign(c);
            c = c -t*sig;
            for i = k1:k3
               if i > k1
                  t = norm(h(i:i+2,i-1));
                  sig = -sign(h(i,i-1));
                  c = h(i,i-1) - t*sig;
                  s = h(i+1,i-1);
                  g = h(i+2,i-1);
                  h(i,i-1) = t;
                  h(i+1,i-1) = 0.0;
                  h(i+2,i-1) = 0.0;
               end
               t = norm([c s g]);
               if t ~= 0.0
                  c = c/t;
                  s = s/t;
                  g = g/t;
               end
               z = [c s g]';
               
               %                 apply transformation from left to rows of H
               
               t =  sig*2.0*(z'*h(i:i+2,i:kend));
               h(i:i+2,i:kend) = sig*h(i:i+2,i:kend) - z*t;
               
               %                 apply transformation from right to columns of H
               
               ip3 = i+3;
               if ip3 > kend
                  ip3 = kend;
               end
               t =  sig*2.0*h(1:ip3,i:i+2)*z;
               h(1:ip3,i:i+2) = sig*h(1:ip3,i:i+2) - t*z';
               
               %                 apply transformation from right to columns of V
               
               t =  sig*2.0*v(:,i:i+2)*z;
               v(:,i:i+2) = sig*v(:,i:i+2) - t*z';
               
               %                 accumulate e'  Q so residual can be updated
               %                             k+p
               
               t =  sig*2.0*z'*q(i:i+2);
               q(i:i+2) = sig*q(i:i+2) - z*t;
            end      
         end
         
         %           clean up step with Givens rotation
         
         i = kend-1;
         if i > k1
            t = norm([h(i,i-1) h(i+1,i-1)]);
            if t ~= 0.0
               c = h(i,i-1)/t;
               s = h(i+1,i-1)/t;
            else
               c = 1.0;
               s = 0.0;
            end
            h(i,i-1) = t;
            h(i+1,i-1) = 0.0;
         end
         
         %           apply rotation from left to rows of H
         
         G = [c s ; -s c];
         h(i:i+1,i:kend) = G* h(i:i+1,i:kend);
         
         %           apply rotation from right to columns of H
         
         ip2 = i+2;
         if ip2 > kend
            ip2 = kend;
         end
         h(1:ip2,i:i+1) =  h(1:ip2,i:i+1)*G';
         
         %           apply rotation from right to columns of V
         
         v(:,i:i+1) =  v(:,i:i+1)*G';
         
         %           accumulate e'  Q so residual can be updated
         %                       k+p
         
         q(i:i+1) =  G*q(i:i+1);
      end
   end
end

%  update residual and store in the k+1 -st column of v

k = kend - num;
v(:,kend+1) = v(:,kend+1)*q(k);
if k < size(h,1)
   v(:,k+1) = v(:,kend+1) + v(:,k+1)*h(k+1,k);
end

% ===========================================================

function [v,h,defloc] =  arnold(k1,k2,A,B,v,h,bfactor)

%  ARNOLD  Computes or extends an Arnoldi factorization.
%
%  ARNOLD was designed to be called by EIGS in cases where the
%  shift parameter sigma is non-numeric.  This routine uses power
%  methods and not inverse iteration. 
%
%  Let V and H be the Arnoldi factors from a k1-step Arnoldi
%  process such that
%
%     AV - BVH = re' ,           where r is stored in
%            k1            the k1+1-st column of V.
% 
%  [V,H,defloc] = ARNOLD(k1,k2,A,B,V,H,bfactor)
%
%  extends the existing factorization by k2 - k1 steps, resulting in
%
%          AV - V H = r e' .          where r is stored in
%            +   + +   + k2         the k2+1-st column of V.
%
%  On input and output, the columns of V should be B-orthogonal,
%  i.e. V'BV = I.
%
%  If B can be Cholesky factored as B = LL', where L is lower
%  triangular, then the factor L should be passed to ARNOLD in
%  the B position and the parameter bfactor should be set to true.
%
%  See also EIGS.

%  Dan Sorensen and Richard J. Radke, 11/95

if nargin < 7
   bfactor = 0;
end

if bfactor
   BL = B;
   B = [];
end

defloc = 0;

for j = k1+1:k2,
   jm1 = j-1; 
   jp1 = j+1;
   if bfactor
      vtBL = v(:,j)' * BL;
      beta = norm(vtBL);
    else
      beta = sqrt(v(:,j)'*B*v(:,j));
   end
   h(j,jm1) = beta;
   if (beta <= 10*eps*norm(h(1:jm1,jm1)))     % If beta is "small"
      v(:,j) = rand(size(v,1),1);             % we deflate by
      if bfactor
         s = v(:,1:jm1)'*BL*(v(:,j)'*BL)';
      else
         s = v(:,1:jm1)'*B*v(:,j);               % setting the
      end
      v(:,j) = v(:,j) - v(:,1:jm1)*s;         % corresponding
      if bfactor
         s = v(:,1:jm1)'*BL*(v(:,j)'*BL)';
      else
         s = v(:,1:jm1)'*B*v(:,j);               % subdiagonal of H
      end
      v(:,j) = v(:,j) - v(:,1:jm1)*s;         % equal to 0 and
      if bfactor
         vtBL = v(:,j)'*BL;
         beta = norm(vtBL);
      else
         beta = sqrt(v(:,j)'*B*v(:,j));          % starting the
      end
      beta = 1.0/beta;                        % basis for a new 
      h(j,jm1) = 0.0;                         % invariant subspace.
      defloc = j;
   else
      beta = 1.0/beta;
   end
   v(:,j) = v(:,j)*beta;
   
   %     Compute w = Av and store w in j+1 -st col of V
   
   if bfactor
      v2 = BL \ v(:,j);
   else
      v2 = v(:,j);
   end
   
   if isstr(A)
      w = feval(A,v2);
   else
      w = A*v2;
   end
   
   if bfactor
      v(:,jp1) = BL' \ w;
   else
      v(:,jp1) = B \ w;
   end
   
   %     Compute the next (j-th) column of H
   
   if bfactor
      h(1:j,j) = v(:,1:j)'*BL*(v(:,jp1)'*BL)';
   else
      h(1:j,j) = v(:,1:j)'*B*v(:,jp1);
   end
   
   %     Compute the residual in the j+1 -st col of V
   
   v(:,jp1) = v(:,jp1) - v(:,1:j)*h(1:j,j);
   
   %     Perform one step of iterative refinement to correct the orthogonality
   
   if bfactor
      s = v(:,1:j)'*BL*(v(:,jp1)'*BL)';
     else
      s = v(:,1:j)'*B*v(:,jp1);
   end
   v(:,jp1) = v(:,jp1) - v(:,1:j)*s;
   h(1:j,j) = h(1:j,j) + s;
end

% ==============================================================

function [v,h,defloc] =  arnoldnob(k1,k2,A,v,h)

%  ARNOLDNOB  Computes or extends an Arnoldi factorization.
%
%  ARNOLDNOB is exactly the same as ARNOLD when B = I.
%

defloc = 0;

for j = k1+1:k2,
   jm1 = j-1; 
   jp1 = j+1;
   beta = norm(v(:,j)); 
   h(j,jm1) = beta;
   if (beta <= 10*eps*norm(h(1:jm1,jm1)))    % If beta is "small"
      v(:,j) = rand(size(v,1),1);            % we deflate by
      s = v(:,1:jm1)'*v(:,j);                % setting the
      v(:,j) = v(:,j) - v(:,1:jm1)*s;        % corresponding
      s = v(:,1:jm1)'*v(:,j);                % subdiagonal of H
      v(:,j) = v(:,j) - v(:,1:jm1)*s;        % equal to 0 and
      beta = norm(v(:,j));                   % starting the
      beta = 1.0/beta;                       % basis for a new 
      h(j,jm1) = 0.0;                        % invariant subspace.
      defloc = j;
   else
      beta = 1.0/beta;
   end
   v(:,j) = v(:,j)*beta;
   
   %     Compute w = Av and store w in j+1 -st col of V
   
   v2 = v(:,j);
   
   if isstr(A)
      v(:,jp1) = feval(A,v2);
   else
      v(:,jp1) = A*v2;
   end
   
   %     Compute the next (j-th) column of H
   
   h(1:j,j) = v(:,1:j)'*v(:,jp1);
   
   %     Compute the residual in the j+1 -st col of V
   
   v(:,jp1) = v(:,jp1) - v(:,1:j)*h(1:j,j);
   
   %     Perform one step of iterative refinement to correct the orthogonality
   
   s = v(:,1:j)'*v(:,jp1);
   v(:,jp1) = v(:,jp1) - v(:,1:j)*s;
   h(1:j,j) = h(1:j,j) + s;
end

% ===============================================================

function [v,h,defloc] =  arnold2(k1,k2,L,U,B,v,h,tol);

%  ARNOLD2  Computes or extends an Arnoldi factorization.
%
%  ARNOLD2 was designed to be called by EIGS in cases where the
%  shift parameter sigma is numeric.  This routine uses inverse
%  iteration and the LU factors of (A - sigma*B).
%
%  Let V and H be the Arnoldi factors from a k1-step Arnoldi
%  process such that
%
%    inv(L*U) B V - V H = re'  ,     where r is stored in
%                 k1       the k1+1-st column of V.
% 
%  [V,H,defloc] = ARNOLD2(k1,k2,L,U,B,V,H,tol)
%  extends the existing factorization by k2 - k1 steps, resulting in
%
%         inv(L*U) B V - V H = r e' .     where r is stored in
%                     +   + +   + k2      the k2+1-st column of V.
%
%  On input and output, the columns of V should be B-orthogonal,
%  i.e. V'BV = I.
%
%  The parameter tol provides a tolerance for deciding when to
%  deflate the problem.  If the problem can be deflated, a pointer
%  to the appropriate location in H of the start of the active basis
%  is returned in the output argument defloc.
%
%  If the L|U factors of A-sigma*B are ill-conditioned and return
%  unusable U\(L\v), then v and h remain unchanged and defloc is
%  assigned -1.
%
%  Dan Sorensen and Richard J. Radke, 11/95

defloc = 0;
vold = v;
hold = h;

for j = k1+1:k2,
   jm1 = j-1; 
   jp1 = j+1;
   beta = sqrt(v(:,j)'*B*v(:,j));
   h(j,jm1) = beta;
   if (beta <= tol)
      v(j,j) = 1.0;
      s = v(:,1:jm1)'*B*v(:,j);
      v(:,j) = v(:,j) - v(:,1:jm1)*s;
      s = v(:,1:jm1)'*B*v(:,j); 
      v(:,j) = v(:,j) - v(:,1:jm1)*s;
      beta = sqrt(v(:,j)'*B*v(:,j));
      beta = 1.0/beta;
      defloc = j;
   else
      beta = 1.0/beta;
   end
   v(:,j) = v(:,j)*beta;
   
   %     Compute w = Av and store w in j+1 -st col of V
   
   v(:,jp1) = U \ (L \ (B*v(:,j)));
   
   if sum(~isfinite(v(:,jp1)))
      v = vold;
      h = hold;
      defloc = -1;
      return
   end
   
   %     Compute the next (j-th) column of H
   
   h(1:j,j) = v(:,1:j)'*B*v(:,jp1);
   
   %     Compute the residual in the j+1 -st col of V
   
   v(:,jp1) = v(:,jp1) - v(:,1:j)*h(1:j,j);
   
   %     One step of iterative refinement to correct the orthogonality
   
   s = v(:,1:j)'*B*v(:,jp1);
   v(:,jp1) = v(:,jp1) - v(:,1:j)*s;
   h(1:j,j) = h(1:j,j) + s;
end

% ================================================================

function [v,h,defloc] =  arnold2nob(k1,k2,L,U,v,h,tol);

%  ARNOLD2NOB  Computes or extends an Arnoldi factorization.
%
%  ARNOLD2NOB is the same as ARNOLD2 when B = I.
%

defloc = 0;
vold = v;
hold = h;

for j = k1+1:k2,
   jm1 = j-1; 
   jp1 = j+1;
   beta = norm(v(:,j));
   h(j,jm1) = beta;
   if (beta <= tol)
      v(j,j) = 1.0;
      s = v(:,1:jm1)'*v(:,j);
      v(:,j) = v(:,j) - v(:,1:jm1)*s;
      s = v(:,1:jm1)'*v(:,j); 
      v(:,j) = v(:,j) - v(:,1:jm1)*s;
      beta = sqrt(v(:,j)'*v(:,j));
      beta = 1.0/beta;
      defloc = j;
   else
      beta = 1.0/beta;
   end
   v(:,j) = v(:,j)*beta;
   
   %     Compute w = Av and store w in j+1 -st col of V
   
   v(:,jp1) = U \ (L \ v(:,j));
   
   if sum(~isfinite(v(:,jp1)))
      v = vold;
      h = hold;
      defloc = -1;
      return
   end
   
   %     Compute the next (j-th) column of H
   
   h(1:j,j) = v(:,1:j)'*v(:,jp1);
   
   %     Compute the residual in the j+1 -st col of V
   
   v(:,jp1) = v(:,jp1) - v(:,1:j)*h(1:j,j);
   
   %     One step of iterative refinement to correct the orthogonality
   
   s = v(:,1:j)'*v(:,jp1);
   v(:,jp1) = v(:,jp1) - v(:,1:j)*s;
   h(1:j,j) = h(1:j,j) + s;
end

% ==============================================================

function w = accpoly(v)

%   ACCPOLY applies an accelerant polynomial in the operator
%   specified by the global variable Chebystruct.op 
%       to the input vector v.
%
%   ACCPOLY was designed to be called by EIGS when the input argument
%       A to EIGS is an m-file specifying a symmetric matrix-vector
%       product.
%
%   Approximate lower and upper bounds on the eigenvalues of op
%   should be present in the global variables lbd and ubd.
%   The global variable sig is 'SO' or 'LO' depending on whether
%   the smallest or largest eigenvalues of op are desired.
%   sig can also be a numeric shift.
%
%   The accelerant polynomial in ChebyStruct.op is either a scaled 
%       Chebyshev polynomial over [lbd,ubd] or an equiripple single-peak
%       polynomial obtained by the filter design code contained in bp_fer.

%   Richard J. Radke, 3/96.

global ChebyStruct

op = ChebyStruct.op;
lbd = ChebyStruct.lbd;
ubd = ChebyStruct.ubd;
sig = ChebyStruct.sig;
filtpoly = ChebyStruct.filtpoly;
iter = ChebyStruct.iter;

p = .25;  %  1-p is the fraction of the spectrum that gets mapped
%  into the equiripple region (sig = string)
m = 10;

if ~isstr(sig)
   slope = 2/(ubd-lbd);
   intercept = 1 - (ubd*slope);
   if sig > ubd
      sig = ubd-.01*(ubd-lbd);
   elseif sig < lbd
      sig = lbd+.01*(ubd-lbd);
   end
elseif strcmp(sig,'SO') | strcmp(sig,'SR')
   slope = 2/(ubd-lbd-p*(ubd-lbd));slope = 2/(ubd-p*(ubd-lbd)-lbd);
   intercept = -1 - (lbd*slope);
   intercept = 1 - (ubd*slope);
elseif strcmp(sig,'LO') | strcmp(sig,'LR')
   slope = 2/(ubd-p*(ubd-lbd)-lbd);
   intercept = -1 - (lbd*slope);
end

if isstr(sig)
   
   w0 = v;
   w1 = slope*feval(op,v)+intercept*v;
   
   for jj = 2:m
      w = 2*(slope*feval(op,w1)+intercept*w1)-w0;
      w0 = w1;
      w1 = w;
   end
   
else
   
   if iter == 1
      
      wp = acos(slope*sig+intercept);
      m = 10;
      N = 2*m+1;
      L = 4;
      del = 0.01;
      h = bp_fer(N,L,wp,del,-del,2^7);
      ChebyStruct.filtpoly = h2x(h);
      filtpoly = ChebyStruct.filtpoly;
      iter = 2;
      
   end
   
   w = filtpoly(1)*v;
   
   for i=2:(m+1)
      w = (slope*feval(op,w)+intercept*w) + filtpoly(i)*v;
   end
   
end

% ========================================================

function [w,qq] = shftit(h, kstart, kend, sigma)

%  SHFTIT  Calculate shifts to update an Arnoldi factorization.
%
%  SHFTIT was designed to be called by EIGS.
%
%  Syntax: [w,q] = SHFTIT(H, kstart, kend, sigma) where
%
%  H is an upper Hessenberg matrix (from the Arnoldi factorization
%     AV = VH + fe_k'),
%
%  kstart points to the start of the active block in H,
%
%  kend points to the end of the active block in H,
%
%  sigma is a numeric shift or one of 'LR','SR','LM','SM','BE'.
%     (If sigma is a number, sigma = 'LM' is used, since
%     the problem has already been shifted by sigma.  The
%     real eigenvalues are recovered later.)
%
%  SHFTIT calculates [q,w] = eig(The active block of H), where
%  the eigenvalues and eigenvectors are reordered according to sigma,
%  with the eigenvalues to use as shifts put first.

%  Dan Sorensen and Richard J. Radke, 7/95 

[q,ww] = eig(h(kstart:kend,kstart:kend));

w = diag(ww);

k = kend-kstart+1;

%  select filter mechanism by activating appropriate choice below

if ~isstr(sigma)
   
   [s,ir] = sort(abs(w));
   
elseif strcmp(sigma,'BE')
   
   %     sort on real part, alternating from high to low end of
   %     the spectrum
   
   par = rem(k,2);
   
   [s,ir] = sort(-real(w));
   
   ix(1:2:k-1+par) = 1:ceil(k/2);
   ix(2:2:k-par) = k:-1:(ceil(k/2)+1);
   
   ir = flipud(ir(ix));
   
elseif strcmp(sigma,'LO')
   
   %     keep the k-1 largest and 1 smallest e-vals
   %     at the end of the sort
   
   [s,ir] = sort(real(w));
   ir = ir([2:k,1]);
   
elseif strcmp(sigma,'SO')
   
   %     keep the k-1 smallest and 1 largest e-vals
   %     at the end of the sort
   
   [s,ir] = sort(-real(w));
   ir = ir([2:k,1]);
   
elseif strcmp(sigma,'LR')
   
   %     sort for largest real part (shifts are smallest real part)
   
   [s,ir] = sort(real(w));
   
elseif strcmp(sigma,'SR')
   
   %     sort for smallest real part (shifts are largest real part)
   
   [s,ir] = sort(-real(w));
   
elseif strcmp(sigma,'SM')
   
   %     sort for smallest absolute val (shifts are largest abs val)
   
   [s,ir] = sort(-abs(w));
   
else
   
   %     sort for largest absolute val (shifts are smallest abs val)
   
   [s,ir] = sort(abs(w));
   
end

%  apply sort to w

w = w(ir);
qq = q(:,ir);

% =====================================================

function p = add_poly(p1,p2)
%
%   function p = add_poly(p1,p2);
%   Adds 2 polynomials.
%
l1 = length(p1); l2 = length(p2);
if l1 == 0, p = rlz(p2); break; end
if l2 == 0, p = rlz(p1); break; end
p1 = rlz(p1(:)); p2 = rlz(p2(:));
l1 = length(p1); l2 = length(p2);
if l1 > l2
   p2 = [zeros(l1-l2,1); p2];
else
   p1 = [zeros(l2-l1,1); p1]; 
end
p = rlz((p1 + p2).');

% =======================================================

function [h,h2,rs] = bp_fer(N,L,wp,us,ls,g)
% A program for the design of linear phase bandpass FIR filters with a
% Flat monotonically decreasing Pass band (on either side of wp)
% and  an EquiRipple Stop band.
%       With this program, the user can specify the stop band ripple size
%       but not the passband frequency of flatness.
%
% N  : length of total filter
% L  : degree of flatness
% wp : passband frequency of flatness
% us : upper bound in stop band
% ls : lower boudn in stop band
% g  : grid size
%
% EXAMPLE
%    N  =  25;
%    L  =  8;
%    wp  = .4*pi;
%    us =  0.01;
%    ls = -us;
%    g  =  2^8;
%    [h,h2,rs] = bp_fer(N,L,wp,us,ls,g);
% or
%    N = 31; L = 16; 
%    bp_fer(N,L,wp,us,ls,g);


if (rem(N,2) == 0) | (rem(L,4) ~= 0)
   disp('N must be odd and L must be divisible by 4')
   return
else
   SN = 1e-9;                                      % SN : SMALL NUMBER
   q  = (N-L+1)/2;
   w  = [0:g]'*pi/g;
   ws1 = wp*0.8;
   a = 5; b = 1;
   ws2 = (a*wp+b*pi)/(a+b);
   
   d = ws1/(pi-ws2);            % q1 : number of ref. set freq. in 1st stopband
   q1 = round(q/(1+1/d));       % q2 : number of ref. set freq. in 2nd stopband
   if q1 == 0
      q1 = 1;
   elseif q1 == q
      q1 = q - 1;
   end
   
   q2 = q - q1;
   
   if q1 == 1;
      rs1 = ws1;
   else
      rs1 = [0:q1-1]'*(ws1/(q1-1));
   end
   if q2 == 1
      rs2 = ws2;
   else
      rs2 = [0:q2-1]'*(pi-ws2)/(q2-1)+ws2;
   end
   rs = [rs1; rs2];
   
   Y1 = [ls*(1-(-1).^(q1:-1:1))/2 + us*((-1).^(q1:-1:1)+1)/2]';
   Y2 = [ls*(1-(-1).^(1:q2))/2 + us*((-1).^(1:q2)+1)/2]';
   Y = [Y1; Y2];
   
   n  = 0:q-1;
   Z   = zeros(2*(g-q)+1,1);
   % A1  = (-1)^(L/2) * (sin(w/2-wp/2).*sin(w/2+wp/2)).^(L/2);
   % A1r = (-1)^(L/2) * (sin(rs/2-wp/2).*sin(rs/2+wp/2)).^(L/2);
   A1  = (-1)^(L/2) * ((cos(wp)-cos(w))/2).^(L/2);
   A1r = (-1)^(L/2) * ((cos(wp)-cos(rs))/2).^(L/2);
   it = 0;
   while 1 & (it < 35)
      if length(rs) ~= length(n)
         rs ,n
         error('Filter design code error')
      end
      a2 = cos(rs*n)\[(Y-1)./A1r];
      A2 = real(fft([a2(1);a2(2:q)/2;Z;a2(q:-1:2)/2])); A2 = A2(1:g+1);
      A  = 1 + A1.*A2;
      %   figure(1), plot(w/pi,A),
      %        hold on, plot(rs/pi,Y,'o'), hold off, % pause(.1)
      %   figure(2), plot(w/pi,20*log10(abs(A))),
      %        hold on, plot(rs/pi,20*log10(abs(Y)),'o'), hold off, % pause(.1)
      %        pause(.5)
      ri = sort([localmax(A); localmax(-A)]);
      lri = length(ri);
      % ri(1:length(ri)-q) = [];
      if lri ~= q+1
         if abs(A(ri(lri))-A(ri(lri-1))) < abs(A(ri(1))-A(ri(2)))
            ri(lri) = [];
         else
            ri(1) = [];
         end
      end
      rs = (ri-1)*pi/g;
      [temp, k] = min(abs(rs - wp)); rs(k) = [];
      q1 = sum(rs < wp);
      q2 = sum(rs > wp);
      
      Y1 = [ls*(1-(-1).^(q1:-1:1))/2 + us*((-1).^(q1:-1:1)+1)/2]';
      Y2 = [ls*(1-(-1).^(1:q2))/2 + us*((-1).^(1:q2)+1)/2]';
      Y = [Y1; Y2];
      
      % rs = refine2(a2,L/2,rs);
      % A1r = (-1)^(L/2) * (sin(rs/2-wp/2).*sin(rs/2+wp/2)).^(L/2);
      A1r = (-1)^(L/2) * ((cos(wp)-cos(rs))/2).^(L/2);
      Ar = 1 + (cos(rs*n)*a2) .* A1r;
      Err = max([max(Ar)-us, ls-min(Ar)]);
      %   fprintf(1,'    Err = %18.15f\n',Err);
      if Err < SN, break, end
      it = it + 1;
   end
   
   h2 = [a2(q:-1:2)/2; a2(1); a2(2:q)/2];
   h = h2;
   for k = 1:L/2
      h = conv(h,[1 -2*cos(wp) 1]')/4;
   end
   h((N+1)/2) = h((N+1)/2) + 1;
   
end

% ========================================================

function C = chebpoly(n)

% C = cheb_poly(n)
% Chebychev polynomial
%

if n == 0
   C = 1;
elseif n == 1
   C = [1 0];
else
   A = 1;
   B = [1 0];
   for k = 2:n
      C = 2*[B 0] - [0 0 A];
      A = B;
      B = C;
   end
end

% =========================================================

function b = cos2x(a)
%
% converts the cos polynomial, 
%       a(f) = a(1) + a(2)*cos(w) + ... + a(m+1)*cos(m*w)
% over [0,1/2],
% to the polynomial
%       b(x) = b(m+1) + b(m)*x + ... + b(1)*x^m
% over [-1,1]
%
% the transformation is : x = cos(w)
%
% (x == 1) and (w == 0) are mapped together
% (x == -1) and (w == pi) are mapped together
%

m = length(a)-1;
c = a(1);
for k = 2:m+1
   c = add_poly(c,a(k)*chebpoly(k-1));
end
b = zeros(size(a));
b(:) = c;

% =============================================

function a = h2cos(h)
% a = h2cos(h);

N = length(h);
if sum(abs(h-h(N:-1:1))) > 0.00001
   disp('for symmetric h only')
   a = [];
   return
end
if rem(N,2)==1
   % N even
   a = 2*h((N+1)/2:N);
   a(1) = a(1)/2;
else
   disp('for odd length only')
   return
end

% ==================================================

function x = h2x(h)
%
% x = h2x(h)
%
x = cos2x(h2cos(h));

% ====================================================

function k = localmax(x)
% k = localmax(x)
% finds location of local maxima
%
s = size(x); x = [x(:)].'; N = length(x);
b1 = x(1:N-1)<=x(2:N); b2 = x(1:N-1)>x(2:N);
k = find(b1(1:N-2)&b2(2:N-1))+1;
if x(1)>x(2), k = [k, 1]; end
if x(N)>x(N-1), k = [k, N]; end
k = sort(k); if s(2) == 1, k = k'; end

% =============================================

function p = rlz(p)

if isempty(p)
   break
end

while p(1) == 0 
   p(1) = [];
   if isempty(p)
      break
   end
end

% =================================================
