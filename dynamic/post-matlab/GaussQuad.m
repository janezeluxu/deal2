function [xx, ww] = GaussQuad(n, a, b)
   narginchk(1, 3);

   % Assign default values to missing arguments.
   switch nargin
      case 1                    % GAUSSQUAD(N)
         b = 1;
         a = -1;
      case 2                    % GAUSSQUAD(N, C)
         b = a;
         a = 0;
   end

   u = 1 : n-1;
   u = u ./ sqrt(4*u.^2 - 1);

   % Same as A = diag(u, -1) + diag(u, 1), but faster (no addition).
   A = zeros(n, n);
   A( 2 : n+1 : n*(n-1) ) = u;
   A( n+1 : n+1 : n^2-1 ) = u;

   % Find the base points X and weight factors W for the interval [-1,1].
   [v, x] = eig(A);
   [x, k] = sort(diag(x));
   w = 2 * v(1,k)'.^2;

   % Linearly transform from [-1,1] to [a,b].
   x = (b - a) / 2 * x + (a + b) / 2;   % transform base points X
   w = (b - a) / 2 * w;                 % adjust weigths

   % If output arguments are given, return output and exit.
   if nargout
      xx = x;
      ww = w;
      return
   end
end