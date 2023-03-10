function B = calcInternalForceMatrix(opts)

%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   beta : thin plate energy (second order)
%   gamma : Step Size (Time)

aa = opts.alpha * opts.dt/opts.ds/opts.ds;
bb = opts.beta  * opts.dt/opts.ds/opts.ds;

b(1) = bb;

b(2) = -(aa + 4*bb);

b(3) = 2*aa + 6*bb;

b(4) = b(2);

b(5) = b(1);

idMatrix = eye(opts.numPoints);

A =     b(1)*circshift(idMatrix,  2);
A = A + b(2)*circshift(idMatrix,  1);
A = A + b(3)*circshift(idMatrix,  0);
A = A + b(4)*circshift(idMatrix, -1);
A = A + b(5)*circshift(idMatrix, -2);

A(1:2, :)       = 0;
A(end-1:end, :) = 0;

A(1,1) = 0 + aa;
A(1,2) =    -aa;

A(2,1) =      -aa;
A(2,2) = 0 + 2*aa;
A(2,3) =      -aa;

A(end,   end)   = 0 + aa;
A(end,   end-1) =    -aa;

A(end-1, end)   =      -aa;
A(end-1, end-1) = 0 + 2*aa;
A(end-1, end-2) =      -aa;


B = inv(A + opts.dt * idMatrix);

