function [ ynew, improvement, msgTR ] = ...
          bcdfo_find_new_yj( QZ, RZ, Y, j, Delta, eps_L, xbase, lSolver, ...
          whichmodel, xl, xu, indfree, stratLam, scale, shift_Y )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Computes a point which is best to replace yj, the j-th (j>1) column of Y (the
%  base for the j-th polynomial) in a ball of radius Delta centered at the
%  first column of Y.  This is achieved by maximizing the absolute value of
%  the j-th Lagrange polynomial in that ball.  

%  For conventions on how polynomals are represented, see the documentation of
%  evalZ.

%  INPUT:

%  QZ          : the Q matrix of the QR decomposition of Z(Y)
%  RZ          : the R matrix of the QR decomposition of Z(Y)
%  Y           : a matrix whose columns contain the current interpolation points
%  j           : the index of the interpolation point one wishes to replace (j > 1)
%  Delta       : the radius of the ball centered at Y(:,1) in which the
%                replacement point must be found
%  eps_L       : the relative accuracy on the trust-region constraint for
%                maximization of the Lagrange polynomials
%  xbase       : the current base point
%  lSolver     : linear solver used for the minimization of the model
%  whichmodel  : kind of model/Lagrange polynomial to compute
%  scale       : the current interpolation set scaling
%  shift_Y     : 0 if no shift in interpolation points, 1 otherwise

%  OUTPUT:

%  ynew        : the best replacement for Y(:,j)
%  improvement : the improvement in poisedness obtained by the update, which
%                is equal to |L_j(new y)|. If this value is smaller than the
%                threshold input parameter, L and X are unchanged by the
%                procedure. 

%  PROGRAMMING: Ph. Toint, February 2009. (This version 22 VI 2009)

%  USES: bcdfo_gradP, bcdfo_hessP, bcdfo_solve_TR_MS

%  TEST:
%  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
%  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 0 );
%  [ ynew, improvement ] = bcdfo_find_new_yj_bc( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
%       whichmodel, [-1e10,-1e10], [1e10, 1e10], [1:2], 1, scale, 0 )
%  should give
%  ynew =
%
%    3.2808
%   -0.9598
%
%  improvement =
%
%  314.8825
%
%  the same must be obtained by the shifted and scaled version:
%  Y = [ 3 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; whichmodel = 0;
%  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y , whichmodel, 1 );
%  [ ynew, improvement ] = bcdfo_find_new_yj_bc( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, ...
%       whichmodel, [-1e10,-1e10], [1e10, 1e10], [1:2], 1, scale, 1 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose     = 0;              % 1 for debug
n           = size(Y,1);
ynew        = zeros(1,n);
improvement = 0;

if ( verbose )
   disp('--------- enter find_new_yj_bc ')
end

if ( j < 2 )           % never attempt to replace the current iterate.   
   return;
end

%  Get the j-th Lagrange polynomial 

Lj = bcdfo_computeLj( QZ, RZ, j, Y, whichmodel, scale, shift_Y );

if ( length( find( isnan(Lj) ) ) ~= 0 | length( find( ~isreal(Lj) ) ) ~= 0 | ...
     length( find( isinf(Lj) ) ) ~= 0)

   msgTR = 'Error0: Lagrange polynomial contains NaN or Inf or nonreal components!!';
   if ( verbose )
      disp(msgTR);
   end
   return
end

%  Maximize Lj in a larger 2-norm TR if using infty-norm in the local solver (CG)

if ( lSolver == 2 )
	Delta = sqrt(n)*Delta;
end

lb = xl(indfree) - Y(:,1);
ub = xu(indfree) - Y(:,1);

%  Get the polynomial's gradient and Hessian at the current iterate.

if ( shift_Y )

   %  When shifted, the origin in the scaled variables corresponds 
   %  to Y(:,1) in the unscaled space

   g  = bcdfo_gradP( Lj, zeros(n,1), xbase, scale, 0 );
   H  = bcdfo_hessP( Lj, zeros(n,1), xbase, scale, 0 );
   
   %  Minimize this polynomial and its opposite.
   
   [ pstep, lambda, norms, pvalue, gplus, nfact, neigd, msgTR ] = ...
         bcdfo_solve_TR_MS_bc( g, H, lb.*scale(2), ub.*scale(2), Delta*scale(2), ...
            eps_L, stratLam );
   pstep = pstep / scale(2);
   [ mstep, lambda, norms, mvalue, gplus, nfact, neigd, msgTR ] = ...
         bcdfo_solve_TR_MS_bc( -g, -H, lb.*scale(2), ub.*scale(2), Delta*scale(2), ...
            eps_L, stratLam );
   mstep = mstep / scale(2);
 
else

   %  When no shift occurs, the current iterate is Y(:,1)

   g  = bcdfo_gradP( Lj, Y(:,1), xbase, scale, 0 );
   H  = bcdfo_hessP( Lj, Y(:,1), xbase, scale, 0 );

   %  Minimize this polynomial and its opposite.

   [ pstep, lambda, norms, pvalue, gplus, nfact, neigd, msgTR ] = ...
         bcdfo_solve_TR_MS_bc( g, H, lb, ub, Delta, eps_L, stratLam );
   [ mstep, lambda, norms, mvalue, gplus, nfact, neigd, msgTR ] = ...
         bcdfo_solve_TR_MS_bc( -g, -H, lb, ub, Delta, eps_L, stratLam );
end


if ( verbose )
   disp([ ' === find_new_yj_bc: j = ' int2str(j), ' positive value = ', num2str(pvalue),' step:' ])
   pstep'
   disp([ ' === find_new_yj_bc: j = ' int2str(j), ' negative value = ', num2str(mvalue),' step:' ])
   mstep'
end

%  Select the maximum in absolute value.

if ( mvalue < pvalue )
   improvement = abs( mvalue );
   ynew        = Y(:,1) + mstep;
else
   improvement = abs( pvalue );
   ynew        = Y(:,1) + pstep;
end
if ( verbose )
   disp('--------- exit find_new_yj_bc ')
end
