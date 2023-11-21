# irbleigs
A MATLAB program for computing a few eigenvalues and associated eigenvectors located anywhere in spectrum of a large sparse Hermitian matrix. 

% IRBLEIGS will find a few eigenvalues and eigenvectors for either the 
% standard eigenvalue problem A*x = lambda*x or the generalized eigenvalue 
% problem A*x = lambda*M*x, where A is a sparse Hermitian matrix and the 
% matrix M is positive definite.
%
% [V,D,PRGINF] = IRBLEIGS(A,OPTIONS) 
% [V,D,PRGINF] = IRBLEIGS(A,M,OPTIONS)
% [V,D,PRGINF] = IRBLEIGS('Afunc',N,OPTIONS) 
% [V,D,PRGINF] = IRBLEIGS('Afunc',N,M,OPTIONS)
%
% The first input argument into IRBLEIGS can be a numeric matrix A or an M-file 
% ('Afunc') that computes the product A*X, where X is a (N x blsz) matrix. If A is 
% passed as an M-file then the (N x blsz) matrix X is the first input argument, the 
% second input argument is N, (the size of the matrix A), and the third input argument 
% is blsz, i.e. Afunc(X,N,blsz). For the generalized eigenvalue problem the matrix M is 
% positive definite and passed only as a numeric matrix. IRBLEIGS will compute the 
% Cholesky factorization of the matrix M by calling the internal MATLAB function CHOL. 
% M may be a dense or sparse matrix. In all the implementations IRBLEIGS(A,...) can
% be replaced with IRBLEIGS('Afunc',N,...).
%
% NOTE: If the M-file 'Afunc' requires additional input parameters, e.g. a data structure, 
%       use the option FUNPAR to pass any additional parameters to your function (X,N,blsz,FUNPAR).
% 
% OUTPUT OPTIONS:
% ---------------
%
% I.)   IRBLEIGS(A) or IRBLEIGS(A,M)
%       Displays the desired eigenvalues.    
%
% II.)  D = IRBLEIGS(A) or D = IRBLEIGS(A,M)
%       Returns the desired eigenvalues in the vector D. 
%
% III.) [V,D] = IRBLEIGS(A) or [V,D] = IRBLEIGS(A,M)  
%       D is a diagonal matrix that contains the desired eigenvalues along the 
%       diagonal and the matrix V contains the corresponding eigenvectors, such 
%       that A*V = V*D or A*V = M*V*D. If IRBLEIGS reaches the maximum number of
%       iterations before convergence then V = [] and D = []. Use output option
%       IV.) to get approximations to the Ritz pairs.
%
% IV.)  [V,D,PRGINF] = IRBLEIGS(A) or [V,D,PRGINF] = IRBLEIGS(A,M)
%       This option returns the same as (III) plus a two dimensional array PRGINF 
%       that reports if the algorithm converges and the number of matrix vector 
%       products. If PRGINF(1) = 0 then this implies normal return: all eigenvalues have 
%       converged. If PRGINF(1) = 1 then the maximum number of iterations have been 
%       reached before all desired eigenvalues have converged. PRGINF(2) contains the 
%       number of matrix vector products used by the code. If the maximum number of 
%       iterations are reached (PRGINF(1) = 1), then the matrices V and D contain any 
%       eigenpairs that have converged plus the last Ritz pair approximation for the 
%       eigenpairs that have not converged.
%
% INPUT OPTIONS:
% --------------
%                                   
%       ... = IRBLEIGS(A,OPTS) or  ... = IRBLEIGS(A,M,OPTS)
%       OPTS is a structure containing input parameters. The input parameters can
%       be given in any order. The structure OPTS may contain some or all of the 
%       following input parameters. The string for the input parameters can contain
%       upper or lower case characters. 
%       
%  INPUT PARAMETER      DESCRIPTION                     
%   
%  OPTS.BLSZ         Block size of the Lanczos tridiagonal matrix.        
%                    DEFAULT VALUE    BLSZ = 3
%                        
%  OPTS.CHOLM        Indicates if the Cholesky factorization of the matrix M is available. If
%                    the Cholesky factorization matrix R is available then set CHOLM = 1 and 
%                    replace the input matrix M with R where M = R'*R.
%                    DEFAULT VALUE   CHOLM = 0
%
%  OPTS.PERMM        Permutation vector for the Cholesky factorization of M(PERMM,PERMM). 
%                    When the input matrix M is replaced with R where M(PERMM,PERMM)=R'*R
%                    then the vector PERMM is the permutation vector. 
%                    DEFAULT VALUE   PERMM=1:N
%
%  OPTS.DISPR        Indicates if K Ritz values and residuals are to be displayed on each 
%                    iteration. Set positive to display the Ritz values and residuals on 
%                    each iteration.
%                    DEFAULT VALUE   DISPR = 0 
%
%  OPTS.EIGVEC       A matrix of converged eigenvectors.        
%                    DEFAULT VALUE  EIGVEC = []
%                                                      
%  OPTS.ENDPT        Three letter string specifying the location of the interior end- 
%                    points for the dampening interval(s). 
%                    'FLT' - Let the interior endpoints float.                  
%                    'MON' - Interior endpoints are chosen so that the size of the 
%                            dampening interval is increasing. This creates a nested
%                            sequence of intervals. The interior endpoint will approach 
%                            the closest Ritz value in the undampened part of the spectrum 
%                            to the dampening interval.
%                    DEFAULT VALUE   ENDPT = 'MON' (If SIGMA = 'LE' or 'SE'.)
%                    DEFAULT VALUE   ENDPT = 'FLT' (If SIGMA = a numeric value NVAL.)
%
%  OPTS.FUNPAR       If A is passed as a M-file then FUNPAR contains any additional parameters 
%                    that the M-file requires in order to compute the matrix vector product. 
%                    FUNPAR can be passed as any type, numeric, character, data structure, etc. The
%                    M-file must take the input parameters in the following order (X,n,blsz,FUNPAR).
%                    DEFAULT VALUE  FUNPAR =[]                   
%
%  OPTS.K            Number of desired eigenvalues.             
%                    DEFAULT VALUE  K = 3
%
%  OPTS.MAXIT        Maximum number of iterations, i.e. maximum number of block Lanczos restarts.                           
%                    DEFAULT VALUE  MAXIT = 100
%
%  OPTS.MAXDPOL      Numeric value indicating the maximum degree of the dampening 
%                    polynomial allowed.  
%                    DEFAULT VALUE   MAXDPOL = 200     (If SIGMA = 'LE' or 'SE'.)
%                    DEFAULT VALUE   MAXDPOL = N       (If SIGMA = a numeric value NVAL.)
%
%  OPTS.NBLS         Number of blocks in the Lanczos tridiagonal matrix. The program may increase
%                    NBLS to ensure certain requirements in [1] are satisfied. A warning message
%                    will be displayed if NBLS increases.                          
%                    DEFAULT VALUE    NBLS = 3
%
%  OPTS.SIGMA        Two letter string or numeric value specifying the location 
%                    of the desired eigenvalues.            
%                    'SE'  Smallest Real eigenvalues.                
%                    'LE'  Largest Real eigenvalues.                 
%                    NVAL  A numeric value. The program searches for the K closest
%                          eigenvalues to the numeric value NVAL. 
%                    DEFAULT VALUE   SIGMA = 'LE'
%                                                         
%  OPTS.SIZINT       Size of the dampening interval. Value of 1 indicates consecutive
%                    Ritz values are used to determine the endpoints of the dampening
%                    interval. Value of 2 indicates endpoints are chosen from Ritz
%                    values that are seprated by a single Ritz value. A value of 3
%                    indicates endpoints are chosen from Ritz values that are seprated 
%                    by two Ritz values. Etc. The minimum value is 1 and the maximum 
%                    value is (NBLS-1)*BLSZ-K. The program may modify SIZINT without
%                    notification to ensure certain requirements in [1] are satisfied. 
%                    DEFAULT VALUE    SIZINT = 1
%                                 
%  OPTS.TOL          Tolerance used for convergence. Convergence is determined when             
%                    || Ax - lambda*x ||_2 <= TOL*||A||_2. ||A||_2 is approximated by
%                    largest absolute Ritz value.  
%                    DEFAULT VALUE    TOL = 1d-6
%                                                              
%  OPTS.V0           A matrix of starting vectors.       
%                    DEFAULT VALUE  V0 = randn
%
%  OPTS.ZERTYP       Two letter string to indicate which type of zeros to apply.              
%                    'WL' - Weighted fast Leja points. The weight functions are used to help
%                           increase convergence.  
%                    'ML' - Mapped fast Leja points. Fast Leja points are computed on [-2,2] 
%                           and mapped to the dampening interval. This option is not available
%                           when sigma is a numeric value NVAL.
%                    DEFAULT VALUE  ZERTYP = 'ML' (If SIGMA = 'LE' or 'SE'.)
%

%  DATE MODIFIED: 04/20/2004
%  VER:  1.0

%  AUTHORS:
%  James Baglama     University of Rhode Island, E-mail: jbaglama@math.uri.edu
%  Daniela Calvetti  Case Western Reserve University,  E-mail: dxc57@po.cwru.edu
%  Lothar Reichel    Kent State University, E-mail: reichel@mcs.kent.edu
   
% REFERENCES:
%   1.) "IRBL: An Implicitly Restarted Block Lanczos Method for large-scale Hermitian 
%        eigenproblems", J. Baglama, D. Calvetti, and L. Reichel, SIAM J. Sci. Comput., 
%        in press, 2003.
%   2.) "irbleigs: A MATLAB program for computing a few eigenpairs of a large sparse
%        Hermitian matrix", J. Baglama, D. Calvetti, and L. Reichel, Technical Report
%        submitted for publication (2001).
%   3.) "Dealing With Linear Dependence during the Iterations of the Restarted
%        Block Lanczos Methods", J. Baglama, Num. Algs., 25, (2000) pp. 23-36. 
%   4.) "Fast Leja Points", J. Baglama, D. Calvetti, and L. Reichel, ETNA, 
%        Vol. 7 (1998), pp. 124-140.
%   5.) "Computation of a few close eigenvalues of a large matrix with 
%        application to liquid crystal modeling", J. Baglama, D. Calvetti, 
%        L. Reichel, and A. Ruttan, J. of Comp. Phys., 146 (1998), pp. 203-226.
%   6.) "Iterative Methods for the Computation of a Few Eigenvalues of a Large 
%        Symmetric Matrix", J. Baglama, D. Calvetti, and L. Reichel, BIT, 36 
%        (1996), pp. 400-421.
