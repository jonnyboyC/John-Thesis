// JGradient.c
// Gradient along a dimension
// Y = JGradient(X, Dim, Spacing, Method, bnd_idx)
// INPUT:
//   X:   Real DOUBLE array.
//   Spacing: Scalar or vector of the length SIZE(X, Dim).
//        A scalar value is the distance between all points, while a vector
//        contains all coordinates, such that DIFF(Spacing) are the distances.
//        For equally spaced input a scalar Spacing is much faster.
//        Optional, default: 1.0
//   Dim: Dimension to operate on.
//        Optional, default: [] (1st non-singelton dimension).
//   Method: String, order of the applied method for unevenly spaced X:
//        '1stOrder', faster centered differences as in Matlab's GRADIENT.
//        '2ndOrder', 2nd order accurate centered differences.
//        On the edges forward and backward difference are used.
//        Optional, default: '1stOrder'.
//   Bnd_idx: Matrix of values indicating where boundaries exist in the flow
//
// OUTPUT:
//   Y:   Gradient of X, same size as X.
//
// EXAMPLES:
//   t = cumsum(rand(1, 100)) + 0.01;  t = 2*pi * t ./ max(t);
//   x = sin(t);
//   dx1 = JGradient(x, t, 2, '1stOrder');
//   dx2 = JGradient(x, t, 2, '2ndOrder');
//   dx  = cos(t);          % Analytic solution
//   h = plot(t, dx, t, dx1, 'or', t, dx2, 'og');  axis('tight');
//   title('cos(x) and JGradient(sin(x))');
//   legend(h, {'analytic', '1st order', '2nd order'}, 'location', 'best');
//
// NOTES:
// - There are a lot of other derivation tools in the FEX. This function is
//   faster, e.g. 25% faster than dqdt.c and 10 to 16 times faster than Matlab's
//   GRADIENT. In addition it works with multi-dim arrays, on a speicifc
//   dimension only and can use a 2nd order method for unevenly spaced data.
// - This function does not use temporary memory for evenly spaced data and if
//   a single vector is processed. Otherwise the 1st-order method needs one and
//   the 2nd-order method 3 temporary vectors of the length of the processed
//   dimension.
// - Matlab's GRADIENT processes all dimensions ever, while JGradient operates on
//   the specified dimension only.
// - 1st order centered difference:
//     y(i) = (x(i+1) - x(i-1) / (s(i+1) - s(i-1))
// - 2nd order centered difference:
//     y(i) = ((x(i+1) * (s(i)-s(i-1)) / (s(i+1)-s(i))) -
//             (x(i-1) * (s(i+1)-s(i)) / (s(i)-s(i-1)))) / (s(i+1)-s(i-1))
//            + x(i) * (1.0 / (s(i)-s(i-1)) - 1.0 / (s(i+1)-s(i)))
//   For evenly spaced X, both methods reply equal values.
//
// COMPILE:
//   mex -O JGradient.c
// Consider C99 comments on Linux:
//   mex -O CFLAGS="\$CFLAGS -std=c99" JGradient.c
// Pre-compiled Mex: http://www.n-simon.de/mex
// Run the unit test uTest_JGradient after compiling.
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
// Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
// Author: Jan Simon, Heidelberg, (C) 2011 matlab.THISYEAR(a)nMINUSsimon.de
//
// See also GRADIENT, DIFF.
// FEX: central_diff (#12 Robert A. Canfield)
//      derivative (#28920, Scott McKinney)
//      movingslope (#16997, John D'Errico)
//      diffxy (#29312, Darren Rowland)
//      dqdt (#11965, Geoff Wawrzyniak)

// Todo: SSE2 instructions
//       Newton polynomials for 2nd order edges
//       Multi-threading

/*
% $JRev: R0d V:004 Sum:sHhGcnzMMNaA Date:02-Jan-2008 17:46:05 $
% $License: BSD $
% $File: Tools\Mex\Source\JGradient.c $
% History:
% 001: 30-Dec-2010 22:42, First version published under BSD license.
*/

// Includes:
#include "mex.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T               // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// There is an undocumented method to create a shared data copy. This is much
// faster, if the replied object is not changed, because it does not duplicate
// the contents of the array in the memory.
mxArray *mxCreateSharedDataCopy(const mxArray *mx);
#define COPY_ARRAY mxCreateSharedDataCopy
// #define COPY_ARRAY mxDuplicateArray    // slower, but documented

// Disable the /fp:precise flag to increase the speed on MSVC compiler:
#ifdef _MSC_VER
#pragma float_control(except, off)    // disable exception semantics
#pragma float_control(precise, off)   // disable precise semantics
#pragma fp_contract(on)               // enable contractions
// #pragma fenv_access(off)           // disable fpu environment sensitivity
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_HEAD "*** JGradient[mex]: "
#define ERR_ID   "JChabot:JGradient:"

// Prototypes: ---------Added By John Chabot-----------------------------------
void WrapSpaceN(double *X, const mwSize Step, const mwSize nX,
    const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y);
void CoreDim1Space1(double *X, const mwSize Step, const mwSize nX,
	const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y);
void CoreDimNSpace1(double *X, const mwSize Step, const mwSize nX,
	const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y);
void CoreDim1SpaceN(double *X, const mwSize Step, const mwSize nX,
	const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y);
void CoreDimNSpaceN(double *X, const mwSize Step, const mwSize nX,
	const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y);
void GetFactor(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
void GetFactor2(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
void GetFactor3(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
void GetFactor4(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
// ------------------------------------------------------------------------

mwSize FirstNonSingeltonDim(const mwSize Xndim, const mwSize *Xdim);
mwSize GetStep(const mwSize *Xdim, const mwSize N);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X, *Y, Nd, *Space, UnitSpace = 1.0;
    int8_T *Bnd_idx;
    mwSize nX, ndimX, N, nDX, Step, nSpace, nBnd_idx;
    const mwSize *dimX;
    int Order2, Order4;
  
    // Check number and type of inputs and outputs: ------------------------------
    if (nrhs == 0 || nrhs > 4) {
        mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                        ERR_HEAD "1 or 4 inputs required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                        ERR_HEAD "1 output allowed.");
    }
  
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput1",
                        ERR_HEAD "Input must be a full real double array.");
    }
  
    // Pointers and dimension to input array: ------------------------------------
    X     = mxGetPr(prhs[0]);
    nX    = mxGetNumberOfElements(prhs[0]);
    ndimX = mxGetNumberOfDimensions(prhs[0]);
    dimX  = mxGetDimensions(prhs[0]);
  
    // Return fast on empty input matrix:
    if (nX == 0) {
        plhs[0] = COPY_ARRAY(prhs[0]);
        return;
    }
  
    // Check 2th input, get spacing if defined: --------------------------------------------------
    if (nrhs < 2) {  // No 2nd input defined - scalar unit spacing:
        nSpace = 1;
        Space  = &UnitSpace;
     
    } else {         // Get pointer to spacing vector:
        if (!mxIsDouble(prhs[1])) {
            mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput2",
                        ERR_HEAD "2nd input [Spacing] must be a DOUBLE.");
        }
        Space  = mxGetPr(prhs[1]);
        nSpace = mxGetNumberOfElements(prhs[1]);
        if (nSpace == 0) {
            nSpace = 1;
            Space  = &UnitSpace;
        }
    }
  
    // Check 3rd input, determine dimension to operate on: ----------------------------------------
    if (nrhs < 3) {
        N    = FirstNonSingeltonDim(ndimX, dimX);  // Zero based
        Step = 1;
        nDX  = dimX[N];
     
    } else if (mxIsNumeric(prhs[2]))  {  // 3rd input used:
        switch (mxGetNumberOfElements(prhs[2])) {
            case 0:  // Use 1st non-singelton dim if 3rd input is []:
                N    = FirstNonSingeltonDim(ndimX, dimX);
                Step = 1;
                nDX  = dimX[N];
                break;
           
            case 1:  // Numerical scalar:
                Nd = mxGetScalar(prhs[2]);
                N  = (mwSize) Nd - 1;
                if (Nd < 1.0 || Nd != floor(Nd)) {
                    mexErrMsgIdAndTxt(ERR_ID   "BadValueInput3",
                        ERR_HEAD "Dimension must be a positive integer scalar.");
                }
           
                if (N < ndimX) {
                    Step = GetStep(dimX, N);
                    nDX  = dimX[N];
                } else {
                    // Treat imaginated trailing dimensions as singelton, as usual in
                    // Matlab:
                    Step = nX;
                    nDX  = 1;
                }
                break;
           
                default:
                    mexErrMsgIdAndTxt(ERR_ID   "BadSizeInput3",
                         ERR_HEAD "3rd input [Dim] must be scalar index.");
        }
     
    } else {  // 2nd input is not numeric:
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput3",
                ERR_HEAD "3rd input must be scalar index.");
    }
  
    // Check matching sizes of X and Spacing:
    if (nSpace != 1 && nSpace != nDX) {
        mexErrMsgIdAndTxt(ERR_ID   "BadSizeInput2",
                ERR_HEAD "2nd input [Spacing] does not match the dimensions.");
    }
  
    // Check 4th input: ----------------------------------------------------------
    // -----------------------------TODO-----------------------------------------
    if (nrhs < 4) {
        Bnd_idx = &nX;
        Bnd_idx = (int8_T *) mxMalloc(nX, sizeof(int8_T));
        for(int i = 0; i < nX; ++i) {
            *(Bnd_idx + 1) = 1;
        }
     
    } else {
        if (!mxIsDouble(prhs[1])) {
            mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput5",
                  ERR_HEAD "5th input [Boundary Index] must be a DOUBLE.");
        }
        Bnd_idx = mxGetPr(prhs[3]);
        nBnd_idx = mxGetNumberOfElements(prhs[3]);
        if (nBnd_idx != nX) {
            Bnd_idx = (int8_T *) mxMalloc(nX, sizeof(int8_T));
            for(int i = 0; i < nX; ++i) {
                *(Bnd_idx + i) = 1;
            }
        }
    }
    // -----------------------------TODO-----------------------------------------
  
    // Create output matrix: -----------------------------------------------------
    plhs[0] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
    Y      = mxGetPr(plhs[0]);

    // Reply ZEROS, if the length of the processed dimension is 1:
    if (nDX == 1) {
        return;
    }
  
    // Calculate the gradient: ---------------------------------------------------
    if (nSpace == 1) {         // Scalar spacing
        if (Step == 1) {        // Operate on 1st dimension
			CoreDim1Space1(X, nX, nDX, *Space, Bnd_idx, Y);
        } else {                // Step >= 1, operate on any dimension
            CoreDimNSpace1(X, Step, nX, nDX, *Space, Y);
        }
        
    //--------------------------- TODO -------------------------------------------
    } else {                   // Spacing defined as vector, 1st order method:
        if (nX == nDX) {        // Single vector only - dynamic spacing factors:
			CoreDim1SpaceN(X, nX, nDX, Space, Bnd_idx, Y);
        } else {
			WrapSpaceN(X, Step, nX, nDX, Space, Bnd_idx, Y);
        }
    }
    
    return;
}

// Subroutines: ================================================================
mwSize FirstNonSingeltonDim(const mwSize Xndim, const mwSize *Xdim)
{
    // Get first non-singelton dimension - zero based.
    
    mwSize N;
    
    for (N = 0; N < Xndim; N++) {
        if (Xdim[N] != 1) {
            return (N);
        }
    }
    
    return (0);  // Use the first dimension if all dims are 1
}

// =============================================================================
mwSize GetStep(const mwSize *Xdim, const mwSize N)
{
    // Get step size between elements of a subvector in the N'th dimension.
    // This is the product of the leading dimensions.
    
    const mwSize *XdimEnd, *XdimP;
    mwSize       Step;
    
    Step    = 1;
    XdimEnd = Xdim + N;
    for (XdimP = Xdim; XdimP < XdimEnd; Step *= *XdimP++) ; // empty loop
    
    return (Step);
}

// =============================================================================
void CoreDim1Space1(double *X, const mwSize nX, const mwSize nDX, double Space,
                    int8_T *Bnd_idx, double *Y)
{
    // Operate on first dimension, scalar spacing, 1st order method.
    
    double x0, x1, x2, *Xf, *Xc, fac1, fac2;
    
	// Will pass in factors array
    
    Xf = X + nX;                 // End of input array
    while (X < Xf) {
		switch (method) {
			case 1:
			case 2:
			case 3:
			case 4:
			//... will pass in method pointer
		}
    }
    return;
}

// =============================================================================
void CoreDim1SpaceN(double *X, const mwSize nX, const mwSize nDX,
	double *Space, int8_T *Bnd_idx double *Y)
{
	// Operate on the first dimension, spacing is a vector, order 1 method.
	// The spacing factors are calculated dynamically. This is efficient for
	// a single vector, but slower for matrices. No temporary vector is needed.

	double x0, x1, x2, *Xf, *Xg, *Sp, s0, s1, s2;

	// Will pass in factors array

	Xf = X + nX;
	while (X < Xf) {
		switch (method) {
		case 1:
		case 2:
		case 3:
		case 4:
			//... will pass in method pointer
		}
	}
	return;
}

// =============================================================================
void CoreDimNSpace1(double *X, const mwSize Step, const mwSize nX,
                    const mwSize nDX, double Space, int8_T *Bnd_idx, double *Y)
{
    // Operate on any dimension, scalar spacing, 1st order method.
    // Column oriented approach: Process contiguous memory blocks of input and
    // output.
    
    double *Xf, *X1, *X2, *Xc, fac1, fac2;
    mwSize nDXStep;
    
	// Will pass in factors array
    
    // Distance between first and last element of X in specified dim:
    nDXStep = nDX * Step;
    
	Xf = X + nX;                 // End of input array
	while (X < Xf) {
		switch (method) {
		case 1:
		case 2:
		case 3:
		case 4:
			//... will pass in method pointer
		}
	}    
    return;
}

// =============================================================================
void CoreDimNSpaceN(double *X, const mwSize Step, const mwSize nX,
					const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y)
{
	// Operate on any dimension, scalar spacing, 1st order method.
	// Column oriented approach: Process contiguous memory blocks of input and
	// output.

	double *Xf, *X1, *X2, *Xc, fac1, fac2;
	mwSize nDXStep;

	// Will pass in factors array

	// Distance between first and last element of X in specified dim:
	nDXStep = nDX * Step;

	Xf = X + nX;                 // End of input array
	while (X < Xf) {
		switch (Method) {
		case 1:
		case 2:
		case 3:
		case 4:
			//... will pass in method pointer
		}
	}
	return;
}

// =============================================================================
void WrapSpace1(double *X, const mwSize Step, const mwSize nX,
	const mwSize nDX, double Space, int8_T *Bnd_idx, double *Y)
{
	// Call different methods depending of the dimensions ofthe input.
	// X has more than one vector. Therefore it is cheaper to calculate the
	// spacing factors once only.

	double *A, *B, *C, *D, *E;
	int8_T *Method;

	// Predetermine Method
	Method = (int8_T *)mxMalloc(nX *sizeof(int8_T));

	if (A == NULL || B == NULL || C == NULL) {
		mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
			ERR_HEAD "No memory for Factor vectors.");
	}

	GetFactor(Space, nDX, nX, A, B, C, D, E, Bnd_idx, Method);
	GetMethod(nDX, nX, Bnd_idx);

	if (Step == 1) {        // Operate on first dimension:
		//CoreDim1FactorNOrder2(X, nX, nDX, A, B, C, Y);
		CoreDim1SpaceNOrder2(X, nX, nDX, Space, Bnd_idx, Y);
	}
	else {                // Operate on any dimension:
		//CoreDimNFactorNOrder2(X, Step, nX, nDX, A, B, C, Y);
		CoreDimNSpaceNOrder2(X, Step, nX, nDX, Space, Bnd_idx, Y);
	}

	mxFree(Method);

	return;
}

// =============================================================================
void WrapSpaceN(double *X, const mwSize Step, const mwSize nX,
				const mwSize nDX, double *Space, int8_T *Bnd_idx, double *Y)
{
	// Call different methods depending of the dimensions ofthe input.
	// X has more than one vector. Therefore it is cheaper to calculate the
	// spacing factors once only.

	double *A, *B, *C, *D, *E;
	int8_T *Method;

	// Precalculate spacing factors & methods:
	A = (double *)mxMalloc(nDX * sizeof(double));
	B = (double *)mxMalloc(nDX * sizeof(double));
	C = (double *)mxMalloc(nDX * sizeof(double));
	D = (double *)mxMalloc(nDX * sizeof(double));
	E = (double *)mxMalloc(nDX * sizeof(double));
	Method = (int8_T *)mxMalloc(nX *sizeof(int8_T));

	if (A == NULL || B == NULL || C == NULL) {
		mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
			ERR_HEAD "No memory for Factor vectors.");
	}

	GetFactor(Space, nDX, nX, A, B, C, D, E, Bnd_idx, Method);
	GetMethod(nDX, nX, Bnd_idx);

	if (Step == 1) {        // Operate on first dimension:
		//CoreDim1FactorNOrder2(X, nX, nDX, A, B, C, Y);
		CoreDim1SpaceNOrder2(X, nX, nDX, Space, Bnd_idx, Y);
	}
	else {                // Operate on any dimension:
		//CoreDimNFactorNOrder2(X, Step, nX, nDX, A, B, C, Y);
		CoreDimNSpaceNOrder2(X, Step, nX, nDX, Space, Bnd_idx, Y);
	}

	mxFree(A);
	mxFree(B);
	mxFree(C);
	mxFree(D);
	mxFree(E);
	mxFree(Method);

	return;
}

// =============================================================================
void GetFactor(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
{

}

// =============================================================================
void GetFactor2(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
				double *C, double *D, double *E);
{

}

// =============================================================================
void GetFactor3(double *Space, mwSize nDX, int8_T *Bnd_idx, double double *A, double *B,
	double *C, double *D, double *E);
{

}

// =============================================================================
void GetFactor4(double *Space, mwSize nDX, int8_T *Bnd_idx, double *A, double *B,
	double *C, double *D, double *E);
{

}

