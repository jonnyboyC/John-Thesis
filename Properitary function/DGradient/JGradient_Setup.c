// JGradient_Setup.c
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
#define ERR_HEAD "*** JGradient_Setup[mex]: "
#define ERR_ID   "JChabot:JGradient_Setup:"

void GetMethodDimCol(mwSize *Bnd_idx, mwSize *Method, const mwSize nX,
	const mwSize nDX);
void GetMethodDimRow(mwSize *Bnd_idx, mwSize *Method, const mwSize nX,
	const mwSize *dimX, const mwSize Step);
mwSize FirstNonSingeltonDim(const mwSize Xndim, const mwSize *Xdim);
mwSize GetStep(const mwSize *Xdim, const mwSize N);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *A, *B, *C, *D, *E, *Method_Double, Nd;
	mwSize *Bnd_idx, *Method, nX, ndimX, N, nDX, Step;
	const mwSize *dimX;

	// Check number and type of inputs and outputs: ------------------------------
	if (nrhs == 0 || nrhs > 2) {
		mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
			ERR_HEAD "1 or 2 inputs required.");
	}

	if (nlhs != 6) {
		mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
			ERR_HEAD "6 outputs required.");
	}

	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0])) {
		mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput1",
			ERR_HEAD "Input must be a full real double array.");
	}

	// Pointers and dimension to input array: ------------------------------------
	Bnd_idx = (mwSize *)mxGetData(prhs[0]);
	nX = mxGetNumberOfElements(prhs[0]);
	ndimX = mxGetNumberOfDimensions(prhs[0]);
	dimX = mxGetDimensions(prhs[0]);

	// Return fast on empty input matrix:
	if (nX == 0) {
		plhs[0] = COPY_ARRAY(prhs[0]);
		plhs[1] = COPY_ARRAY(prhs[0]);
		plhs[2] = COPY_ARRAY(prhs[0]);
		plhs[3] = COPY_ARRAY(prhs[0]);
		plhs[4] = COPY_ARRAY(prhs[0]);
		plhs[5] = COPY_ARRAY(prhs[0]);
		return;
	}

	// Check 2nd input, determine dimension to operate on: ----------------------------------------
	if (nrhs < 1) {
		N = FirstNonSingeltonDim(ndimX, dimX);  // Zero based
		Step = 1;
		nDX = dimX[N];

	}
	else if (mxIsNumeric(prhs[1]))  {  // 3rd input used:
		switch (mxGetNumberOfElements(prhs[1])) {
		case 0:  // Use 1st non-singelton dim if 3rd input is []:
			N = FirstNonSingeltonDim(ndimX, dimX);
			Step = 1;
			nDX = dimX[N];
			break;

		case 1:  // Numerical scalar:
			Nd = mxGetScalar(prhs[1]);
			N = (mwSize)Nd - 1;
			if (Nd < 1.0 || Nd != floor(Nd)) {
				mexErrMsgIdAndTxt(ERR_ID   "BadValueInput3",
					ERR_HEAD "Dimension must be a positive integer scalar.");
			}

			if (N < ndimX) {
				Step = GetStep(dimX, N);
				nDX = dimX[N];
			}
			else {
				// Treat imaginated trailing dimensions as singelton, as usual in
				// Matlab:
				Step = nX;
				nDX = 1;
			}
			break;

		default:
			mexErrMsgIdAndTxt(ERR_ID   "BadSizeInput3",
				ERR_HEAD "3rd input [Dim] must be scalar index.");
		}

	}
	else {  // 2nd input is not numeric:
		mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput3",
			ERR_HEAD "3rd input must be scalar index.");
	}

	// Allocate output matrices: --------------------------------------------------
	plhs[0] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
	plhs[3] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
	plhs[4] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);
	plhs[5] = mxCreateNumericArray(ndimX, dimX, mxDOUBLE_CLASS, mxREAL);

	A = mxGetPr(plhs[0]);
	B = mxGetPr(plhs[1]);
	C = mxGetPr(plhs[2]);
	D = mxGetPr(plhs[3]);
	E = mxGetPr(plhs[4]);
	Method_Double = mxGetPr(plhs[5]);

	Method = (mwSize *)mxMalloc(sizeof(mwSize)*nX);


	// Determine methods and factors: -----------------------------------------------
	if (Step == 1) {
		GetMethodDimCol(Bnd_idx, Method, nX, nDX);
	}
	else {
		GetMethodDimRow(Bnd_idx, Method, nX, dimX , Step);
	}

	for (mwSize i = 0; i < nX; ++i) {
		*Method_Double++ = *Method++;
	}

	mxFree(Method);
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
mwSize GetStep(const mwSize *dimX, const mwSize N)
{
	// Get step size between elements of a subvector in the N'th dimension.
	// This is the product of the leading dimensions.

	const mwSize *XdimEnd, *XdimP;
	mwSize       Step;

	Step = 1;
	XdimEnd = dimX + N;
	for (XdimP = dimX; XdimP < XdimEnd; Step *= *XdimP++); // empty loop

	return (Step);
}

void GetMethodDimCol(mwSize *Bnd_idx, mwSize *Method, const mwSize nX,
	const mwSize nDX)
{
	mwSize *Bnd_idxF = Bnd_idx + nX;
	mwSize *Bnd_idxC = Bnd_idx + nDX;
	int gap = 0;

	while (Bnd_idx < Bnd_idxF) {

		// Check if Points has moved off column
		if (Bnd_idx == Bnd_idxC) {
			Bnd_idxC += nX;
		}

		// If in boundary set method to 0 and continue
		if (*Bnd_idx == -1) {
			*Method++ = 0;
			Bnd_idx++;
			continue;
		}

		// Count consecutive spaces in open flow
		while (*Bnd_idx++ != -1 && Bnd_idx < Bnd_idxC){
			gap++;
		}

		// Assign Method based on gap size
		switch (gap) {
		case 1:
			*Method++ = 1;
			break;
		case 2:
			*Method++ = 2;
			*Method++ = 3;
			break;
		case 3:
			*Method++ = 4;
			*Method++ = 5;
			*Method++ = 6;
			break;
		case 4:
			*Method++ = 7;
			*Method++ = 8;
			*Method++ = 9;
			*Method++ = 10;
			break;
		default:
			*Method++ = 11;
			*Method++ = 12;
			for (int i = 0; i < gap - 4; ++i) {
				*Method++ = 13;
			}
			*Method++ = 14;
			*Method++ = 15;
		}
		gap = 0;
	}
}

void GetMethodDimRow(mwSize *Bnd_idx, mwSize *Method, const mwSize nX,
	const mwSize *Xdim, const mwSize Step)
{
	mwSize rows = Xdim[0], cols = Xdim[1];
	mwSize row = 0, col = 0, col_method = 0;
	mwSize *Bnd_idxF = Bnd_idx + nX;
	int gap = 0;

	while ((Bnd_idx + row + Step*col) < Bnd_idxF) {

		// Check if Points has moved off column
		if (col == cols) {
			row++;
			col = 0;
		}

		col_method = col;

		// If in boundary set method to 0 and continue
		if (*(Bnd_idx+row + Step*col) == -1) {
			*(Method + row + Step*col_method++) = 0;
			col++;
			continue;
		}

		// Count consecutive spaces in open flow
		while (*(Bnd_idx + row + Step*col++) != -1 && col < cols){
			gap++;
		}

		// Assign Method based on gap size
		switch (gap) {
		case 1:
			*(Method + row + Step*col_method++) = 1;
			break;
		case 2:
			*(Method + row + Step*col_method++) = 2;
			*(Method + row + Step*col_method++) = 3;
			break;
		case 3:
			*(Method + row + Step*col_method++) = 4;
			*(Method + row + Step*col_method++) = 5;
			*(Method + row + Step*col_method++) = 6;
			break;
		case 4:
			*(Method + row + Step*col_method++) = 7;
			*(Method + row + Step*col_method++) = 8;
			*(Method + row + Step*col_method++) = 9;
			*(Method + row + Step*col_method++) = 10;
			break;
		default:
			*(Method + row + Step*col_method++) = 11;
			*(Method + row + Step*col_method++) = 12;
			for (int i = 0; i < gap - 4; ++i) {
				*(Method + row + Step*col_method++) = 13;
			}
			*(Method + row + Step*col_method++) = 14;
			*(Method + row + Step*col_method++) = 15;
		}
		gap = 0;
	}
}