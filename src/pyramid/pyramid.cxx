// DO NOT MODIFY THIS FILE EXCEPT WHERE EXPLICITLY NOTED!!!

#include "pyramid.h"

////////////////////////////////////////////////////////////////
//          The pyramid class constructor routines            //
////////////////////////////////////////////////////////////////


pyramid::pyramid(const vil_image_view<vxl_byte>& im)
{
	// set the default value for the alpha parameter to 0.4
	a_ = 0.4;

	build(im);
}

pyramid::pyramid(const vil_image_view<vxl_byte>& im, double a)
{
	a_ = a;

	build(im);
}

//
// This is the main routine for building the Gauss and the
// Laplacian pyramid
//
// It takes as input an image (of arbitrary dimensions MxP)
// and builds a pyramid whose 0 level is of size (2^N+1)x(2^N+1)
// and N is the smallest power of 2 such that 2^N+1 > max(M,P)
//
void pyramid::build(const vil_image_view<vxl_byte>& im)
{
	int ni, nj;
	int l;
	int n;
	int i,j, plane;

	//
	// initialize the smoothing kernel
    //

	// the kernel always has length 5 pixels
	w_hat_ = new double[5];
	// shift pointer by two so that kernel indices are in
	// the range [-2,2]
	w_hat_ += 2;
	// set the values of the kernel
	w_hat_[2] = w_hat_[-2] = 0.25 - a_/2;
	w_hat_[1] = w_hat_[-1] = 0.25;
	w_hat_[0] = a_;

	// 
	// Allocate the arrays holding the Gauss and Laplacian pyramid images
	// 

	// compute the number of levels, N_
	ni = (int)ceil(log((double)im.ni()-1)/log(2));
	nj = (int)ceil(log((double)im.nj()-1)/log(2));
	N_ = vcl_max(ni,nj);

	// 
    // a pyramid is represented as an array of N_+1 images, 0,...,N_
	// where the n-th image is of dimension (2^n+1)x(2^n+1)
	// 

	// allocate a temporary array holding the Gaussian pyramid images
	vil_image_view<vxl_byte>* g_temp = new vil_image_view<vxl_byte>[N_+1];

	// allocate the array holding the Laplacian pyramid images
	L_ = new vil_image_view<int>[N_];

	// 
	// Build the pyramids
	//

	// Step 1: Copy the original image into level 0 of the Gauss pyramid

	// allocate space and initialize the level-0 image with zeros, since 
	// the input image may not occupy all (2^N_+1)x(2^N+1) pixels
	g_temp[0].set_size((int)pow(2,N_)+1, (int)pow(2,N_)+1, im.nplanes());
	g_temp[0].fill(0);

	// copy the image so that pixel im(0,0) goes to pixel g_temp[0](0,0)
	vil_copy_to_window(im, g_temp[0], 0, 0);


	// Step 2: Compute levels 1,...,N_ of the Gauss pyramid

	for (l=1; l<=N_; l++) 
		// each level is a reduced version of the image immediately below it
		reduce(g_temp[l-1], w_hat_, g_temp[l]);

	// Step 3: Compute levels 0,...,N_-1 of the Laplacian pyramid

	for (l=0; l<N_; l++) {
		// allocate space for the l-th level of the Laplacian pyramid
		L_[l] = vil_image_view<int>(g_temp[l].ni(), g_temp[l].nj(), 
			                              g_temp[l].nplanes());
		vil_image_view<vxl_byte> gtmp;

		// each level is computed by the difference
		//   L_[l] = g_[l] - expand(g_[l+1])
		// i.e., it represents all the "details" in g_[l] that are
		// "lost" when the image is reduced from g_[l] to g_[l+1]

		expand(g_temp[l+1], w_hat_, gtmp);
		vil_math_image_difference(g_temp[l], gtmp, L_[l]);
	}

	// Store just the N_ level image of the Gauss pyramid
	g_N_ = g_temp[N_];
}

// Create a pyramid data structure from a sequence of
// Laplacian pyramid levels and the N-th level Gauss image
pyramid::pyramid(const vil_image_view<int>* L, 
				 vil_image_view<vxl_byte> g_N, 
				 int N, double a)
{
	int l;

	N_ = N;
	a_ = a;

	// the kernel always has length 5 pixels
	w_hat_ = new double[5];
	// shift pointer by two so that kernel indices are in
	// the range [-2,2]
	w_hat_ += 2;
	// set the values of the kernel
	w_hat_[2] = w_hat_[-2] = 0.25 - a_/2;
	w_hat_[1] = w_hat_[-1] = 0.25;
	w_hat_[0] = a_;

	L_ = new vil_image_view<int>[N_];
	for (l=0; l<N; l++)
		vil_copy_deep(L[l], L_[l]);

	vil_copy_deep(g_N, g_N_);
}

//
// Pyramid accessor routines
// 


// Return the total number of levels in the pyramid
int pyramid::N() const
{
	return N_;
}

// Return the alpha value of the w_hat kernel
double pyramid::a() const
{
	return a_;
}

// Copy level l of the Laplacian pyramid into the supplied
// image; the routine returns FALSE if (l is not in the 
// range [0,...,N_-1]
bool pyramid::L(int l, vil_image_view<int>& L_l) const
{
	if ((l >=0) && (l < N_)) {
		vil_copy_deep(L_[l], L_l);
		return true;
	} else
		return false;
}

// Copy level l1 of the Laplacian pyramid into the supplied
// image, and expend it to level l2; the routine returns 
// FALSE if (l1 is not in the  range [0,...,N_-1]) or
// if (l1 < l2)
bool pyramid::L(int l1, int l2, vil_image_view<int>& L_l) const
{
	if ((l1 >=0) && (l1 < N_) && (l1 >= l2)) {
		vil_image_view<int> temp1, temp2;

		// get level l1 of the pyramid
		L(l1, temp1);

		// expand the level to level l2
		for (int l=l1; l>l2; l--) {
			expand(temp1, w_hat_, temp2);
			temp1 = temp2;
		}
		vil_copy_deep(temp1, L_l);

		return true;
	} else
		return false;
}

// Copy all levels of the Laplacian pyramid into an array of
// images and return a pointer to that array

vil_image_view<int>* pyramid::L() const
{
	// allocate the array holding the images
	vil_image_view<int>* L_temp = new vil_image_view<int>[N_];

	// copy the corresponding levels into the array
	for (int l=0; l<N_; l++) 
		vil_copy_deep(L_[l], L_temp[l]);
	
	return L_temp;
}


////////////////////////////////////////////////////////////////
//        Routines for computing the Gauss pyramid            //
//                                                            //
//   The gauss pyramid is not stored explicitly, so we        //
//   need to reconstruct it from the Laplacian images upon    //
//   request                                                  //
//                                                            //
////////////////////////////////////////////////////////////////

// a helper function that adds two images and also performs intensity bound
// checking to avoid overflows when assigning numbers to unsigned byte
// images
static vil_image_view<vxl_byte> add_g_and_L(vil_image_view<vxl_byte>& g, 
						                    vil_image_view<int>& L)
{
	vil_image_view<vxl_byte> result(g.ni(), g.nj(), g.nplanes());
	
	for (int p=0; p<g.nplanes(); p++)
		for (int i=0; i<g.ni(); i++)
			for (int j=0; j<g.nj(); j++) {
				int value = g(i,j,p) + L(i,j,p);
				result(i,j,p) = vcl_min(vcl_max(value, 0), 255);
			}

	return result;
}


// Return an array that holds the entire Gauss pyramid
// The routine reconstructs the Gauss pyramid from the Laplacian
// pyramid images

vil_image_view<vxl_byte>* pyramid::g() const
{
	// allocate the array holding the images
	vil_image_view<vxl_byte>* g_temp = new vil_image_view<vxl_byte>[N_+1];

	// the N_-th level of the Gauss pyramid is only level stored in the
	// pyramid data structure
	g_temp[N_] = g_N_;

	// the remaining levels have to be computed
	for (int l=N_; l>0; l--) {
		vil_image_view<vxl_byte> gtmp;
		// compute the l-th level of the Gauss pyramid using the
		// formula g_(l-1) = expand(g_l) + L_[l-1]
		expand(g_temp[l], w_hat_, gtmp);
		g_temp[l-1] = add_g_and_L(gtmp, L_[l-1]);
	}

	return g_temp;
}


// Compute the l-th level of the Gauss pyramid from the 
// Laplacian pyramid images; the routine
// returns FALSE if l is not in the range [0,...,N_]
bool pyramid::g(int l, vil_image_view<vxl_byte>& g_l) const
{
	return g(l, l, g_l);
}

// Compute the l1-th level of the Gauss pyramid from the 
// Laplacian pyramid images. The routine then expands 
// the l1-th level image repeatedly so that its size 
// becomes (2^l2 + 1)x(2^l2 + 1).
//
// The routine returns FALSE if l1,l2 are outside the 
// range [0,...,N_] and/or if l2 > l1.
// 
// If (l1=l2), the routine simply computes the l1-level image of the
// Gauss pyramid


bool pyramid::g(int l1, int l2, vil_image_view<vxl_byte>& g_l) const
{
	int l;
	vil_image_view<vxl_byte> temp;

    if ((l1 >= 0) && (l1 <= N_) &&
		(l2 >= 0) && (l2 <= N_) &&
		(l2 <= l1)) {

        // the N_-th level of the Gauss pyramid is only level stored in the
        // pyramid data structure
		vil_image_view<vxl_byte> g_l1;
		g_l1 = g_N_;

		// the remaining levels have to be computed
		for (l=N_; l>l1; l--) {
			vil_image_view<vxl_byte> gtmp;
			// compute the l-th level of the Gauss pyramid using the
			// formula g_(l-1) = expand(g_l) + L_[l-1]
			expand(g_l1, w_hat_, gtmp);
			temp = add_g_and_L(gtmp, L_[l-1]);
			g_l1 = temp;

		}

		// create gl_2 by expanding until level l2
		for (l=l1, g_l=g_l1; l>l2; l--) {
			// implement the formula g_l = expand(g_l)			
			expand(g_l, w_hat_, temp);
			g_l = temp;
		}
		return true;
	} else
		return false;
}

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING ABOVE THIS LINE
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//               The expand/reduce routines                   //
////////////////////////////////////////////////////////////////

//
// The REDUCE() routine
// 

// Given image of dimensions (2^N+1)x(2^N+1), return an image
// of dimensions (2^(N-1)+1)x(2^(N-1)+1) that is a smoothed and
// subsampled version of the original
// 
// Parameters:
//    im:     the input image
//    w_hat:  the 5-pixel smoothing kernel (w_hat[-2],...,w_hat[2])
//
// Boundary treatment: When the kernel extends beyond the image 
//                     boundary, you should compute the convolution
//                     (1) by taking into account only the kernel
//                     elements that overlap with the image
//                     and (2) by dividing the result by the sum of 
//                     these kernel elements

void pyramid::reduce(const vil_image_view<vxl_byte> im, 
			         const double* w_hat,
			         vil_image_view<vxl_byte>& im_red)
{
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////

    // Smooth every even pixel on 1-dimension from input image "im",
    // and save to reduced imaged im_red
    // Need to loop thru each plane of the image

    // Allocate mem for im_red
    int imSize = (im.ni() / 2) + 1;
    im_red.set_size(imSize, imSize, im.nplanes());
    // A temp intermediate image after 1-dimension horizontally reduced
    vil_image_view<vxl_byte> tmpIm = vil_image_view<vxl_byte>(imSize, im.nj(), im.nplanes());

    int maxI = im.ni();
    int maxJ = im.nj();
    const int maxIminus1 = maxI - 1;
    const int maxJminus1 = maxJ - 1;
    double value = 0;
    double edgeWeight = w_hat[0] + w_hat[1] + w_hat[2]; //assume symetrical kernel
    double fullWeight = 0;
    for (int i = -2; i < 3; i++)
        fullWeight += w_hat[i];

    for (int p = 0; p < im.nplanes(); p++)
    {
        //Reduce horizontally
        for (int j = 0; j < maxJ; j++)
        {
            for (int i = 0; i < maxI; i += 2)
            {
                value = 0;
                //Check for edge cases
                if (i == 0)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        value += (w_hat[k] * im(i+k, j, p));
                    }
                    //Scale it back to max of 1
                    value = value / edgeWeight;
                    tmpIm(0, j, p) = (vxl_byte)floor(value+0.49999);
#ifdef DEBUG320
                    double dbgCheckPt = tmpIm(0, j, p);
                    dbgCheckPt = 0;
#endif
                }
                else if (i == maxIminus1)
                {
                    for (int k = -2; k < 1; k++)
                    {
                        value += (w_hat[k] * im(i+k, j, p));
                    }
                    //Scale it back to max of 1
                    value = value / edgeWeight;
                    tmpIm(i/2, j, p) = (vxl_byte)floor(value+0.49999);
#ifdef DEBUG320
                    double dbgCheckPt = tmpIm(0, j, p);
                    dbgCheckPt = 0;
#endif
                }
                else
                { // non-edge cases
                    for (int k = -2; k < 3; k++)
                    {
                        value += (w_hat[k] * im(i+k, j, p));
                    }
                    // Kernel is already weighed = 1
                    tmpIm(i/2, j, p) = (vxl_byte)floor((value/fullWeight)+0.49999);
#ifdef DEBUG320
                    double dbgCheckPt = tmpIm(i/2, j, p);
                    dbgCheckPt = 0;
#endif
                }
            }
        }

        //Reduce vertically
        for (int j = 0; j < maxJ; j += 2)
        {
            for (int i = 0; i < imSize; i++)
            {
                value = 0;
                //Check for edge cases
                if (j == 0)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        value += (w_hat[k] * tmpIm(i, j+k, p));
                    }
                    //Scale it back to max of 1
                    value = value / edgeWeight;
                    im_red(i, j/2, p) = (vxl_byte)floor(value+0.49999);
#ifdef DEBUG320
                    double dbgCheckPt = im_red(i, j/2, p);
                    dbgCheckPt = 0;
#endif
                }
                else if (j == maxJminus1)
                {
                    for (int k = -2; k < 1; k++)
                    {
                        value += (w_hat[k] * tmpIm(i, j+k, p));
                    }
                    //Scale it back to max of 1
                    value = value / edgeWeight;
                    im_red(i, j/2, p) = (vxl_byte)floor(value+0.49999);
                }
                else
                { // non-edge cases
                    for (int k = -2; k < 3; k++)
                    {
                        value += (w_hat[k] * tmpIm(i, j+k, p));
                    }
                    // Kernel is already weighed = 1
                    im_red(i, j/2, p) = (vxl_byte)floor((value/fullWeight)+0.49999);
#ifdef DEBUG320
                    double dbgCheckPt = im_red(i, j/2, p);
                    dbgCheckPt = 0;
#endif
                }
            }
        }
    }

}

//
// The EXPAND() routine
// 

// Given image of dimensions (2^(N-1)+1)x(2^(N-1)+1), return an image
// of dimensions (2^N+1)x(2^N+1) that is an interpolated version
// of the original using the w_hat as the interpolating kernel
// 
// Parameters:
//    im:     the input image
//    w_hat:  the 5-pixel smoothing kernel (w_hat[-2],...,w_hat[2])
//
//
// Boundary treatment: When the kernel extends beyond the image 
//                     boundary, you should compute the convolution
//                     (1) by taking into account only the kernel
//                     elements that overlap with the image
//                     and (2) by dividing the result by the sum of 
//                     these kernel elements

void pyramid::expand(const vil_image_view<vxl_byte> im, 
		             const double* w_hat, 
		             vil_image_view<vxl_byte>& im_exp)
{
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////

    // Just convert type and call the byte version of "expand" 
    vil_image_view<int> tmpIm = vil_image_view<int>(im.ni(), im.nj(), im.nplanes());
    vil_convert_cast(im, tmpIm);

    vil_image_view<int> tmpIm_exp;
    expand(tmpIm, w_hat, tmpIm_exp);
    vil_convert_cast(tmpIm_exp, im_exp);
    //int_to_ubyte(tmpIm_exp, im_exp); 
    // ^^^^ that function's description was NOT CLEAR, wasted me half a day damnit
    // Would be good to mention all values are increased by 127!!! before casting to ubyte!
 }

void pyramid::expand(const vil_image_view<int> im, 
		             const double* w_hat, 
		             vil_image_view<int>& im_exp)
{
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////

    // Make copies of kernel for even and odd pixel
    // So no need to reweight value of every pixel
    // Assume kernel is symmetrical, so kernel for odd pixel is just {0.5, 0.5}
    double *evenKern;
    // "Even" is the one of size 3, since we start index=0
    evenKern = new double[3];
    double weight = 0;
    for (int i = -2, j = 0; i < 3; i+=2, j++)
    {
        evenKern[j] = w_hat[i];
        weight += w_hat[i];
    }
    for (int i = 0; i< 3; i++)
    {
        evenKern[i] /= weight;
    }
    // Change evenKern index from -1 to 1
    evenKern = evenKern + 1;

    //Allocate mem to im_exp
    int imSize = (im.ni() - 1) * 2 + 1;
    im_exp.set_size(imSize, imSize, im.nplanes());

    double value = 0;
    bool isEven = true;
    int imSizeMinus1 = imSize - 1;

    for (int p=0; p < im.nplanes(); p++)
    {
        //Expand horizontally
        // So alternating ROWS will be fully filled in the expanded image
        for (int j = 0; j < imSize ; j += 2)
        {
            //Deal with the edge cases first
            // first and last pixels, both are even
            weight = evenKern[0] + evenKern[1];
            value = (evenKern[0] * im(0, j/2, p)) + (evenKern[1] * im(1, j/2, p));
            value = value / weight;
            im_exp(0, j, p) = (int)floor(value+0.49999);

            value = (evenKern[-1] * im(im.ni() - 2, j/2, p)) +
                (evenKern[0] * im(im.ni() - 1, j/2, p));
            value = value / weight;
            im_exp(imSizeMinus1, j, p) = (int)floor(value+0.49999);

            isEven = false;
            // Below are the non-edge cases
            for (int i = 1; i < imSizeMinus1; i++)
            {
                value = 0.0;
                // The current i,j are index of the expanded image
                // Need index of the original image when getting value        
                if (isEven)
                {
                    //Intepolate from 3 pixels
                    for (int k = -1; k < 2; k++)
                    {
                        value += (im((i/2)+k, j/2, p) * evenKern[k]);
                    }
                }
                else //odd index
                {
                    //Intepolate from 2 pixels only
                    value = (im(i/2, j/2, p) + im((i/2)+1, j/2, p)) / 2;
                }

                isEven = !isEven;
                // Kernel is already re-weighted, should not need to normalize here
                im_exp(i, j, p) = (int)floor(value+0.49999);
#ifdef DEBUG320
                double dbgCheckPt = im_exp(i, j, p);
                dbgCheckPt = 0;
#endif
            }
        }

        //Expand vertically
        // ONLY odd index rows to be filled, each of the pixel interpolate vertically
        for (int j = 1; j < imSize; j += 2)
        {
            for (int i = 0; i < imSize; i++)
            {
                // No edge case, 
                // since first and last rows (both even index) should be filled already 
                im_exp(i, j, p) = (im_exp(i, j-1, p) + im_exp(i, j+1, p)) / 2 ;
#ifdef DEBUG320
                double = dbgCheckPt = im_exp(i, j, p);
                dbgCheckPt = 0;
#endif
            }
        }
    }
}

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING BELOW THIS LINE
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//           Pyramid class utility routines                   //
////////////////////////////////////////////////////////////////

//
// Crop an input image so that it becomes square,
// has dimensions 2^N + 1, and 2^N+1 is both smaller
// and closest to the largest dimension of the image
//

void pyramid::crop_to_power_of_2plus1(
		const vil_image_view<vxl_byte>& im,
		vil_image_view<vxl_byte>& im_crop)
{
	int ni, nj, N;

	// compute the number of levels, N
	ni = (int)floor(log((int)im.ni())/log(2));
	nj = (int)floor(log((int)im.nj())/log(2));
	N = vcl_max(ni,nj);

	im_crop = vil_crop(im, 0, (int)(pow(2,N)+1), 0, (int)(pow(2,N)+1));
}

// 
// Write levels l1 to l2 of the Gaussian pyramid to 
// disk. If expand_to_l2=true, all images are expanded 
// to the size of the level l2 images. 
//
// The routine also outputs a packed version of the pyramid
// into a single image
void pyramid::dump_gauss(int l1, int l2, bool expand_to_l2, const char* basename)
{
	int i;

	for (int l=l2; l<=l1; l++) {
		vil_image_view<vxl_byte> im;
		char fname[256];
		vcl_ostringstream pyr_fname(fname);

		if (expand_to_l2 == true) {
			pyr_fname << basename << ".g_exp." << l << ".jpg" << vcl_ends;
			g(l, l2, im);
		} else {
			pyr_fname << basename << ".g." << l << ".jpg" << vcl_ends;
			g(l, im);
		}

		 
		// to save, we need to access a (char *) representation
		// of the output string stream
		vcl_cerr 
			<< "Saving Gauss pyramid level " << l << " to file "
			<< (pyr_fname.str()).c_str() << "\n";
		vil_save(im, (pyr_fname.str()).c_str());
	}

	// save a packed version of the pyramid
	char fname[256];
	vcl_ostringstream pack_fname(fname);
	vil_image_view<vxl_byte> im;

	pack_fname << basename << ".g.pack.jpg" << vcl_ends;
	vcl_cerr 
		<< "Saving packed Gauss pyramid to file "
		<< (pack_fname.str()).c_str() << "\n";
	pack_gauss(im);
	vil_save(im, (pack_fname.str()).c_str());
}

void pyramid::dump_laplacian(int l1, int l2, bool expand_to_l2, const char* basename)
{
	int i;

	if (l1 >= N_)
		l1 = N_-1;
	if (l1 < 0) 
		l1 = 0;
	if (l2 >= N_)
		l2 = N_;
	if (l2 < 0)
		l2 = 0;

	for (int l=l2; l<=l1; l++) {
		vil_image_view<int> im;
		vil_image_view<vxl_byte> imb;
		char fname[256];
		vcl_ostringstream pyr_fname(fname);

		if (expand_to_l2 == true) {
			pyr_fname << basename << ".L_exp." << l << ".jpg" << vcl_ends;
			L(l, l2, im);
		} else {
			pyr_fname << basename << ".L." << l << ".jpg" << vcl_ends;
			L(l, im);
		}

		 
		// to save, we need to access a (char *) representation
		// of the output string stream
		vcl_cerr 
			<< "Saving Laplacian pyramid level " << l << " to file "
			<< (pyr_fname.str()).c_str() << "\n";
		int_to_ubyte(im, imb);
		vil_save(imb, (pyr_fname.str()).c_str());
	}

	// save a packed version of the pyramid
	char fname[256];
	vcl_ostringstream pack_fname(fname);
	vil_image_view<vxl_byte> im;

	pack_fname << basename << ".L.pack.jpg" << vcl_ends;
	vcl_cerr 
		<< "Saving packed Laplacian pyramid to file "
		<< (pack_fname.str()).c_str() << "\n";
	pack_laplacian(im);
	vil_save(im, (pack_fname.str()).c_str());

}

// 
// Return an image that contains all levels of the Laplacian pyramid
// This routine is useful for pyramid visualization purposes
// 
void pyramid::pack_laplacian(vil_image_view<vxl_byte>& imb) const
{
	// allocate space for the output image; the packed image has
	// 1/2 more columns and the same number of rows as the original image 
	vil_image_view<int> im(2*(L_[0].ni()-1), 
		                   L_[0].nj(), 
						   L_[0].nplanes());
	// Note: a more efficient packing is possible:
	// allocate space for the output image; the packed image has
	// 1/2 more columns and the same number of rows as the original image 
	//vil_image_view<int> im(L_[0].ni()+(L_[0].ni()-1)/2+1, 
	//	                   L_[0].nj(), 
	//					   L_[0].nplanes());

	// fill it with zeros
	im.fill(0);

	// call the recursive packing routine
	pack(L_, N_, 0, 0, im);

	// since the Laplacian images contain signed byte values, we first convert it to an
	// unsigned byte image by mapping the range [-127..128] to [0..255]
	int_to_ubyte(im, imb);
}

// convert an int-type image (used for storing signed Laplacian images)
// to a ubyte image
void pyramid::int_to_ubyte(const vil_image_view<int>& imi, vil_image_view<vxl_byte>& imb)
{
	imb.set_size(imi.ni(), imi.nj(), imi.nplanes());
	for (int p=0; p<imi.nplanes(); p++)
		for (int i=0; i<imi.ni(); i++)
			for (int j=0; j<imi.nj(); j++)
				imb(i,j,p) = vcl_min<int>(imi(i,j,p)+127,255);

}

// 
// Return an image that contains all levels of the Gauss pyramid
// This routine is useful for pyramid visualization purposes. 
//
void pyramid::pack_gauss(vil_image_view<vxl_byte>& im) const
{
	// Since the Gaussian pyramid is not stored explicitly, we first
	// compute & store the Gaussian pyramid in a temporary array
	// of images
	vil_image_view<vxl_byte>* g_temp = g();

	// allocate space for the output image; the packed image has
	// the same number of rows as the original image and twice the 
	// columns
	im.set_size(2*(g_temp[0].ni()-1), g_temp[0].nj(), g_temp[0].nplanes());
	// note: a tigher packing is also possible
	//im.set_size(g_temp[0].ni()+(g_temp[0].ni()-1)/2+1, 
	//	                        g_temp[0].nj(), 
	//                           g_temp[0].nplanes());

	// fill it with zeros
	im.fill(0);
	// call the recursive pyramid-packing routine
	pack(g_temp, N_, 0, 0, im);
}

// Recursive packing routine for unsigned byte images
void pyramid::pack(vil_image_view<vxl_byte>* pyr, int l, int i, int j, 
				   vil_image_view<vxl_byte>& im) const
{
	if (l > 0) {
		vil_copy_to_window(pyr[0], im, i, j);
		if (l > 1) {
			vil_copy_to_window(pyr[1], im, i+pyr[0].ni(), j);
			if (l > 2) {
				vil_copy_to_window(pyr[2], im, i+pyr[0].ni(), j+pyr[1].nj());
				if (l >= 3) 
					pack(pyr+3, l-3, i+pyr[0].ni()+pyr[2].ni(), j+pyr[1].nj(), im);
			}
		}
	}
}

// Recursive packing routine for signed byte images (useful for 
// packing Laplacian images). 
void pyramid::pack(vil_image_view<int>* pyr, int l, int i, int j, 
				   vil_image_view<int>& im) const
{
	if (l > 0) {
		vil_copy_to_window(pyr[0], im, i, j);
		if (l > 1) {
			vil_copy_to_window(pyr[1], im, i+pyr[0].ni(), j);
			if (l > 2) {
				vil_copy_to_window(pyr[2], im, i+pyr[0].ni(), j+pyr[1].nj());
				if (l >= 3) 
					pack(pyr+3, l-3, i+pyr[0].ni()+pyr[2].ni(), j+pyr[1].nj(), im);
			}
		}
	}
}


