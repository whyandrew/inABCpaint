// 
//  YOU NEED TO MODIFY THIS FILE FOR YOUR INPAINTING
//  IMPLEMENTATION
//
//  DO NOT MODIFY THIS FILE ANYWHERE EXCEPT WHERE 
//  EXPLICITLY NOTED!!!!
//

// See the file inpainting_eval.h for a detailed explanation of
// the input and output parameters of the routines you must 
// implement

#include "inpainting_eval.h"

// When beginning to write your code, I suggest you work in 
// three steps:
// 
//   - First, implement the routine for computing the
//     confidence term. To help with debugging, you should
//     let compute_D() always return 1, so that patches
//     are selected purely based on their confidence term.
//     (correctness of the confidence computation routine
//     should be relatively easy to check since if it 
//     incorrect patches with fewer filled pixels will end up
//     having higher priorities and will be chosen first)
//
//  - Second, implement the routine that does the patch lookup.
//    The correctness of the lookup routine should also be fairly
//    easy to debug, since if it is incorrect the routine will
//    be choosing patches that look nothing like the patch on
//    the fill front.
//    Together, these two steps will allow you to get somewhat
//    reasonable inpaintings done.
//
//  - Third, implement the data term computation routine. You
//    should also try to do this in steps to help with 
//    debugging: (1) make your compute_C() function return 1
//    always, so that patch priorities are computed entirely
//    according to the data term. (2) make your compute_D() 
//    function return just the magnitude of the 
//    gradient computed by the compute_gradient() routine
//    ---this will help you debug gradient computatoins, since an 
//    incorrect computatoin will cause patches that contain 
//    very low intensity variation to be selected before 
//    patches with lots of intensity variations. (3) once
//    you are satisfied that gradient computatoins are 
//    correct, move on to the normal computation routine, etc.
// 
//  - Only when you are satisfied that the above routines
//    are correct individually should you try to compute the
//    priorities as the product C()*D(). Otherwise, if the
//    patch selections 'don't look right' you won't know
//    what is causing this behavior
//
//  

double compute_D(psi& PSI, 
				 const vil_image_view<vxl_byte> im,
				 const vil_image_view<bool> unfilled, 
				 const vil_image_view<bool> fill_front, 
				 double alpha,
				 vnl_double_2& gradient, 
				 vnl_double_2& front_normal)
{
	// holds the perpendicular to the gradient
	vnl_double_2 grad_normal;    

	// compute the gradient at a filled pixel in the 
	// direction of the front normal
	if (compute_gradient(PSI, im, unfilled, gradient)) {

		grad_normal(0) = -gradient(1);
		grad_normal(1) = gradient(0);

		// now compute the normal of the fill front
		if (compute_normal(PSI, fill_front, front_normal)) {
			double dotp;
            
			//dotp = fabs(dot_product(grad_normal, front_normal))/alpha;
            dotp = grad_normal.magnitude();
			return dotp;
		}
		if (alpha > 0) 
			return 1/alpha;
		else {
			return 0;
		}
	} else {
		// if we cannot compute a normal, the fill boundary consists
		// of exactly one pixel; the data term in this case is meaningless
		// so we just return a default value
		return 0;
	}
}

///////////////////////////////////////////////////////////
//     DO NOT CHANGE ANYTHING ABOVE THIS LINE            //
///////////////////////////////////////////////////////////


//
// YOU NEED TO IMPLEMENT THE ROUTINES BELOW
//



double compute_C(psi& PSI, const vil_image_view<double>& C, 
				 const vil_image_view<bool>& unfilled)
{
	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING ABOVE THIS LINE            //
	///////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////

    // Sum confidence of all pixels in the patch and divide by pixel#

    double sum = 0.0;
    int i, j;
    PSI.begin();
    unsigned counter = 0;

    do
    {
        PSI.image_coord(i, j);
        // UNfilled pixel has C==0 anyway, no need to check unfill
        sum += C(i, j);

#ifdef DEBUG320
        counter ++;
#endif

    } while (PSI.next());

    double totalSize = PSI.sz() * PSI.sz();
    return sum / totalSize;

	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING BELOW THIS LINE            //
	///////////////////////////////////////////////////////////
}


bool compute_normal(psi& PSI,
					vil_image_view<bool> fill_front, 
					vnl_double_2& normal)
{
	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING ABOVE THIS LINE            //
	///////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////
	/*	fit a 2nd-order polynomial to each of the curve's coordinate functions,
		using weighted least squares with a Gaussian weight function */
	
    int winRad = PSI.w();

    if (winRad <= 0)
        return false;

	// get local copies of fill_front data first
    vnl_matrix<int> front_mat;
    vnl_matrix<int> valid_mat;
	// This need to be consistent in computing front-normal and patch-gradient
	// DO NOT transpose matrix from get_pixels, use as is,
	//	since only the scalar dot.product is needed at the end anyway.
	// So as long as the relative direction of the front-normal and 
	//	patch-gradient is correct, it's good. Save some operations
    PSI.get_pixels(fill_front, front_mat, valid_mat);
	//front_mat.transpose(); // flip col and row

    int matSize = front_mat.rows();
    int dir = 0;
    // Start at center pixel of patch, need to find it's derivative
    int rowIndex = winRad;
    int colIndex = winRad;
    // vectors for x, y coordinates in paramater t=index
    vnl_vector<int> XCoord(matSize * matSize);
    vnl_vector<int> YCoord(matSize * matSize);
    // Center pixel must be a fill-front by design
    int vecSize = 1;
    XCoord[0] = colIndex;
    YCoord[0] = rowIndex;

    while ((dir = fillFrontDirection(rowIndex, colIndex, front_mat)) >= 0)
    {

    }
	// Construct weighting function matrix
	vnl_matrix<double> weight_mat(3, 3, 0); // only need 2nd-order so 3x3
	weight_mat(0, 0) = exp(-1);
	weight_mat(1, 1) = 1.0;
	weight_mat(2, 2) = weight_mat(0, 0);
	
    

    return true;

	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING BELOW THIS LINE            //
	///////////////////////////////////////////////////////////
}

// return the gradient with the strongest magnitude inside the
// patch of radius w or return false if no gradients can be computed
bool compute_gradient(psi& PSI,
					  const vil_image_view<vxl_byte>& inpainted_grayscale, 
					  const vil_image_view<bool>& unfilled, 
					  vnl_double_2& grad)
{
	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING ABOVE THIS LINE            //
	///////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////
	//              PLACE YOUR CODE HERE                     //
	///////////////////////////////////////////////////////////

	// Use a sliding 3x3 pixel window
	// get max/min pixel to calculate gradient direction & magnitude
	
	if (PSI.w() == 0) return false;

	int winRad = 1;
	int winSize = 2 * winRad + 1; // 2*1 + 1 = 3
	int patSize = PSI.sz();
    int lowestCoord[2];
    int highestCoord[2];
	int highestDiff = -1;

	// Make copies of values
	vnl_matrix<int> valid_mat;
	vnl_matrix<int> unfilled_mat;
	vnl_matrix<int> grayscale_mat;

	PSI.get_pixels(inpainted_grayscale, grayscale_mat, valid_mat);
	PSI.get_pixels(unfilled, unfilled_mat, valid_mat);
 
	grayscale_mat.transpose();
	unfilled_mat.transpose();

	// loop through patch with sliding window
	for (int x = 0; x < patSize - 2 * winRad; x++)
	{
		for (int y = 0; y < patSize - 2 * winRad; y++)
		{
			// Sliding window
            int lowCoord[2];
            int highCoord[2];
            int lowValue = 300;
	        int highValue = -1;
            int pixelValue;
			for (int winX = 0; winX < winSize; winX++)
			{
				for (int winY = 0; winY < winSize; winY++)
				{
                    if (!(bool)unfilled_mat(y+winY, x+winX))
                    {
                        pixelValue = grayscale_mat(y+winY, x+winX);
                        if (pixelValue < lowValue)
                        {
                            lowValue = pixelValue;
                            lowCoord[0] = x+winX;
                            lowCoord[1] = y+winY;
                        }
                        // >= to detect single-pixel only patch vs uniform patch
                        if (pixelValue >= highValue) 
                        {
                            highValue = pixelValue;
                            highCoord[0] = x+winX;
                            highCoord[1] = y+winY;
                        }
                    } // finish 1 pixel of a window
				}
			} // finish one sliding window

            if ((highValue - lowValue) > highestDiff)
            {
                highestDiff = (highValue - lowValue);
                lowestCoord[0] = lowCoord[0];
                lowestCoord[1] = lowCoord[1];
                highestCoord[0] = highCoord[0];
                highestCoord[1] = highCoord[1];
            }
		}
	}// finish whole patch
	
    //Now compute dir and magnitude
    if (highestCoord[0] == lowestCoord[0] &&
        highestCoord[1] == lowestCoord[1]) 
    {  // only 1 single unfilled pixel in patch
        return false;
    }
    else if (highestDiff == 0) 
    { // uniform patch, return arbitary small non-0 value
        grad(0) = 0.01;
        grad(1) = 0;
    }
    else
    {
		int deltaX = highestCoord[0] - lowestCoord[0];
		int deltaY = highestCoord[1] - lowestCoord[1];

        // angle of directrion
        double theta = atan2(deltaY, deltaX);
		// use distance between highest & lowest to inverse-scale diff
		highestDiff *= ( 1 / sqrt(deltaX * deltaX + deltaY * deltaY));

        // set grad to be horizontal with desired mag then rotate
        grad(0) = highestDiff * cos(theta);
        grad(1) = highestDiff * sin(theta);
    }
	return true;

	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING BELOW THIS LINE            //
	///////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////
//     PLACE ANY ADDITIONAL CODE (IF ANY) HERE           //
///////////////////////////////////////////////////////////

/*
    Find connecting pixel of fill_front. and store x,y coord to vect
    Return number of fill_front pixels, i.e. usable size of vector
    Direction:
        0  1  2
        7  x  3
        6  5  4
	Vector index 0 is the center pixel (i.e. t = 0)
	Look for connecting pixel from alternating sides of center pixel
		so it wont' biased all to 1 side, especially in case where 
		the fill_front is closed within the patch.
	Assume t<0 for pixels to left of center pixel, and vice versa.
	Odd index for pixels of t < 0, and even index for pixels of t > 0
	If index is not used (unbalanced # of pixels on +/- sides of t), 
		then set coord to -1.
*/
static inline int fillFrontDirection(
    int row, 
    int col, 
    vnl_vector<int>& rowVect,
    vnl_vector<int>& colVect,
    vnl_matrix<int> fill_mat)
{
    int matSize = fill_mat.rows();
    int vecSize;
    bool haveMorePixel = true;
	bool haveLeftPixel = true;
	bool haveRightPixel = true;
	bool isOdd = true;
	int nextRow;
	int nextCol;

    //sanity check
    if (row >= matSize || row >= matSize || !fill_mat(row, col)) 
        return -1;

    // Starting pixel must be fill_front
    rowVect[0] = row;
    colVect[0] = col;
    vecSize = 1;
	// Remove pixel as it parses so coords won't overlap
	fill_mat(row, col) = 0;

	// Set up the first pair of connecting pixels first to establish left/right
	if (!(haveLeftPixel = checkLeftPixel(
		row, col, rowVect+vecSize, colVect+vecSize, fill_mat)))
	{
		if (!(haveLeftPixel = checkTopBotPixel( row, col, 
			rowVect+vecSize, colVect+vecSize, fill_mat)))
		{
			if (!(haveLeftPixel = checkRightPixel(
				row, col, rowVect+vecSize, colVect+vecSize, fill_mat)))
			{ // single pixel fill_front
				return vecSize; // == 1
			}
		}
	}


	vecSize++;
	haveRightPixel = checkRightPixel(row, col, 
		rowVect+vecSize, colVect+vecSize, fill_mat);




    while (haveMorePixel)
    {
		if (vecSize % 2) // odd index, look for left adjacent pixel
		{
			haveLeftPixel = checkLeftPixel(row, col, 
				rowVect+vecSize, colVect+vecSize, fill_mat);
			if (!haveLeftPixel)
			{
		}
		else // even index, look for right side
		{
			haveRightPixel = checkRightPixel(row, col,
				rowVect+vecSize, colVect+vecSize, fill_mat);
		}

		vecSize++;



        // Check for direct adjacent pixel first
        if (row-1 >= 0 && fill_mat(row-1, col)) // up
        {
            row --;
        }
        else if (row+1 < matSize && fill_mat(row+1, col)) // down
        {
            row ++;
        }
        else if (col+1 < matSize && fill_mat(row, col+1)) //right
        {
            col++;
        }
        else if (col-1 >= 0 && fill_mat(row, col-1)) // left
        {
            col--;
        }
        // Check for diagonal
        else if (row-1 >=0 && col-1 >= 0 && fill_mat(row-1, col-1)) // top-left
        {
            row--;
            col--; 
        }
        else if (row-1 >=0 && col+1 >= 0 && fill_mat(row-1, col+1)) // top-right
        {
            row--;
            col++;
        }
        else if (row+1 >=0 && col-1 >= 0 && fill_mat(row+1, col-1)) // bottom-left
        {
            row++;
            col--;
        }
        else if (row+1 >=0 && col+1 >= 0 && fill_mat(row+1, col+1)) // bottom-right
        {
            row++;
            col++;
        }
        else
        {
            haveMorePixel = false;
        }

        if (haveMorePixel)
        {
        }
    }
    
    return vecSize;
}

// Check if any left adjacent pixel, including upper & lower left corners
// Set return pixel coord to -1 if none is found
static inline bool checkLeftPixel(
	const int row,
	const int col,
	int& retRow,
	int& retCol,
	vnl_matrix<int>& front_mat)
{
	return true;
}

// Check if any right adjacent pixel, including upper & lower right
static inline bool checkRightPixel(
	const int row,
	const int col, 
	int& retRow,
	int& retCol,
	vnl_matrix<int>& front_mat)
{
	return true;
}

// check for direct adjacent pixel on top and bottom
static inline bool checkTopBotPixel(
	const int row,
	const int col,
	int& retRow,
	int& retCol,
	vnl_matrix<int>& front_mat)
{
	return true;
}

/////////////////////////////////////////////////////////

