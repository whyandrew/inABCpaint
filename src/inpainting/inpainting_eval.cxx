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
        //compute_normal(PSI, fill_front, front_normal);
			double dotp;
            
			dotp = fabs(dot_product(grad_normal, front_normal))/alpha;
            //dotp = grad_normal.magnitude();
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
    count_conf++;
    TIMER_START;

    double sum = 0.0;
    int i, j;
    PSI.begin();

    do
    {
        PSI.image_coord(i, j);
        // UNfilled pixel has C==0 anyway, no need to check unfill
        sum += C(i, j);
    } while (PSI.next());

    double totalSize = PSI.sz() * PSI.sz();

    time_conf += TIMER_ELLAPSED;

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
    count_normal++;
    TIMER_START;

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

#ifdef DEBUG320
    		vcl_cerr 
			<< "------------------------------------------" << vcl_endl
			<< "compute_normal patch: center (" << PSI.p()(0) << "," << PSI.p()(1) << ") "
			<< "width=" << PSI.w() << vcl_endl 
			<< "fill_front:" << vcl_endl
			<< front_mat << vcl_endl
			<< "Valid pixels:" << vcl_endl
			<< valid_mat << vcl_endl
			<< "------------------------------------------" << vcl_endl;
#endif

	//front_mat.transpose(); // flip col and row

    int matSize = front_mat.rows();
    int dir = 0;
    // Start at center pixel of patch, need to find it's derivative
    int centIndex = winRad;

    // Restrict the search to 5x5 central of the patch
    // Since weighing-function will be used, anything beyond will get
    //  t-value of (x(t),y(t)) that is too far from t_center and contributes 
    //  little to the curve estimated at center. i.e. waste of computation
    int stepInward = ((winRad - 2) > 0? winRad - 2: 0);
    // vectors for x, y coordinates in paramater t=index
    // After parsing the coord, again we restrict to at most 7 pairs to save computation
    // Ideally from t = -3 to 3, but then here we use t = 0 to 7 for simplicity
    vnl_vector<int> XCoord(7);
    vnl_vector<int> YCoord(7);
    // Largest value of parameter t, and value of t for center pixel
    int maxT, centerT;

    maxT = parseCoordVer2(front_mat, XCoord, YCoord, centerT, stepInward);

#ifdef DEBUG320
    vcl_cerr
        << "------------------------------------------" << vcl_endl
        << "Parsed coordinates:" << vcl_endl
        << "X: " << XCoord << vcl_endl
        << "Y: " << YCoord << vcl_endl
        << "------------------------------------------" << vcl_endl;
#endif

    if (maxT == 0) 
    { // only center pixel is fill_front, can't compute normal
        return false; 
    }
    
	// Construct weighting function matrix
    vnl_matrix<double> weight_mat(maxT+1, maxT+1, 0.0);
    int tvalue = 0;
    for (int i = 0; i<= maxT; i++)
    {
        tvalue = i - centerT;
        weight_mat(i, i) = exp(- (tvalue * tvalue));
    }

#ifdef DEBUG320
    vcl_cerr
        << "------------------------------------------" << vcl_endl
        << "Weighing Matrix:" << vcl_endl
        << weight_mat << vcl_endl;
#endif

    // Construct coefficient matrix
    vnl_matrix<double> coeff_mat(maxT+1, 3);
    coeff_mat.set_column(0, 1);
    for (int i = 0; i <= maxT; i++)
    { 
        tvalue = i - centerT;
        coeff_mat(i, 1) = tvalue;
        coeff_mat(i, 2) = tvalue * tvalue / 2.0;
    }

#ifdef DEBUG320
    vcl_cerr
        //<< "------------------------------------------" << vcl_endl
        << "Coefficient Matrix:" << vcl_endl
        << coeff_mat << vcl_endl;
        //<< "------------------------------------------" << vcl_endl;
#endif

    // Construct (trim down) X, Y coord vector
    vnl_matrix<double> Xvalue(maxT+1, 1);
    vnl_matrix<double> Yvalue(maxT+1, 1);
    for (int i = 0; i <= maxT; i++)
    {
        Xvalue(i, 0) = XCoord(i);
        Yvalue(i, 0) = YCoord(i);
    }

#ifdef DEBUG320
    vcl_cerr
        //<< "------------------------------------------" << vcl_endl
        << "X-value Matrix:" << Xvalue.transpose() << vcl_endl
        << "Y-value Matrix:" << Yvalue.transpose() << vcl_endl;
        //<< "------------------------------------------" << vcl_endl;
#endif

    // Construct solution vector
    // ( x(0), dx(0)/dt, d2x(0)/dt2 ) and ( y(0), dy(0)/dt, d2y(0)/dt2 )
    vnl_matrix<double> Xsoln(3, 1);
    vnl_matrix<double> Ysoln(3, 1);
    vnl_matrix_inverse<double> inv_mat(weight_mat * coeff_mat);   
    Xsoln = inv_mat.solve(weight_mat * Xvalue);
	Ysoln = inv_mat.solve(weight_mat * Yvalue);
    double magnitude = sqrt(Ysoln(1,0) * Ysoln(1,0) + Xsoln(1,0) * Xsoln(1,0));
    
    if (magnitude == 0)
        return false;

#ifdef DEBUG320
    vcl_cerr
        //<< "------------------------------------------" << vcl_endl
        << "Xsoln: " << Xsoln.transpose() << vcl_endl
        << "Ysoln: " << Ysoln.transpose() << vcl_endl
        << "------------------------------------------" << vcl_endl;
#endif

    normal(0) = -Ysoln(1,0) / magnitude;
    normal(1) = Xsoln(1,0) / magnitude;
    
#ifdef DEBUG320
    vcl_cerr << "normal = " << normal(0) << ", " << normal(1) << vcl_endl;
#endif

    time_normal += TIMER_ELLAPSED;

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
	
    count_gradient++;
    TIMER_START;

	if (PSI.w() == 0) return false;

	int winRad = 1;
	int winSize = 2 * winRad + 1; // 2*1 + 1 = 3
	int patSize = PSI.sz();

	// Make copies of values
	vnl_matrix<int> valid_mat;
	vnl_matrix<int> unfilled_mat;
	vnl_matrix<int> grayscale_mat;

	PSI.get_pixels(inpainted_grayscale, grayscale_mat, valid_mat);
	PSI.get_pixels(unfilled, unfilled_mat, valid_mat);
 
	//grayscale_mat.transpose();
	//unfilled_mat.transpose();

    /* ***********************************
        This uses Sobel operator
    ************************************* */
    /*
    double XcurrGrad = 0;
    double YcurrGrad = 0;
    double XhighestGrad = 0;
    double YhighestGrad = 0;
    double highestGrad = 0;
    double currentGrad = 0;

    vnl_matrix<int> Xsobel(winSize, winSize, 1);
    vnl_matrix<int> Ysobel(winSize, winSize);

    Xsobel(0,0) = -1;
    Xsobel(0,1) = 0;
    Xsobel(0,2) = 1;
    Xsobel(1, 0) = -2;
    Xsobel(1, 1) = 0;
    Xsobel(1, 2) = 2;
    Xsobel(2, 0) = -1;
    Xsobel(2, 1) = 0;
    Xsobel(2, 2) = 1;

    Ysobel(0,0) = 1;
    Ysobel(0,1) = 2;
    Ysobel(0,2) = 1;
    Ysobel(1, 0) = 0;
    Ysobel(1, 1) = 0;
    Ysobel(1, 2) = 0;
    Ysobel(2, 0) = -1;
    Ysobel(2, 1) = -2;
    Ysobel(2, 2) = -1;

	// loop through patch with sliding window
	for (int x = 0; x < patSize; x++)
	{
		for (int y = 0; y < patSize; y++)
		{
            // Sobel to get gradient
            XcurrGrad = 0;
            YcurrGrad = 0;
            // Each sliding window (winRad*2 + 1)
            for (int i = -winRad; i <= winRad; i++)
            {
                for (int j = -winRad; j <= winRad; j++)
                {
                    if (x+i >= 0 && x+i < patSize && y+j >= 0 && y+j < patSize)
                        //&& valid_mat(x+i, y+j) && !unfilled_mat(x+i, y+j))
                    {
                        XcurrGrad += 
                            (Xsobel(i+winRad, j+winRad) * grayscale_mat(i+x, j+y));       
                        YcurrGrad += 
                            (Ysobel(i+winRad, j+winRad) * grayscale_mat(i+x, j+y));
                    }

                }
            } //finish individual sliding window

            currentGrad = sqrt(XcurrGrad*XcurrGrad + YcurrGrad*YcurrGrad);
            if (currentGrad > highestGrad)
            {
                highestGrad = currentGrad;
                XhighestGrad = XcurrGrad;
                YhighestGrad = YcurrGrad;
            }
		}
	}// finish whole patch

    grad(0) = XhighestGrad;
    grad(1) = YhighestGrad;
    */

    /*************************************************
        This simply get highest and lowest values and calculate gradient
    **************************************/

    int lowestCoord[2];
    int highestCoord[2];
	double highestDiff = -1;

    int lowCoord[2];
    int highCoord[2];
    int lowValue = 300;
	int highValue = -1;
    int pixelValue;
	// loop through patch with sliding window
	for (int x = 0; x < patSize - 2 * winRad; x++)
	{
		for (int y = 0; y < patSize - 2 * winRad; y++)
		{
			// Sliding window
            lowValue = 300;
	        highValue = -1;
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


    time_gradient += TIMER_ELLAPSED;
	return true;

	///////////////////////////////////////////////////////////
	//     DO NOT CHANGE ANYTHING BELOW THIS LINE            //
	///////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////
//     PLACE ANY ADDITIONAL CODE (IF ANY) HERE           //
///////////////////////////////////////////////////////////

// Get next direct adjacent pixel, store in same row/col parameters
// Return false if none
bool getAdjacentPixel(
    const vnl_matrix<int>& front_mat,
    int& row, int& col)
{
    int matDim = front_mat.rows();

    if (row-1 >= 0 && front_mat(row-1, col)) // up
        row --;
    else if (row+1 < matDim && front_mat(row+1, col)) // down
        row ++;
    else if (col+1 < matDim && front_mat(row, col+1)) //right
        col++;
    else if (col-1 >= 0 && front_mat(row, col-1)) // left
        col--;
    else
        return false;

    return true;
}

// Get next diagonal pixel, store in same row/col parameters
// Return false if none
bool getDiagonalPixel(
    const vnl_matrix<int>& front_mat,
    int& row, int& col)
{
    int matDim = front_mat.rows();

        // Check for diagonal
    if (row-1 >= 0 && col-1 >= 0 && front_mat(row-1, col-1)) // top-left
    {
        row--;
        col--; 
    }
    else if (row-1 >= 0 && col+1 < matDim && front_mat(row-1, col+1)) // top-right
    {
        row--;
        col++;
    }
    else if (row+1 < matDim && col-1 >= 0 && front_mat(row+1, col-1)) // bottom-left
    {
        row++;
        col--;
    }
    else if (row+1 < matDim && col+1 < matDim && front_mat(row+1, col+1)) // bottom-right
    {
        row++;
        col++;
    }
    else
        return false;

    return true;
}

// Parse x(t), y(t), return largest available index
// vector coordinate is stored from t = -2, -1, 0, 1, 2, etc
// index of t=0 for center pixel is returned in "tCenter"
int parseCoordVer2(
    const vnl_matrix<int>& front_mat,
    vnl_vector<int>& XCoord,
    vnl_vector<int>& YCoord,
    int& tCenter,
    int stepInward)
{
    // Different approach:
    // Start from center pixel, and find 2 outlet pixels
    // then trace outward from the outlet pixels
    int center;
    // Actualy central regoin of matrix that we cares..
    int matDim = front_mat.rows() - (2 * stepInward); 
    int maxIndex = -1;
    bool cont = true;
    int t_ctr = 0;
    int row1, col1;
    
    // local copy of matrix & storage vectors to mess around
    vnl_matrix<int> FrontMat(matDim, matDim);
    center = matDim / 2;
    vnl_vector<int> Xtemp(3);
    vnl_vector<int> Ytemp(3);

    // Copy the matrix values we need
    for (int i = 0; i < matDim; i++)
    {
        for (int j = 0; j < matDim; j++)
        {
            FrontMat(i, j) = front_mat(i+stepInward, j+stepInward);
        }
    }

#ifdef DEBUG320
    vcl_cerr
        << "Trimmed fill_front matrix: " << vcl_endl<< FrontMat << vcl_endl;
#endif
    // "Remove" center pixel and start tracing for max of 3 pixels
    FrontMat(center, center) = 0;
    row1 = center;
    col1 = center;
    for (int i = 0; cont && i < 3 ; i++)
    {
        if (!getAdjacentPixel(FrontMat, row1, col1))
        {
            cont = getDiagonalPixel(FrontMat, row1, col1);
        }
        if (cont)
        {
            Xtemp(i) = col1;
            Ytemp(i) = row1;
            maxIndex++;
            // "Remove" this pixel so wont' back-trace
            FrontMat(row1, col1) = 0;
        }
    } // done first 3 (or less) coord

    //copy value to main vectors
    for (int i = 0, tempIndex = maxIndex; i <= maxIndex; i++, tempIndex--)
    {
        XCoord(i) = Xtemp(tempIndex);
        YCoord(i) = Ytemp(tempIndex);
    }

    // Add center pixel coord
    maxIndex++;
    tCenter = maxIndex;
    XCoord(tCenter) = center;
    YCoord(tCenter) = center;

    // Get 3 more, with fresh front-matrix minus the previous outlet from center
    for (int i = 0; i < matDim; i++)
    {
        for (int j = 0; j < matDim; j++)
        {
            FrontMat(i, j) = front_mat(i+stepInward, j+stepInward);
        }
    }
    FrontMat(Ytemp(0), Xtemp(0)) = 0;
    FrontMat(center, center) = 0;
    row1 = center;
    col1 = center;
    cont = true;
    for (int i = 0; cont && i < 3 ; i++)
    {
        if (!getAdjacentPixel(FrontMat, row1, col1))
        {
            cont = getDiagonalPixel(FrontMat, row1, col1);
        }
        if (cont)
        {
            maxIndex++;
            XCoord(maxIndex) = col1;
            YCoord(maxIndex) = row1;
            
            // "Remove" this pixel so wont' back-trace
            FrontMat(row1, col1) = 0;
        }
    } 

    return maxIndex;
}

/////////////////////////////////////////////////////////

