#include <iostream>
#include <math.h>
#include "mex.h"
using namespace std;
double STlinkage(double *loc_list, int maxdist, mwSize n , mwSize m, double *finalarr )
{
  //Assumed input: loc_list, consisting of frame, x, y columns
  //
  // Output: particles array
  // columns: frame, Xpos, Ypos, Zps


  int ii = 0,i = 0 ,j = 0;

  double maxdist_secondaryDir = 1.0;
  double maxDistMultiplier = 1.5;
  int tryeverything = 0;
  int counter = 0;
  double *acclist =new double [n];
  double max1;
  int hcounter = 0;
  int dim1;
  int curframe;
  int modif1,modif2;
  double c;

  int i2 = 0;
  int j2 = 0;
  int s1;
  int s2;
  int s3;
  int s4;
  int correspondingrow ;
  int allassigned = 0;
  int  counterfr = 0;
  bool stuck = false;
 double *test = new double [n * m];

  int qw;

  for(ii = 0; ii < n ; ii++)
  {

      qw = 0;
      for(i = 0; i < n; i++)
        {
            j = 0;

          c = (ii + 1) / loc_list[i];
          if(c == 1.00)
            {
              qw = qw + 1;
            }


        }
        //cout<<qw<<endl;
        acclist[ii] = qw;


  }


   max1 = loc_list[0];

  for(i = 0; i < n; i++)
  {
    j = 0;
    if(loc_list[i + n * j] > max1)
    {
      max1 = loc_list[i + n * j];
    }
  }


  for(curframe = 0; curframe < max1 ; curframe++)
  {
     dim1 = acclist[curframe];
     double *DistanceX = new double [dim1*dim1];
     double *DistanceY = new double [dim1*dim1];
     int  *validcombinations = new int [dim1*dim1];
     double  *availablearr = new double[dim1*dim1];
     double  *availablearrstart = new double[dim1*dim1];
     int *v1 = new int [dim1*dim1];
     int *v2 = new int [dim1*dim1];
     int *v3 = new int [dim1*dim1];
     int *v4 = new int [dim1*dim1];
     double *tarr = new double [dim1 * m];
     double *framearr = new double [dim1 * (m + 1)];
     double *finalmidpointsforctPhasor = new double [3 * dim1];

    if(curframe == 0)
    {
      for(j = 0; j < m; j++)
      {
        for(i = 0; i < dim1; i++)
        {
          tarr[i + j * dim1] = loc_list[i + j * n];

        }

      }

    }
    if(curframe != 0)
    {
      modif1 = 0;
      modif2 = 0;
      for(i = 0; i < curframe ; i++)
        modif1 = modif1 + acclist[i];
      for(i = 0; i < curframe + 1 ; i++)
        modif2 = modif2 + acclist[i];
      for(j = 0  ; j < m  ; j++)
      {
        for(i = modif1 ; i < modif2; i++)
        {
          tarr[(i - (modif1)) + j * dim1] = loc_list[i + j * n];
        }
      }
    }


    for(j = 0; j < m + 1; j++)
    {
      for(i = 0; i < dim1  ; i++)
      {
        if( j == m)
        {
            framearr[i + dim1 * j] = 1;
        }
        if( j != m)
        {
            framearr[i  + dim1 * j] = tarr[i  + dim1 * j];
        }
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  //cout << "Length of array = " << (sizeof(framearr)/sizeof(*framearr)) <<endl;

  for(j2 = 0; j2 < dim1 ; j2++)
    {
      for(i2 = 0; i2 <dim1 ; i2++)
      {
          DistanceX[i2 + dim1 * j2] = fabs(framearr[i2 +  dim1 * 1] - framearr[j2 + dim1 * 1]) ;
          DistanceY[i2 + dim1 * j2] = fabs(framearr[i2 +  dim1 * 2] - framearr[j2 + dim1 * 2]) ;
          availablearr[i2 + dim1 * j2] = 1;
          if(i2 == j2)
          {
            availablearr[i2 + dim1 * j2] = 1;
          }
      }
    }
    allassigned = 0;
    counterfr = 0;
    stuck = false;
    for(j = 0; j < m; j++)
    {
      for(i = 0; i < dim1 ; i++)
      {
        finalmidpointsforctPhasor[i + dim1 * j] = 0;
      }
    }

  while (allassigned == 0 && dim1 > 0 && stuck == false)
  {
    for (j2 = 0; j2 < dim1; j2++)
    {
      for (i2 = 0; i2 < dim1; i2++)
      {
          if (DistanceX[i2 + dim1 * j2] < maxdist_secondaryDir)
          {
            v1[i2 + dim1 * j2] = 1.00;
          }
          else
          {
            v1[i2 + dim1 * j2] = 0.00;
          }
          if (DistanceX[i2 + dim1 * j2] < maxdist)
          {
            v2[i2 + dim1 * j2] = 1.00;
          }
          else
          {
            v2[i2 + dim1 * j2] = 0.00;
          }
          if (DistanceY[i2 + dim1 * j2] < maxdist_secondaryDir)
          {
            v3[i2 + dim1 * j2] = 1.00;
          }
          else
          {
            v3[i2 + dim1 * j2] = 0.00;
          }
          if (DistanceY[i2 + dim1 * j2] < maxdist)
          {
            v4[i2 + dim1 * j2] = 1.00;
          }
          else
          {
            v4[i2 + dim1 * j2] = 0.00;
          }
          validcombinations[i2 + dim1 * j2] = ((v1[i2 + dim1 * j2]*(v4[i2 + dim1 * j2]*maxDistMultiplier)) + (v3[i2 + dim1 * j2] *(v2[i2 + dim1 * j2] *maxDistMultiplier)))*availablearr[i2 + dim1 * j2];

      }

    }

  for (i = 0; i < dim1; i++)
    {
              validcombinations[i + dim1 * i] = 0;
    }
  //loop over columns
  for (j2 = 0; j2 < dim1; j2++)
  {
    for (i2 = 0; i2 < dim1; i2++)
    {
        availablearrstart[i2 + dim1 * j2] = availablearr[i2 + dim1 * j2];
        //cout<<availablearr[i2 * dim1 + j2]<<" ;";
    }
    //cout<<endl;
  }
  //cout<<"////////////////////////////////"<<endl;
  for (int combcol = 0; combcol < dim1 ; combcol++)
  {
    s1 = 0;
    for (i = 0; i < dim1; i++)
    {
      s1 = s1 + validcombinations[i + dim1 * combcol];
    }
    //Only one combination, mix them, and set them as
    //unavailable for next iteration
  if(s1 == 1)
    {

      for (i = 0; i < dim1; i++)
      {

        if (validcombinations[i + dim1 * combcol] == 1)
        {
          correspondingrow = i;

        }

      }
      //Add the middle point between them for phasor

          finalmidpointsforctPhasor[counterfr + dim1 * 0] = curframe + 1;

          finalmidpointsforctPhasor[counterfr + dim1 * 1] = (framearr[combcol + dim1 * 1] * 0.5) + (framearr[correspondingrow + dim1 * 1] * 0.5 );

          finalmidpointsforctPhasor[counterfr + dim1 * 2] = (framearr[combcol + dim1 * 2] * 0.5) + (framearr[correspondingrow + dim1 * 2] * 0.5 );

          counterfr = counterfr + 1;
          //Set them as unavailable

          for(i = 0; i < dim1; i++)
          {
            availablearr[i + dim1 * combcol] = 0;
            availablearr[i + dim1 * correspondingrow] = 0;
            availablearr[combcol + dim1 * i] = 0;
            availablearr[correspondingrow + dim1 * i] = 0;
          }
          break;
    }
  else if(s1 == 0)
    {
    s2 = 0;
    for (i2 = 0; i2 < dim1; i2++)
    {
      s2 = s2 + availablearr[i2 + dim1 * combcol];
    }
    if (s2 > 0)
    {
        //Add it for phasor
        finalmidpointsforctPhasor[counterfr + dim1 * 0] = curframe + 1;
        finalmidpointsforctPhasor[counterfr + dim1 * 1] = framearr[combcol + dim1 * 1];
        finalmidpointsforctPhasor[counterfr + dim1 * 2] = framearr[combcol + dim1 * 2];
        counterfr = counterfr+1;
        //Set it as unavailable
        for(i = 0; i < dim1; i++)
        {
          availablearr[i + dim1 * combcol] = 0;
          availablearr[combcol + dim1 * i] = 0;
        }


    }

  }
  //else
  //{

  //}
  s3 = 0;
  for (j2 = 0; j2 < dim1; j2++)
  {
    for (i2 = 0; i2 < dim1; i2++)
    {
        if (availablearr[i2 + dim1 * j2] == 0)
        {
          s3++;
        }
    }
  }
  if (s3 == dim1 * dim1)
  {
    allassigned = 1;
  }

  }
  s4 = 0;
  for (j2 = 0; j2 < dim1; j2++)
  {
    for (i2 = 0; i2 < dim1; i2++)
    {
        if (availablearr[i2 + dim1 * j2] == availablearrstart[i2 + dim1 * j2])
        {
          s4++;
        }
    }
  }
  if (s4 == dim1 * dim1)
  {
    stuck = true;
  }
  else
    stuck = false;


  }

  counter = counter + counterfr;
  for(j = 0; j < m; j++)
  {
    for( i = hcounter ; i < counter ; i++)
    {
      finalarr[i + n * j] = finalmidpointsforctPhasor [(i - hcounter )+ dim1 * j];
    }
  }
  hcounter = counter;
  delete [] tarr;
  delete [] framearr;
  }
  //finalarr(finalarr(:,1)==0,:) = [];
  // end function




}
//////////////////////////////////////////////////////////////////////////////////////
//Matlab Part/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *mat1;              /* input scalar */
    int md;               /* 1xN input matrix */
    size_t c1;                   /* size of matrix */
    size_t r1;
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1)
    {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    /* make sure the second input argument is scalar */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
      }


    /* check that number of rows in second input argument is 1 */
    //if(mxGetM(prhs[1])!=1) {
        //mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    //}




    /* create a pointer to the real data in the input matrix  */
    mat1 = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    c1 = mxGetN(prhs[0]);
    r1 = mxGetM(prhs[0]);
    /* get the value of the scalar input  */
    md = mxGetScalar(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)r1,(mwSize)c1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    STlinkage(mat1, md, (mwSize)r1, (mwSize)c1, outMatrix);
}
