/*
 * motionblur.cc
 * Submission for assignment 3
 * CPSC 4310 - Image Processing
 *
 * Usage:
 *
 *   motionblur input.pgm output.pgm a b T
 *
 * Joshua Enns
 *
 */

#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "fft.h"

const double PI = acos(-1.0)/2;

using namespace std;

tuple **read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &newpam, tuple **array);
tuple **pad_image(tuple **array, const pam &inpam, pam &newpam);
Coeff **translate_image(Coeff **complexArray, const pam &newpam);
Coeff **image_to_coeff(tuple **array, const pam &newpam);
tuple **coeff_to_image(Coeff **complexArray, const pam &newpam);
Coeff **fourier_transform(Coeff **complexArray, const pam &newpam);
Coeff **inverse_fourier_transform(Coeff **complexArray, const pam &newpam);
Coeff **blur_image(Coeff **complexArray,const pam &newpam, double a, double b, double T);
tuple **crop_image(tuple **array, const pam &inpam, pam &newpam);
int nearest_pow_two(int current);

int main(int argc, char *argv[])
{
  /* variables for storing input parameters */
  double a = atof(argv[3]);
  double b = atof(argv[4]);
  double T = atof(argv[5]);

  /* structures for images */
  pam inpam;
  pam newpam;

  /* dynamic two dimensional arrays to store images */
  tuple **array;
  tuple **newarray;
  tuple **croppedarray;

  /* array of coefficients */
  Coeff **complexArray;

  /* initializes the library */
  pm_init(argv[0], 0);

  /* read the image */
  array = read_image(argv[1], inpam);

  /* set new image properties */
  newpam = inpam;

  /* pad the image */
  newarray = pad_image(array, inpam, newpam);

  /* convert image to coefficients */
  complexArray = image_to_coeff(newarray, newpam);

  /* translate image to bring corner pixels to the center */
  complexArray = translate_image(complexArray, newpam);

  /* perform fourier on coefficient array */
  complexArray = fourier_transform(complexArray, newpam);

  /* blur the image */
  complexArray = blur_image(complexArray, newpam, a, b, T);

  /* perform inverse fourier on coefficient array */
  complexArray = inverse_fourier_transform(complexArray, newpam);

  /* translate image to bring corner pixels to the center */
  complexArray = translate_image(complexArray, newpam);

  /* convert coefficients to image */
  newarray = coeff_to_image(complexArray, newpam);

  /* crop image back to original size */
  croppedarray = crop_image(newarray, inpam, newpam);

  /* write the output */
  write_image(argv[2], inpam, croppedarray);

  /* clean up */
  pnm_freepamarray(array, &inpam);
  pnm_freepamarray(newarray, &newpam);
  pnm_freepamarray(croppedarray, &inpam);
  return 0;
}

/* reads the image into the netpbm structure */
tuple **read_image(char *filename, pam &inpam)
{
  FILE *f;
  tuple **A;

  if ((f = pm_openr(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for reading." << endl;
    exit(1);
  }

  if ((A = pnm_readpam(f, &inpam, PAM_STRUCT_SIZE(tuple_type))) == NULL) {
    cerr << "Cannot read image \"" << filename << "\"." << endl;
    exit(1);
  }
  
  pm_close(f);
  return A;
}

/* writes the image to the given file */
void write_image(char *filename, const pam &newpam, tuple **array)
{
  FILE *f;
  pam outpam = newpam;

  if ((f = pm_openw(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for writing." << endl;
    exit(1);
  }

  outpam.file = f;
  pnm_writepam(&outpam, array);

  pm_close(f);
}

/* pads the image */
tuple **pad_image(tuple **array, const pam &inpam, pam &newpam)
{
  int row, col;

  //set new height and width to be the nearest power of 2
  newpam.height =  nearest_pow_two(inpam.height*2);
  newpam.width = nearest_pow_two(inpam.width*2);

  //create new array to store padded image
  tuple **newarray;
  newarray = pnm_allocpamarray(&newpam);

  //initialize image
  for (row = 0; row < newpam.height; row++) {
    for (col = 0; col < newpam.width; col++) {
        newarray[row][col][0] = 0;
    }
  }

  //pad the image
  for (row = 0; row < inpam.height; row++) {
    for (col = 0; col < inpam.width; col++) {
        newarray[row][col][0] = array[row][col][0];
    }
  }
  return newarray;
}

/* find first power of two that is greater than double the input */
int nearest_pow_two(int current)
{
  int powTwo = 2;
  while (powTwo < (current*2))
  {
    powTwo *= 2;
  }
  return powTwo;
}

/* crop image back to original size */
tuple **crop_image(tuple **array, const pam &inpam, pam &newpam)
{
  int row, col;

  //create new array to store cropped image
  tuple **newarray;
  newarray = pnm_allocpamarray(&inpam);

  //copy padded array to cropped array
  for (row = 0; row < inpam.height; row++) {
    for (col = 0; col < inpam.width; col++) {
      newarray[row][col][0] = array[row][col][0];
    }
  }
  return newarray;
}

/* translate image to center the outside corners */
Coeff **translate_image(Coeff **array, const pam &newpam)
{
  int row, col;
  for (row = 0; row < newpam.height; row++) {
    for (col = 0; col < newpam.width; col++) {
      //multiply each coefficient by -1^(x+y)
      array[row][col] = array[row][col] * pow(-1,(row+col));
    }
  }
  return array;
}

/* converts 2d image array into coefficient array */
Coeff **image_to_coeff(tuple **array, const pam &newpam)
{
  int row, col;

  //create array of length height to store complex numbers
  Coeff **complexArray = new Coeff*[newpam.height];

  for (row = 0; row < newpam.height; row++) {
    //create pointers to array rows
    complexArray[row] = new Coeff[newpam.width];
    for (col = 0; col < newpam.width; col++) {
      //convert existing integers to complex numbers and store them in the array
      complexArray[row][col] = Coeff(static_cast<double>(array[row][col][0]));
    }
  }
  return complexArray;
}

/* converts coefficient array to 2d image array */
tuple **coeff_to_image(Coeff **complexArray, const pam &newpam)
{
  int row, col;
  //create array to hold pixel values
  tuple **array;
  array = pnm_allocpamarray(&newpam);

  for (row = 0; row < newpam.height; row++) {
    for (col = 0; col < newpam.width; col++) {
      //check pixel boundaries
      if(static_cast<int>(complexArray[row][col].real()) < 0){
        //set value to 0 if negative
        array[row][col][0] = 0;
      } else if (static_cast<int>(complexArray[row][col].real()) > 255){
        //set value to 255 if above 255
        array[row][col][0] = 255;
      } else {
        //convert complex number to int and copy to int array
        array[row][col][0] = static_cast<int>(complexArray[row][col].real());
      }
    }
  }
  //clean up memory used by complex array
  for (row = 0; row < newpam.height; row++) {
    delete [] complexArray[row];
  }
  delete [] complexArray;

  return array;
}

/* perform fourier transform on image */
Coeff **fourier_transform(Coeff **complexArray, const pam &newpam)
{
  int row, col;
  //create buffer, set to the max of height and width
  Array buffer = new Coeff[max(newpam.width,newpam.height)];

  //perform 1D fourier transform on rows
  for (row = 0; row < newpam.height; row++) {
    FFT(complexArray[row], newpam.width, buffer);
  }

  //create array to hold columns
  Coeff **complexColArray = new Coeff*[newpam.width];

  //populate column array and perform 1D fourier transform
  for (col = 0; col < newpam.width; col++) {
    //create column array
    complexColArray[col] = new Coeff[newpam.height];
    for (row = 0; row < newpam.height; row++) {
      //populate column array by copying values from row array
      complexColArray[col][row] = complexArray[row][col];
    }
    //perform 1D fourier transform
    FFT(complexColArray[col], newpam.height, buffer);
  }

  //copy column array to return array
  for (row = 0; row < newpam.height; row++) {
    for (col = 0; col < newpam.width; col++) {
      complexArray[row][col] = complexColArray[col][row];
    }
  }
  //clean up memory used by complex column array
  for (row = 0; row < newpam.width; row++) {
    delete [] complexColArray[row];
  }
  delete [] complexColArray;

  return complexArray;
}

/* perform inverse fourier transform on image */
Coeff **inverse_fourier_transform(Coeff **complexArray, const pam &newpam)
{
  int row, col;
  //create buffer, set to the max of height and width
  Array buffer = new Coeff[max(newpam.width,newpam.height)];

  //perform 1D inverse fourier transform on rows
  for (row = 0; row < newpam.height; row++) {
    inverseFFT(complexArray[row], newpam.width, buffer);
  }

  //create array to hold columns
  Coeff **complexColArray = new Coeff*[newpam.width];

  //populate column array and perform 1D inverse fourier transform
  for (col = 0; col < newpam.width; col++) {
    //create column array
    complexColArray[col] = new Coeff[newpam.height];
    for (row = 0; row < newpam.height; row++) {
      //populate column array by copying values from row array
      complexColArray[col][row] = complexArray[row][col];
    }
    //perform 1D inverse fourier transform
    inverseFFT(complexColArray[col], newpam.height, buffer);
  }

  //copy column array to return array
  for (row = 0; row < newpam.height; row++) {
    for (col = 0; col < newpam.width; col++) {
      complexArray[row][col] = complexColArray[col][row];
    }
  }
  //clean up memory used by complex column array
  for (row = 0; row < newpam.width; row++) {
    delete [] complexColArray[row];
  }
  delete [] complexColArray;

  return complexArray;
}

/* blur image by modifying coefficients using degredation function*/
Coeff **blur_image(Coeff **complexArray, const pam &newpam, double a, double b, double T)
{
  int u, v;
  //calculate u and v modifications for H()
  int subu = newpam.height/2;
  int subv = newpam.width/2;
  
  //loop through entire array
  for (u = 0; u < newpam.height; u++) {
    for (v = 0; v < newpam.width; v++) {
      //store demoninator
      //demoninator = pi (ua + vb) 
      double denominator = PI * (((u-subu) * a) + ((v-subv) * b));
      double sinpart = sin(denominator);

      //calculate and store e^-jpi(ua+vb)
      Coeff omega = exp(Coeff(0, -denominator));

      //check to make sure it won't divide by 0
      if (denominator == 0.0){
        denominator = 1;
        sinpart = 1;
        omega = 1;
      }
      
      //calculate new coefficient
      //G(u,v) = F(u,v) * H(u,v) 
      //where F(u,v) is the calculated fourier coefficient
      //and H(u,v) is (T/pi(ua+vb)) * sin(pi(ua+vb)) * e^-jpi(ua+vb)
      complexArray[u][v] = complexArray[u][v] * (T/denominator) * sinpart * omega;
    }
  }
  return complexArray;
}
