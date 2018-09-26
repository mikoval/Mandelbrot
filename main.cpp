#include <iostream>
#include <thread>
#include <gmp.h>
#include <vector>
#include <algorithm>

#define FRACT_X "-0.9299546476636262932607706528493890456352598293371764563158011951410031896505225030586813656409098178717693332879518"
#define FRACT_Y "0.29188677067185720913132284864117522849541418009272893295153958202520477192287123714629376162463356657967494855049625"

double R_TARGET = 0;
#include "lodepng.h"

#include <math.h>

using namespace std;

#define MAX_FRAMES 50000
#define KEY_SIZE 50
#define START 0
#define START_AVG 0

int width2 = 2000, height2 = 2000;
int alias = 1;
int mult = 3;

int width = width2 * alias, height = height2 * alias;
int max_threads = 10;

static double total_r = 0.0;
static double total_r_counter = 0.0;

mpf_t mpf_x0;
mpf_t mpf_y0;
mpf_t mpf_xstart;
mpf_t mpf_ystart;
mpf_t mpf_tmp_x0;
mpf_t mpf_tmp_y0;
mpf_t mpf_tmp_x1;
mpf_t mpf_tmp_y1;
mpf_t mpf_tmp_x2;
mpf_t mpf_tmp_y2;
mpf_t mpf_tmp_x3;
mpf_t mpf_tmp_y3;
mpf_t mpf_tmp_x4;
mpf_t mpf_tmp_y4;

void interpolateImage(std::vector<double> &image1, std::vector<double> &image2, std::vector<double> &image3, float ITR_PARENT, float ITR, float MAX) {
	double count = ITR + (ITR_PARENT - KEY_SIZE);
	//count *= 2.0;
	double count2 = ITR_PARENT;
	double count3 = ITR_PARENT - KEY_SIZE;
	//count2 *= 2.0;
	//count3 *= 2.0;
	double mag = pow( 0.5, count / 60.0 );
	double mag2 = pow( 0.5, count2 / 60.0 );
	double mag3 = pow( 0.5, count3 / 60.0 );

	mpf_set_d(mpf_tmp_x1, mag );
	mpf_set_d(mpf_tmp_y1, mag );
	mpf_set_d(mpf_tmp_x2, 0.25 * width);
	mpf_set_d(mpf_tmp_y2, 0.25 * height);
	mpf_set_d(mpf_tmp_x4, 0.5 * width);
	mpf_set_d(mpf_tmp_y4, 0.5 * height);
	for (int x = 0; x < width; x++){
	    for (int y = 0; y < height; y++){

			mpf_set_d(mpf_tmp_x0, 4.0 * (((long double)x / (long double)width) - 0.5));
			mpf_set_d(mpf_tmp_y0, 4.0 * (((long double)y / (long double)height) - 0.5));
			double red = 0;
			double green = 0;
			double blue = 0;
            
			{


			    mpf_set_d(mpf_tmp_x3, 1.0 / mag2);
			    mpf_set_d(mpf_tmp_y3, 1.0 / mag2);
			    mpf_mul(mpf_tmp_x0, mpf_tmp_x0, mpf_tmp_x1);
			    mpf_mul(mpf_tmp_y0, mpf_tmp_y0, mpf_tmp_y1);

			    mpf_add(mpf_x0, mpf_tmp_x0, mpf_xstart);
			    mpf_add(mpf_y0, mpf_tmp_y0, mpf_ystart);

			    mpf_sub(mpf_x0, mpf_x0, mpf_xstart);
			    mpf_sub(mpf_y0, mpf_y0, mpf_ystart);

			    mpf_mul(mpf_x0, mpf_tmp_x2, mpf_x0);
			    mpf_mul(mpf_y0, mpf_tmp_y2, mpf_y0);


			    mpf_mul(mpf_x0, mpf_tmp_x3, mpf_x0);
			    mpf_mul(mpf_y0, mpf_tmp_y3, mpf_y0);


			    mpf_add(mpf_x0, mpf_tmp_x4, mpf_x0);
			    mpf_add(mpf_y0, mpf_tmp_y4, mpf_y0);

			}
			double ox = mpf_get_d(mpf_x0);
			double oy = mpf_get_d(mpf_y0);
			mpf_set_d(mpf_tmp_x3, 1.0 / mag3);
			mpf_set_d(mpf_tmp_y3, 1.0 / mag3);


			mpf_set_d(mpf_tmp_x0, 4.0 * (((long double)x / (long double)width) - 0.5));
			mpf_set_d(mpf_tmp_y0, 4.0 * (((long double)y / (long double)height) - 0.5));
            
			{


			    mpf_mul(mpf_tmp_x0, mpf_tmp_x0, mpf_tmp_x1);
			    mpf_mul(mpf_tmp_y0, mpf_tmp_y0, mpf_tmp_y1);

			    mpf_add(mpf_x0, mpf_tmp_x0, mpf_xstart);
			    mpf_add(mpf_y0, mpf_tmp_y0, mpf_ystart);

			    mpf_sub(mpf_x0, mpf_x0, mpf_xstart);
			    mpf_sub(mpf_y0, mpf_y0, mpf_ystart);

			    mpf_mul(mpf_x0, mpf_tmp_x2, mpf_x0);
			    mpf_mul(mpf_y0, mpf_tmp_y2, mpf_y0);


			    mpf_mul(mpf_x0, mpf_tmp_x3, mpf_x0);
			    mpf_mul(mpf_y0, mpf_tmp_y3, mpf_y0);


			    mpf_add(mpf_x0, mpf_tmp_x4, mpf_x0);
			    mpf_add(mpf_y0, mpf_tmp_y4, mpf_y0);

			}
			double ox2 = mpf_get_d(mpf_x0);
			double oy2 = mpf_get_d(mpf_y0);



            double xfract = ox - floor(ox);
            double yfract = oy - floor(oy);
            xfract = 1.0 - xfract;
            yfract = 1.0 - yfract;

            int xind = (int) ox;
            int yind = (int) oy;
            int xind2 = (int) ox + 1;
            int yind2 = (int) oy  + 1;

            std::vector<double> *iptr = &image1;

            bool skip = false;
            double i0 = 0.0;
            double r = 0.0;
            double i0_tmp = 0.0;
            double r_tmp = 0.0;
            if(xind < 0 || xind >= width ||
               xind2 < 0 || xind2 >= width || 
               yind < 0 || yind >= height || 
               yind2 < 0 || yind2 >= height) {

                skip = true;
                i0 = -1.0;
            }



			if(!skip){

			    float r1 = (*iptr)[yind * 4 * height + xind * 4 + 0];
			    float r2 = (*iptr)[yind * 4 * height + xind2 * 4 + 0];
			    float r3 = (*iptr)[yind2 * 4 * height + xind * 4 + 0];
			    float r4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 0];
			    float red1 = r1 * xfract + r2 * (1.0 - xfract);
			    float red2 = r3 * xfract + r4 * (1.0 - xfract);
			    float red3 = red1 * yfract + red2 * (1.0 - yfract);
			    red3 = r1;
			    if(r1 < 0 || r2 < 0 || r3 < 0 || r4 < 0){
			        red3 = -1.0;
			    }

			    float g1 = (*iptr)[yind * 4 * height + xind * 4 + 1];
			    float g2 = (*iptr)[yind * 4 * height + xind2 * 4 + 1];
			    float g3 = (*iptr)[yind2 * 4 * height + xind * 4 + 1];
			    float g4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 1];
			    float green1 = g1 * xfract + g2 * (1.0 - xfract);
			    float green2 = g3 * xfract + g4 * (1.0 - xfract);
			    float green3 = green1 * yfract + green2 * (1.0 - yfract);

			    float b1 = (*iptr)[yind * 4 * height + xind * 4 + 2];
			    float b2 = (*iptr)[yind * 4 * height + xind2 * 4 + 2];
			    float b3 = (*iptr)[yind2 * 4 * height + xind * 4 + 2];
			    float b4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 2];
			    float blue1 = b1 * xfract + b2 * (1.0 - xfract);
			    float blue2 = b3 * xfract + b4 * (1.0 - xfract);
			    float blue3 = blue1 * yfract + blue2 * (1.0 - yfract);


			  

			  //////////
			  //
			  //
			  //

			      i0 =  red3; 
			      r = green3; 
			      blue = blue3;
			}
                    /////////////////////////////////////////////////////////////////////

            xfract = ox2 - floor(ox2);
            yfract = oy2 - floor(oy2);
            xfract = 1.0 - xfract;
            yfract = 1.0 - yfract;

            xind = (int) ox2;
            yind = (int) oy2;
            xind2 = (int) ox2 + 1;
            yind2 = (int) oy2  + 1;


            //`cout << xind2 << ", " << yind2 <<  endl;
            

            iptr = &image2;
    
            float r1 = (*iptr)[yind * 4 * height + xind * 4 + 0];
            float r2 = (*iptr)[yind * 4 * height + xind2 * 4 + 0];
            float r3 = (*iptr)[yind2 * 4 * height + xind * 4 + 0];
            float r4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 0];
            float red1 = r1 * xfract + r2 * (1.0 - xfract);
            float red2 = r3 * xfract + r4 * (1.0 - xfract);
            float red3 = red1 * yfract + red2 * (1.0 - yfract);
            if(r1 < 0 || r2 < 0 || r3 < 0 || r4 < 0){
                red3 = -1.0;
            }

            float g1 = (*iptr)[yind * 4 * height + xind * 4 + 1];
            float g2 = (*iptr)[yind * 4 * height + xind2 * 4 + 1];
            float g3 = (*iptr)[yind2 * 4 * height + xind * 4 + 1];
            float g4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 1];
            float green1 = g1 * xfract + g2 * (1.0 - xfract);
            float green2 = g3 * xfract + g4 * (1.0 - xfract);
            float green3 = green1 * yfract + green2 * (1.0 - yfract);

            float b1 = (*iptr)[yind * 4 * height + xind * 4 + 2];
            float b2 = (*iptr)[yind * 4 * height + xind2 * 4 + 2];
            float b3 = (*iptr)[yind2 * 4 * height + xind * 4 + 2];
            float b4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 2];
            float blue1 = b1 * xfract + b2 * (1.0 - xfract);
            float blue2 = b3 * xfract + b4 * (1.0 - xfract);
            float blue3 = blue1 * yfract + blue2 * (1.0 - yfract);



			i0_tmp =  red3; 
			r_tmp = green3; 
			blue = blue3;


            if(xind < 0 || xind >= width ||
               xind2 < 0 || xind2 >= width || 
               yind < 0 || yind >= height || 
               yind2 < 0 || yind2 >= height) {

                    
            }
			double l = (double) ITR / (double) KEY_SIZE;
			//l = 0.0;

			
			
			if(i0 == -1.0  || r != r){
				i0 = i0_tmp;
				r = r_tmp;

 			}
             
            else if (i0_tmp == -1.0 || r_tmp != r_tmp){

            } else {
                 
                 i0 = i0 * l + i0_tmp * (1.0 - l);
                 r = r * l + r_tmp * (1.0 - l);
            }
            
			image3[y * 4 * height2 + x * 4 + 0] = i0; 
			image3[y * 4 * height2 + x * 4 + 1] = r; 
			image3[y * 4 * height2 + x * 4 + 2] = 0; 
			image3[y * 4 * height2 + x * 4 + 3] = 0; 
      }
  }
}

float computeAvg(std::vector<double> &image, float count){
	double total_i = 0;
	float max_itr =(100 + (0.000058 * pow((float)count * 2, 2.25)));
	vector<double> arr;
	for (int x = 0; x < width; x++){
        for (int y = 0; y < height; y++){
        	double i0 = image[y * 4 * height2 + x * 4 + 0];
        	double r = image[y * 4 * height2 + x * 4 + 1];


        	if(i0 != -1 && i0 == i0 && i0 < max_itr && i0 > 0){
			arr.push_back(i0);
        		total_i += i0;
        		count++;
        	}
        	
        }
    }
	sort(arr.begin(), arr.end());
	double val = arr[arr.size() * 0.95];
	return val;


}

void computeImage(std::vector<double> &image_in, std::vector<unsigned char> &image_out, double avg, double count){
	count = count;
	float max_itr =(100 + (0.000058 * pow((float)count * 2, 2.25)));
    for (int x = 0; x < width; x++){
        for (int y = 0; y < height; y++){
        	double i0 = image_in[y * 4 * height2 + x * 4 + 0];
        	double r = image_in[y * 4 * height2 + x * 4 + 1];




        	float red = 0;
        	float green = 0; 
        	float blue = 0;


			{
			    double  v = 1.0 * pow(1.5 + count/1300.0, 3.0);
			    
			 
			    
			   
			    float theta = 0.0;
			   
			    
			    r = r * v;

			    
			    r = 0.5 + 0.5 * sin(1.0 * r + count / 200.0);
			    

			    
			     r = r * mult;


			    int ind = (int((r - 1.0 + mult))) % mult;
			    int ind2 = ( int(r)) % mult;
			    
			    float frac = r - floor(r);
			  

	
		        
		        float reds[] =
		        {   0.5 + 0.5 * sin( 3.0 +  theta  + (i0)/ 90.0  + (count)/ 280.0),
		            0.0,
		            //0.5 + 0.5 * sin(5.0 +  theta + (i)/ 50.0  + (count)/ 90.0),
		            0.5 + 0.5 * sin(8.0 +  theta + (i0)/ 105.0  + (count)/ 300.0)
		        };
		        float greens[] =
		        {   0.5 + 0.5 * sin( 5.0 +  theta + (i0)/ 85.0 + (count)/ 180.0),
		            0.0,
		            //0.5 + 0.5 * sin(4.0 +  theta + (i)/ 75.0 + (count)/ 100.0),
		            0.5 + 0.5 * sin(6.0 +  theta + (i0)/ 95.0 + (count)/ 230.0)
		        };
		        float blues[] =
		        {   0.5 + 0.5 * sin(3.0 +  2.0 * theta + (i0)/ 90.0 + (count)/ 235.0),
		            0.0,
		            //0.5 + 0.5 * sin(7.0 +  2.0 * theta + (i)/ 20.0 + (count)/ 80.0),
		            0.5 + 0.5 * sin(0.0 +  2.0 * theta + (i0)/ 110.0 + (count)/ 205.0)
		        };
		        
		        
		        float r1x = reds[ind]; float r2x = reds[ind2];
		        float g1x = greens[ind]; float g2x = greens[ind2];
		        float b1x = blues[ind]; float b2x = blues[ind2];

		                  
		        red = r1x + (r2x - r1x) * frac;
		        green = g1x + (g2x - g1x) * frac;
		        blue = b1x + (b2x - b1x) * frac;


		   		
		        if(i0 <= 0.0 || i0 > max_itr - 0.0){
		            red = 0;
		            green = 0;
		            blue = 0;
		        }
		        


		        
			}             


			  red = red * 255;
			  green = green * 255;
			  blue = blue * 255;

			  double amt = avg;

			  
			  if(avg < 30){ avg = 30;}
			  if(i0 >  avg ){
			  	red = 0;
			  	green = 0;
			  	blue = 0;

			  }
			  

			
            image_out[y * 4 * height2 + x * 4 + 0] = red; 
            image_out[y * 4 * height2 + x * 4 + 1] = green; 
            image_out[y * 4 * height2 + x * 4 + 2] = blue; 
            image_out[y * 4 * height2 + x * 4 + 3] = 255;
        }
    }

}

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

void draw(std::vector<double>* image, int i, int max_threads, int width, int height,  double count ){

	float max_itr =(100 + (0.000058 * pow((float)(count + KEY_SIZE) * 2, 2.25)));


    mpf_set_default_prec(100 + (int)((float)count / 10.0));
    int start = (int)(  ((float)i / (float)max_threads) * (float)height);
    int end = (int)(((float)(i+1) / (float)max_threads) * (float)height);

    // cout << start << ", " << end << endl;


    mpf_t mpf_x, mpf_y;












    long double red = 0;
    long double green = 0;
    long double blue = 0;
    mpf_t mpf_x0;
    mpf_t mpf_y0;
    mpf_t mpf_xstart;
    mpf_t mpf_ystart;
    mpf_t mpf_tmp_x0;
    mpf_t mpf_tmp_y0;
    mpf_t mpf_tmp_x1;
    mpf_t mpf_tmp_y1;
    mpf_t mpf_zx;
    mpf_t mpf_zy;
    mpf_t mpf_finalx;
    mpf_t mpf_finaly;

    mpf_init (mpf_x0);           
    mpf_init (mpf_y0);
    mpf_init (mpf_xstart);           
    mpf_init (mpf_ystart);
    mpf_init (mpf_tmp_x0);           
    mpf_init (mpf_tmp_y0);
    mpf_init (mpf_tmp_x1);           
    mpf_init (mpf_tmp_y1);
    mpf_init (mpf_zx);
    mpf_init (mpf_zy);
    mpf_init (mpf_finalx);
    mpf_init (mpf_finaly);

    mpf_t mpf_tmpx;
    mpf_t mpf_tmpy;
    mpf_t mpf_2;
    mpf_t mpf_4;
    mpf_t mpf_length;
    mpf_t p1;
    mpf_t p2;
    mpf_t p3;
    mpf_t p5;
    mpf_t p6;

    mpf_init (mpf_tmpx);           
    mpf_init (mpf_tmpy);
    mpf_init (mpf_2);
    mpf_init (mpf_4);
    mpf_init (mpf_length);
    mpf_init (p1);           
    mpf_init (p2);
    mpf_init (p3);           
    mpf_init (mpf_tmpy);
    mpf_init (p5);           
    mpf_init (p6);

    mpf_init (mpf_zy);


    mpf_set_str(mpf_xstart,FRACT_X,10);
    mpf_set_str(mpf_ystart,FRACT_Y,10);

    mpf_set_str (mpf_2, "2.0", 10);
    mpf_set_str (mpf_4, "1000000.0", 10);
    for(int y = start; y < end; y++){
        //cout << y << endl;
        for(int x = 0; x < width; x++){





            //cout << endl<< " base value offset " << endl;
            //mpf_out_str(stdout,  10, 200, mpf_x0);



            mpf_set_d(mpf_tmp_x0, 4.0 * (((long double)x / (long double)width) - 0.5));
            mpf_set_d(mpf_tmp_y0, 4.0 * (((long double)y / (long double)height) - 0.5));


            double sx = 0.360240443437;
            double sy = -0.64131306106480317;
            double cx = sx + 4.0 * (((long double)x / (long double)width) - 0.5);
            double cy = sy + 4.0 * (((long double)y / (long double)height) - 0.5)  ;


            //cout << endl<< " offset " << endl;
            //mpf_out_str(stdout,  10, 200, mpf_tmp_x0);

            mpf_set_d(mpf_tmp_x1, pow( 0.5, count / 60.0 ) );
            mpf_set_d(mpf_tmp_y1, pow( 0.5, count / 60.0 ) );

            //cout << endl<< " scale " << endl;
            //mpf_out_str(stdout,  10, 200, mpf_tmp_x1);

            mpf_mul(mpf_tmp_x0, mpf_tmp_x0, mpf_tmp_x1);
            mpf_mul(mpf_tmp_y0, mpf_tmp_y0, mpf_tmp_y1);

            //cout << endl<< " scaled offset " << endl;
            //mpf_out_str(stdout,  10, 200, mpf_tmp_x0);

            mpf_add(mpf_x0, mpf_xstart, mpf_tmp_x0);
            mpf_add(mpf_y0, mpf_ystart, mpf_tmp_y0);

            //cout << endl<< " total " << endl;
            //mpf_out_str(stdout,  10, 200, mpf_tmp_x0);





            //std::cout << std::fixed << setprecision(40);
            //cout <<  x0 << endl;


            //cout << endl<< "DONE" << endl;










            mpf_set_str(mpf_zx,"0.0", 10);
            mpf_set_str(mpf_zy,"0.0", 10);





            double i;


            float skip = 0.0;// pow(count/100.0, 2.0);
            float sum = 0.0, sum2 = 0.0;

            for(i = 0.0; i < max_itr; i++){



                mpf_mul(p1, mpf_zx, mpf_zx);
                mpf_mul(p2, mpf_zy, mpf_zy);
                mpf_sub(mpf_tmpx, p1, p2);

                mpf_add(mpf_finalx, mpf_tmpx, mpf_x0);


                mpf_mul(p3, mpf_zx, mpf_zy);
                mpf_add(p3, p3, p3);

                mpf_add(mpf_finaly, p3, mpf_y0);

                mpf_set (mpf_zy, mpf_finaly);
                mpf_set (mpf_zx, mpf_finalx);






                double zx = mpf_get_d(mpf_zx);
                double zy = mpf_get_d(mpf_zy);

                sum2 = sum;


                double tx = mpf_get_d(mpf_finalx);
                double ty = mpf_get_d(mpf_finaly);
                float mp=sqrt(tx * tx + ty * ty);

                float m = abs(mp  - sqrt(cx * cx + cy * cy)  );
                float M = mp + sqrt(cx * cx + cy * cy);

                double curve1 = 0.5 + 0.5 * sin(3.0 * atan2(zx, zy));
                int ind = (int) i;
                sum += curve1;









                mpf_mul(p5, mpf_zx, mpf_zx);
                mpf_mul(p6, mpf_zy, mpf_zy);
                mpf_add(mpf_length, p5, p6);
                if( mpf_cmp(mpf_length, mpf_4) > 0){
                    break;
                }



            }


            double zx = mpf_get_d(mpf_zx);
            double zy = mpf_get_d(mpf_zy);



            sum = sum / (i);








            sum2 = sum2 / (i - 1.0);




            




            double i2 = i + 1.0 + 1.0/log(2.0) * log(log(1000.0)/ log(sqrt(zx*zx + zy*zy)));
            
            if(i2 == i2){ i = i2; }

            


            double d = i - floor(i);


            double r = (sum * d + sum2 * (1.0 - d));

            if(r == r){
                total_r += r;
                total_r_counter++;

            }
            if(r != r){
            	r = 0;
            }
            

            r = r - R_TARGET;


            double r0 = r;




            double  v = 1.0 * pow(1.5 + count/500.0, 4.0);




            float theta = atan2(cx - sx, cy - sy);
            theta = 0.0;

            float t2 = 0.5  * sin(2.0 * theta);
            theta = sin(1.0 * theta);


            r = r * v;





            r = 0.5 + 0.5 * sin(1.0 * r + count / 100.0 + t2 );




            r = r * mult;





            int ind = (int((r - 1.0 + mult))) % mult;
            int ind2 = ( int(r)) % mult;

            float frac = r - floor(r);






            float reds[] =
            {   0.5 + 0.5 * sin( 3.0 +  theta  + (i)/ 40.0  + (count)/ 280.0),
                0.0,
                //0.5 + 0.5 * sin(5.0 +  theta + (i)/ 50.0  + (count)/ 90.0),
                0.5 + 0.5 * sin(8.0 +  theta + (i)/ 35.0  + (count)/ 300.0)
            };
            float greens[] =
            {   0.5 + 0.5 * sin( 5.0 +  theta + (i)/ 25.0 + (count)/ 180.0),
                0.0,
                //0.5 + 0.5 * sin(4.0 +  theta + (i)/ 25.0 + (count)/ 100.0),
                0.5 + 0.5 * sin(6.0 +  theta + (i)/ 35.0 + (count)/ 230.0)
            };
            float blues[] =
            {   0.5 + 0.5 * sin(3.0 +  2.0 * theta + (i)/ 30.0 + (count)/ 235.0),
                0.0,
                //0.5 + 0.5 * sin(7.0 +  2.0 * theta + (i)/ 20.0 + (count)/ 80.0),
                0.5 + 0.5 * sin(0.0 +  2.0 * theta + (i)/ 25.0 + (count)/ 205.0)
            };


            float r1 = reds[ind]; float r2 = reds[ind2];
            float g1 = greens[ind]; float g2 = greens[ind2];
            float b1 = blues[ind]; float b2 = blues[ind2];

            red = r1 + (r2 - r1) * frac;
            green = g1 + (g2 - g1) * frac;
            blue = b1 + (b2 - b1) * frac;





            if(r0 != r0){
                r0 = 0;
            }

            (*image)[y*4*height + x*4 + 0] = i;//(int) (255.0 * red);
            (*image)[y*4*height + x*4 + 1] = r0;//(int) (255.0 * green);
            (*image)[y*4*height + x*4 + 2] = (int) (255.0 * blue);
            (*image)[y*4*height + x*4 + 3] = 255;
        }
    }

}

void newCoord(double *ret, double x, double y, double mag, double mag2){
    /*
       double sx = 0.360240443437;
       double sy = -0.64131306106480317;

       double cx = sx + mag * 4.0 * (((long double)x / (long double)width) - 0.5);
       double cy = sy + mag * 4.0 * (((long double)y / (long double)height) - 0.5);

       double ox = 0.25 * width * (cx - sx) / mag2 + 0.5 * width;
       double oy = 0.25 * height * (cy - sy) / mag2 + 0.5 * height;
       */

    mpf_t mpf_x0;
    mpf_t mpf_y0;
    mpf_t mpf_xstart;
    mpf_t mpf_ystart;
    mpf_t mpf_tmp_x0;
    mpf_t mpf_tmp_y0;
    mpf_t mpf_tmp_x1;
    mpf_t mpf_tmp_y1;
    mpf_t mpf_tmp_x2;
    mpf_t mpf_tmp_y2;
    mpf_t mpf_tmp_x3;
    mpf_t mpf_tmp_y3;
    mpf_t mpf_tmp_x4;
    mpf_t mpf_tmp_y4;

    mpf_init (mpf_x0);           
    mpf_init (mpf_y0);
    mpf_init (mpf_xstart);           
    mpf_init (mpf_ystart);
    mpf_init (mpf_tmp_x0);           
    mpf_init (mpf_tmp_y0);
    mpf_init (mpf_tmp_x1);           
    mpf_init (mpf_tmp_y1);
    mpf_init (mpf_tmp_x2);           
    mpf_init (mpf_tmp_y2);
    mpf_init (mpf_tmp_x3);           
    mpf_init (mpf_tmp_y3);
    mpf_init (mpf_tmp_x4);
    mpf_init (mpf_tmp_y4);
    //
    ////////////////////
    mpf_set_str(mpf_xstart,FRACT_X,10);
    mpf_set_str(mpf_ystart,FRACT_Y,10);
    mpf_set_d(mpf_tmp_x0, 4.0 * (((long double)x / (long double)width) - 0.5));
    mpf_set_d(mpf_tmp_y0, 4.0 * (((long double)y / (long double)height) - 0.5));
    mpf_set_d(mpf_tmp_x1, mag );
    mpf_set_d(mpf_tmp_y1, mag );
    mpf_set_d(mpf_tmp_x2, 0.25 * width);
    mpf_set_d(mpf_tmp_y2, 0.25 * height);
    mpf_set_d(mpf_tmp_x3, 1.0 / mag2);
    mpf_set_d(mpf_tmp_y3, 1.0 / mag2);
    mpf_set_d(mpf_tmp_x4, 0.5 * width);
    mpf_set_d(mpf_tmp_y4, 0.5 * height);

    {


        mpf_mul(mpf_tmp_x0, mpf_tmp_x0, mpf_tmp_x1);
        mpf_mul(mpf_tmp_y0, mpf_tmp_y0, mpf_tmp_y1);

        mpf_add(mpf_x0, mpf_tmp_x0, mpf_xstart);
        mpf_add(mpf_y0, mpf_tmp_y0, mpf_ystart);

        mpf_sub(mpf_x0, mpf_x0, mpf_xstart);
        mpf_sub(mpf_y0, mpf_y0, mpf_ystart);

        mpf_mul(mpf_x0, mpf_tmp_x2, mpf_x0);
        mpf_mul(mpf_y0, mpf_tmp_y2, mpf_y0);


        mpf_mul(mpf_x0, mpf_tmp_x3, mpf_x0);
        mpf_mul(mpf_y0, mpf_tmp_y3, mpf_y0);


        mpf_add(mpf_x0, mpf_tmp_x4, mpf_x0);
        mpf_add(mpf_y0, mpf_tmp_y4, mpf_y0);

    }
    double ox = mpf_get_d(mpf_x0);
    double oy = mpf_get_d(mpf_y0);
    ret[0] = ox;
    ret[1] = oy;

}

int main(){
    mpf_set_default_prec(300);


    mpf_init (mpf_x0);           
    mpf_init (mpf_y0);
    mpf_init (mpf_xstart);           
    mpf_init (mpf_ystart);
    mpf_init (mpf_tmp_x0);           
    mpf_init (mpf_tmp_y0);
    mpf_init (mpf_tmp_x1);           
    mpf_init (mpf_tmp_y1);
    mpf_init (mpf_tmp_x2);           
    mpf_init (mpf_tmp_y2);
    mpf_init (mpf_tmp_x3);           
    mpf_init (mpf_tmp_y3);
    mpf_init (mpf_tmp_x4);
    mpf_init (mpf_tmp_y4);
    //
    ////////////////////
    mpf_set_str(mpf_xstart,FRACT_X,10);
    mpf_set_str(mpf_ystart,FRACT_Y,10);


    //////////////////////////////////////////////////
    std::vector<double> image1;
    image1.resize(width * height * 4);


    std::vector<double> image2;
    image2.resize(width * height * 4);

    std::vector<double> image3;
    image3.resize(width * height * 4);

    std::vector<unsigned char> image;
    image.resize(width * height * 4);

    std::thread myThreads[max_threads];


    float max_avg = START_AVG;
    float current_avg = START_AVG;

	total_r = 0;
	total_r_counter = 0;

        for (int i=0; i<max_threads; i++){
            myThreads[i] = std::thread(draw, &image1, i, max_threads, width, height, (double)10000 / 2.0);

        }
        for (int i=0; i<max_threads; i++){
            myThreads[i].join();
        }

	R_TARGET = total_r / total_r_counter;

    for(int j = START; j < MAX_FRAMES; j += KEY_SIZE){


        std::vector<double> tmp = image2;
        image2 = image1;
        image1 = tmp;

        for (int i=0; i<max_threads; i++){
            myThreads[i] = std::thread(draw, &image1, i, max_threads, width, height, (double)j);

        }
        for (int i=0; i<max_threads; i++){
            myThreads[i].join();
        }



        static bool s = false;
        if(s == false){
            s = true;
            continue;
        }

        for(double i = .5; i < KEY_SIZE; i++){
            double iterations = i + (double)(j - KEY_SIZE);
            string filename = "pictures/scene"  + to_string((int)iterations) + ".png";

            interpolateImage(image1, image2, image3, j, i, KEY_SIZE);

            float avg = computeAvg(image3, iterations);

	    if(avg > max_avg) { max_avg = avg; }

	    float diff_avg = max_avg - current_avg;
	    current_avg += diff_avg/50.0;
	    //current_avg = max_avg;
            computeImage(image3, image, current_avg, iterations);


            encodeOneStep((const char*)filename.data(), image, width2, height2);
        }
	cout << "AVG = " << current_avg << endl;

    }




    return 0;
}
