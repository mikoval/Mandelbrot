#include <iostream>
#include <thread>
#include <gmp.h>


#include "lodepng.h"

#include <math.h>

using namespace std;

#define MAX_FRAMES 20000
#define KEY_SIZE 30
#define START 0

int width2 = 1000, height2 = 1000;
int alias = 1;
int mult = 3;

int width = width2 * alias, height = height2 * alias;
int max_threads = 10;

static double total_r = 0.0;
static double total_r_counter = 0.0;
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

void draw(std::vector<double>* image, int i, int max_threads, int width, int height,  double count ){
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
    

        mpf_set_str(mpf_xstart,"0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928",10);
        mpf_set_str(mpf_ystart,"-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991",10);
                
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
        
                
                float max_itr = (80 + (0.00000185 * pow((float)count, 2.51)));
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
                
                
                
                
               
                i = i + 1.0 + 1.0/log(2.0) * log(log(1000.0)/ log(sqrt(zx*zx + zy*zy)));
                double d = i - floor(i);
                
                
                double r = (sum * d + sum2 * (1.0 - d));
                
                if(r == r){
                    total_r += r;
                    total_r_counter++;

                }
                
                r = r - .65414045;


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
            
                

          
    
                if(i > max_itr){
                    red = 0;
                    green = 0;
                    blue = 0;
                    i = -1.0;
                }
                red = red ;
                green = green ;
                blue = blue;
                

                

                if(r != r){
                    i = 0;

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
        mpf_set_str(mpf_xstart,"0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928",10);
        mpf_set_str(mpf_ystart,"-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991",10);
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
		mpf_set_str(mpf_xstart,"0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928",10);
		mpf_set_str(mpf_ystart,"-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991",10);
	  std::vector<double> image1;
	  image1.resize(width * height * 4);
	  

      std::vector<double> image2;
      image2.resize(width * height * 4);

      	  std::vector<unsigned char> image;
      	  image.resize(width * height * 4);
	  
	  std::thread myThreads[max_threads];


        
	  for(int j = START; j < MAX_FRAMES; j += KEY_SIZE){


          std::vector<double> tmp = image2;
          image2 = image1;
          image1 = tmp;


		  for (int i=0; i<max_threads; i++){
			  myThreads[i] = std::thread(draw, &image1, i, max_threads, width, height, (double)j / 2.0);

		  }
		  for (int i=0; i<max_threads; i++){
			  myThreads[i].join();
		  }

          if(j == 0){
            continue;
          }

          for(int i = 0; i < KEY_SIZE; i++){
              int iterations = i + (j - KEY_SIZE);
              string filename = "pictures/scene"  + to_string(iterations) + ".png";



              double count = iterations;
	      count /= 2.0;
              double count2 = j;
              double count3 = j - KEY_SIZE;
	      count2 /= 2.0;
	      count3 /= 2.0;
              double mag = pow( 0.5, count / 60.0 );
              double mag2 = pow( 0.5, count3 / 60.0 );
              double mag3 = pow( 0.5, count2 / 60.0 );


                mpf_set_d(mpf_tmp_x1, mag );
                mpf_set_d(mpf_tmp_y1, mag );
                mpf_set_d(mpf_tmp_x2, 0.25 * width);
                mpf_set_d(mpf_tmp_y2, 0.25 * height);
                mpf_set_d(mpf_tmp_x3, 1.0 / mag2);
                mpf_set_d(mpf_tmp_y3, 1.0 / mag2);
                mpf_set_d(mpf_tmp_x4, 0.5 * width);
                mpf_set_d(mpf_tmp_y4, 0.5 * height);
              for (int x = 0; x < width; x++){
                  for (int y = 0; y < height; y++){
                      //cout << "INDEX: " << x << ", " << y << endl;

                mpf_set_d(mpf_tmp_x0, 4.0 * (((long double)x / (long double)width) - 0.5));
                mpf_set_d(mpf_tmp_y0, 4.0 * (((long double)y / (long double)height) - 0.5));
                      double red = 0;
                      double green = 0;
                      double blue = 0;
                        
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

                        if(xind < 0 || xind >= width ||
                           xind2 < 0 || xind2 >= width || 
                           yind < 0 || yind >= height || 
                           yind2 < 0 || yind2 >= height) {

                                continue;
                        }





			
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


                      

                      //////////
                      //
                      //
                      //

                          double i0 =  red3; 
                          double r = green3; 
                          blue = blue3;
                                /////////////////////////////////////////////////////////////////////

                        xfract = ox2 - floor(ox2);
                        yfract = oy2 - floor(oy2);
                        xfract = 1.0 - xfract;
                        yfract = 1.0 - yfract;

                        xind = (int) ox;
                        yind = (int) oy;
                        xind2 = (int) ox + 1;
                        yind2 = (int) oy  + 1;


                        iptr = &image2;
                
                        r1 = (*iptr)[yind * 4 * height + xind * 4 + 0];
                        r2 = (*iptr)[yind * 4 * height + xind2 * 4 + 0];
                        r3 = (*iptr)[yind2 * 4 * height + xind * 4 + 0];
                        r4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 0];
                        red1 = r1 * xfract + r2 * (1.0 - xfract);
                        red2 = r3 * xfract + r4 * (1.0 - xfract);
                        red3 = red1 * yfract + red2 * (1.0 - yfract);
                        if(r1 < 0 || r2 < 0 || r3 < 0 || r4 < 0){
                            red3 = -1.0;
                        }

                        g1 = (*iptr)[yind * 4 * height + xind * 4 + 1];
                        g2 = (*iptr)[yind * 4 * height + xind2 * 4 + 1];
                        g3 = (*iptr)[yind2 * 4 * height + xind * 4 + 1];
                        g4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 1];
                        green1 = g1 * xfract + g2 * (1.0 - xfract);
                        green2 = g3 * xfract + g4 * (1.0 - xfract);
                        green3 = green1 * yfract + green2 * (1.0 - yfract);

                        b1 = (*iptr)[yind * 4 * height + xind * 4 + 2];
                        b2 = (*iptr)[yind * 4 * height + xind2 * 4 + 2];
                        b3 = (*iptr)[yind2 * 4 * height + xind * 4 + 2];
                        b4 = (*iptr)[yind2 * 4 * height + xind2 * 4 + 2];
                        blue1 = b1 * xfract + b2 * (1.0 - xfract);
                        blue2 = b3 * xfract + b4 * (1.0 - xfract);
                        blue3 = blue1 * yfract + blue2 * (1.0 - yfract);



                          double i0_tmp =  red3; 
                          double r_tmp = green3; 
                          blue = blue3;


                        if(xind < 0 || xind >= width ||
                           xind2 < 0 || xind2 >= width || 
                           yind < 0 || yind >= height || 
                           yind2 < 0 || yind2 >= height) {

                                
                        }
			else {
                          //i0 = i0_tmp;
                          //r = r_tmp;
			
			}

			  
			  ///////////////////////////////////////////////////////////////////////
                
			  {
                        double  v = 1.0 * pow(1.5 + count/1300.0, 4.0);
                        
                     
                        
                       
                        float theta = 0.0;
                        
                        
                        r = r * v;

                        
                        r = 0.5 + 0.5 * sin(1.0 * r + count / 200.0);
                        
                   
                        
                         r = r * mult;

                  
                        int ind = (int((r - 1.0 + mult))) % mult;
                        int ind2 = ( int(r)) % mult;
                        
                        float frac = r - floor(r);
                      

                ///////////////
                ///
                            float max_itr =(80 + (0.00000185 * pow((float)count, 2.51)));
                            float reds[] =
                            {   0.5 + 0.5 * sin( 3.0 +  theta  + (i0)/ 40.0  + (count)/ 280.0),
                                0.0,
                                //0.5 + 0.5 * sin(5.0 +  theta + (i)/ 50.0  + (count)/ 90.0),
                                0.5 + 0.5 * sin(8.0 +  theta + (i0)/ 35.0  + (count)/ 300.0)
                            };
                            float greens[] =
                            {   0.5 + 0.5 * sin( 5.0 +  theta + (i0)/ 25.0 + (count)/ 180.0),
                                0.0,
                                //0.5 + 0.5 * sin(4.0 +  theta + (i)/ 25.0 + (count)/ 100.0),
                                0.5 + 0.5 * sin(6.0 +  theta + (i0)/ 35.0 + (count)/ 230.0)
                            };
                            float blues[] =
                            {   0.5 + 0.5 * sin(3.0 +  2.0 * theta + (i0)/ 30.0 + (count)/ 235.0),
                                0.0,
                                //0.5 + 0.5 * sin(7.0 +  2.0 * theta + (i)/ 20.0 + (count)/ 80.0),
                                0.5 + 0.5 * sin(0.0 +  2.0 * theta + (i0)/ 25.0 + (count)/ 205.0)
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
                            

                      
                            /*
                
                            if(i0 >= max_itr -1.0){
                                red = 0;
                                green = 0;
                                blue = 0;
                            }
                            */
                    }             
                  ///////////////


                      red = red * 255;
                      green = green * 255;
                      blue = blue * 255;

                          image[y * 4 * height2 + x * 4 + 0] = (int)(red); 
                          image[y * 4 * height2 + x * 4 + 1] = (int)(green); 
                          image[y * 4 * height2 + x * 4 + 2] = (int)(blue); 
                          image[y * 4 * height2 + x * 4 + 3] = 255; 
                  }
              }

              encodeOneStep((const char*)filename.data(), image, width2, height2);
          }

	  }

	    

	 
    return 0;
}
