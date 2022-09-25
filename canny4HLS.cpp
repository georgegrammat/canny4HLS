#include <stdio.h>
#include <iostream>
#include <fstream>
#include "ac_int.h"
#include "ac_channel.h"
#include "ac_fixed.h" 
#include "mc_scverify.h"

#include <ac_math/ac_atan2_cordic.h>
#include <ac_math/ac_sqrt.h>

typedef ac_int<11,false> dtype;
typedef ac_fixed<12,9,true> fixed_12_9; // dtype for theta angle
typedef ac_fixed<10,8,false> fixed; //0-255 
typedef ac_fixed<21,21,false> fixed_21_21; // for sum = Yx^2 + Yy^2 with max possible val 2080800
typedef ac_fixed<12,11,true> fixed_12_11;
typedef ac_fixed<11,11,false> fixed_11_11; 

static const fixed pi = 3.145;
static const int M = 366;
static const int N = 574; 

#pragma hls_design
void noiseReduction(ac_channel<dtype> &in, ac_channel<dtype> &out)
{   
    dtype localWindow[3][3], lineBuffer1[M+2], lineBuffer2[M+2]; //pixel cache variables
    dtype calc0,calc1,calc2,calc3,calc4,calc5,calc6,calc7,calc8, result;


    READ_CH:
    for (int j = 1; j < M+1; j++)
    {
        lineBuffer1[j] = 0;
        #ifndef __SYNTHESIS__
        if(in.available(1))
        #endif
        {
            lineBuffer2[j] = in.read(); 
        }
    }

    ROWS:
    for (int i = 1; i < N+1; i++)
    {
        COLS:
        for (int j = 0; j < M+1; j++)
        {    
            if(j==0)
            {   
                #ifndef __SYNTHESIS__
                if(in.available(1))
                #endif
                {
                localWindow[2][1] = in.read();
                }
            }
            else
            {
                    
				if(j==1){
					localWindow[0][0] = 0;
					localWindow[0][1] = lineBuffer1[j];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = 0;
					localWindow[1][1] = lineBuffer2[j];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = 0;
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j] = lineBuffer2[j];
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j] = localWindow[2][1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				else if(j==M)
				{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = 0;
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = 0;
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					localWindow[2][2] = 0;
				}
				else{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
			 
			 //gaussian blur
			 calc0 = localWindow[0][0]/(dtype)16;
			 calc1 = localWindow[0][1]/(dtype)8; //*0125;
			 calc2 = localWindow[0][2]/(dtype)16; //*0.0625
			 calc3 = localWindow[1][0]/(dtype)8; //*0.125
			 calc4 = localWindow[1][1]/(dtype)4; //*0.25;
			 calc5 = localWindow[1][2]/(dtype)8; //*0.125;
			 calc6 = localWindow[2][0]/(dtype)16; //*0.0625;
			 calc7 = localWindow[2][1]/(dtype)8; //*0.125;
			 calc8 = localWindow[2][2]/(dtype)16; //*0.0625;

			result = calc0+calc1+calc2+calc3+calc4+calc5+calc6+calc7+calc8;

			out.write(result);
			}
        }
    }
}

#pragma hls_design
void gradient_th_magn(ac_channel<dtype> &in, ac_channel<dtype> &out, ac_channel<fixed_12_9> &out2)
{
    fixed_12_9 theta_angle, temp;
    fixed rad2deg= 57.295; // 180 / PI
    fixed localWindow[3][3], lineBuffer1[M+2], lineBuffer2[M+2];
    fixed_11_11 sqrt_output;
    fixed_21_21 sum;
    dtype Gradient[N][M];
    fixed_12_11 Yx[N][M], Yy[N][M]; // convolution with sobel operators - possible max val -> 1020 -> 10bits
    
    READ_CH:
    for (int j = 1; j < M+1; j++)
    {
        lineBuffer1[j] = 0;
        #ifndef __SYNTHESIS__
        if(in.available(1))
        #endif
        {
            lineBuffer2[j] = in.read(); 
        }
    }
    
    ROWS:
    for (int i = 1; i < N+1; i++)
    {
        COLS:
        for (int j = 0; j < M + 1; j++)
        {    
            if(j==0)
            {   
                #ifndef __SYNTHESIS__
                if(in.available(1))
                #endif
                {
                localWindow[2][1] = in.read();
                }
            }
            else
            {    
				if(j==1){
					localWindow[0][0] = 0;
					localWindow[0][1] = lineBuffer1[j];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = 0;
					localWindow[1][1] = lineBuffer2[j];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = 0;
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j] = lineBuffer2[j];
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j] = localWindow[2][1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				else if(j==M)
				{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = 0;
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = 0;
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					localWindow[2][2] = 0;
				}
				else{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				
				//convolution - sobel operators
				Yy[i-1][j-1] = localWindow[0][0] + localWindow[0][1] * 2 + localWindow[0][2]  
								- localWindow[2][0]  - localWindow[2][1]*2 - localWindow[2][2];
				Yx[i-1][j-1] = -localWindow[0][0]  + localWindow[0][2] - localWindow[1][0]*2 
								+ localWindow[1][2]*2  - localWindow[2][0] + localWindow[2][2];

				ac_math::ac_atan2_cordic(Yy[i-1][j-1], Yx[i-1][j-1], temp);
				
				theta_angle = temp * rad2deg; // conversion from rad to degrees
				if(theta_angle < 0)
					theta_angle += 180;//0-180
				
				
				sum = Yx[i-1][j-1]*Yx[i-1][j-1] + Yy[i-1][j-1]*Yy[i-1][j-1];
				ac_math::ac_sqrt(sum, sqrt_output);
		   
				Gradient[i-1][j-1] = sqrt_output.to_int();
		   
				out.write(Gradient[i-1][j-1]); 
				out2.write(theta_angle);
				
			}
        }  
    }
}

#pragma hls_design
void nonMaxSuppresion(ac_channel<dtype> &in, ac_channel<dtype> &out, ac_channel<fixed_12_9> &in2)
{
    dtype localWindow[3][3], lineBuffer1[M+2], lineBuffer2[M+2], local; 
    fixed th[N][M];
   
    READ_CH:
    for (int j = 1; j < M+1; j++)
    {
        lineBuffer1[j] = 0;
        #ifndef __SYNTHESIS__
        if(in.available(1))
        #endif
        {
            lineBuffer2[j] = in.read();
            th[0][j-1] = in2.read();
        }
    }

    ROWS:
    for (int i = 1; i < N+1; i++)
    {   
        COLS:
        for (int j = 0; j < M+1; j++)
        {
            if(j==0)
            {   
                #ifndef __SYNTHESIS__
                if(in.available(1))
                #endif
                {
                localWindow[2][1] = in.read();
                th[i][0] = in2.read();
                }
            }
            else
            {    
			     if(j==1){
					localWindow[0][0] = 0;
					localWindow[0][1] = lineBuffer1[j];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = 0;
					localWindow[1][1] = lineBuffer2[j];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = 0;
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
					    
						localWindow[2][2] = in.read();
						th[i][j] = in2.read();
					}
					
					lineBuffer1[j] = lineBuffer2[j];
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j] = localWindow[2][1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				else if(j==M)
				{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = 0;
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = 0;
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					localWindow[2][2] = 0;
				}
				else{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					#ifndef __SYNTHESIS__
					if(in.available(2))
					#endif
					{
						localWindow[2][2] = in.read();
						th[i][j] = in2.read();
					} 
					
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				  
            dtype q = 255;
            dtype r = 255;

            if ((th[i-1][j-1] < 22.5 && th[i-1][j-1] >= 0) || (th[i-1][j-1] <= 180 && th[i-1][j-1] >= 157.5))
            {
                q = localWindow[1][2];
                r = localWindow[1][0];
            }
            else if ((th[i-1][j-1] < 67.5 && th[i-1][j-1] >= 22.5))
            {
                q = localWindow[2][0];
                r = localWindow[0][2];
            }
            else if ((th[i-1][j-1] < 112.5 && th[i-1][j-1] >= 67.5))
            {
                q = localWindow[2][1];
                r = localWindow[0][1];
            }
            else if ((th[i-1][j-1] < 157.5 && th[i-1][j-1] >= 112.5))
            {
                q = localWindow[0][0];
                r = localWindow[2][2];          
            }

            if (localWindow[1][1] >= q && localWindow[1][1] >= r)
            {
                local = localWindow[1][1];
            }
            else{
                local = 0;              
            }
            
            out.write(local);
         }
        } 
    }
}

#pragma hls_design
void threshold(ac_channel<dtype> &in, ac_channel<dtype> &out)
{
    dtype highThresh  = 80;
    dtype lowThresh = 30;
    dtype lineBuffer1[M+2], lineBuffer2[M+2], final, localWindow[3][3];

    READ_CH:
    for (int j = 1; j < M+1; j++)
    {
        lineBuffer1[j] = 0;
        #ifndef __SYNTHESIS__
        if(in.available(1))
        #endif
        {
            lineBuffer2[j] = in.read(); 
        }
    }
    
    ROWS:
    for (int i = 1; i < N+1; i++)
    {   
        COLS:
        for (int j = 0; j < M+1; j++)
        {   
            if(j==0)
            {   
                #ifndef __SYNTHESIS__
                if(in.available(1))
                #endif
                {
                localWindow[2][1] = in.read();
                }
            }
            else
            {            
				if(j==1){
					localWindow[0][0] = 0;
					localWindow[0][1] = lineBuffer1[j];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = 0;
					localWindow[1][1] = lineBuffer2[j];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = 0;
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j] = lineBuffer2[j];
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j] = localWindow[2][1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				else if(j==M)
				{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = 0;
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = 0;
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					localWindow[2][2] = 0;
				}
				else{
					localWindow[0][0] = localWindow[0][1];
					localWindow[0][1] = localWindow[0][2];
					localWindow[0][2] = lineBuffer1[j+1];
					localWindow[1][0] = localWindow[1][1];
					localWindow[1][1] = localWindow[1][2];
					localWindow[1][2] = lineBuffer2[j+1];
					localWindow[2][0] = localWindow[2][1];
					localWindow[2][1] = localWindow[2][2];
					#ifndef __SYNTHESIS__
					if(in.available(1))
					#endif
					{
						localWindow[2][2] = in.read();
					}
					lineBuffer1[j+1] = lineBuffer2[j+1];
					lineBuffer2[j+1] = localWindow[2][2];
				}
				if( localWindow[1][1] > highThresh){
				    final = 255;
				}
				
				else if (localWindow[1][1]< highThresh && localWindow[1][1] > lowThresh && ((localWindow[1][0] > highThresh) || 
					(localWindow[1][2] > highThresh) || (localWindow[0][1] > highThresh) || (localWindow[2][1]> highThresh) ||
					(localWindow[0][0] > highThresh) || (localWindow[2][0] > highThresh) || (localWindow[0][2] > highThresh) || 
					(localWindow[2][2] > highThresh)))
				{
					final = 255;
				}
				else
					final = 0;
                
				out.write(final);//0 or 255
			}         
        }

    }
}

#pragma hls_design top
void CCS_BLOCK(top)(ac_channel<dtype> &in, ac_channel<dtype> &out)
{
  static ac_channel<dtype> noiseOut;
  static ac_channel<dtype> gradientOut;
  static ac_channel<dtype> nonMaxOut;
  static ac_channel<fixed_12_9> theta_channel;
  
  noiseReduction(in, noiseOut);
  gradient_th_magn(noiseOut, gradientOut, theta_channel);   
  nonMaxSuppresion(gradientOut, nonMaxOut, theta_channel);
  threshold(nonMaxOut, out);
}

CCS_MAIN(int argc, char* argv[])
{   
	int gray[N][M];
	FILE *pFile;
	pFile = fopen("emma.yuv", "rb");
	if(pFile == NULL){
	    perror("error while opening file");
	}
	else{
	    for(int i=0; i<N; i++){
            for(int j=0; j<M; j++){
                gray[i][j] = fgetc(pFile);
            }
        }
	   
	    fclose(pFile);
    }
    
    static ac_channel<dtype> input1, output1;
    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            input1.write(gray[i][j]);
        }
    }

    top(input1, output1);
    
    FILE *frame_yuv;
    frame_yuv = fopen("finalResult.yuv", "wb");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            #ifndef __SYNTHESIS__
            if(output1.available(1))
			#endif
			{
                fputc(output1.read(), frame_yuv);
            }
        }
    }
    
    fclose(frame_yuv);   

    CCS_RETURN(0);
}
