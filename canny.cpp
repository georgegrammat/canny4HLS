#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>

static const double pi = 3.14;

#define filename "emma.bmp"

int **ReadBMP(int &height, int &width)
{
    //int i;
    FILE *f = fopen(filename, "rb");

    if (f == NULL)
        throw "Argument Exception";

    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    width = *(int *)&info[18];
    height = *(int *)&info[22];

    // cout << "  Name: " << filename << endl;
    //cout << " Width: " << width << endl;
    // cout << "Height: " << height << endl;

    int row_padded = (width * 3 + 3) & (~3);
    unsigned char *data = new unsigned char[row_padded];

    int **gray = 0;
    gray = new int *[height];

    int max = -1;
    for (int i = 0; i < height; i++)
    {
        gray[i] = new int[width];

        fread(data, sizeof(unsigned char), row_padded, f);
        for (int j = 0; j < width * 3; j += 3)
        {
            // Convert (B, G, R) to (R, G, B)
            int tmp = data[j];
            data[j] = data[j + 2];
            data[j + 2] = tmp;

            //cout << "R: " << (int)data[j] << " G: " << (int)data[j + 1] << " B: " << (int)data[j + 2] << endl;
            gray[i][j / 3] = ((int)data[j] + (int)data[j + 1] + (int)data[j + 2]) / 3;
            //std::cout << gray[i][j / 3] << " ";
            if (gray[i][j / 3] > max)
                max = gray[i][j / 3];
        }
    }

    fclose(f);

    return gray;
}

int **noiseReduction(int **gray, int N, int M)
{
    int **gray2 = 0;
    gray2 = new int *[N];

    for (int i = 0; i < N; i++)
    {
        gray2[i] = new int[M];
        for (int j = 0; j < M; j++)
        {
            if ((i == 0) || (i == N - 1) || (j == 0) || (j == M - 1))
            {
                gray2[i][j] = gray[i][j];
            }
            else
            {
                //gausian blur
                gray2[i][j] = gray[i - 1][j - 1] * 0.0625 + gray[i - 1][j] * 0.125 + gray[i - 1][j + 1] * 0.0625 + gray[i][j - 1] * 0.125 + gray[i][j] * 0.25 +
                              gray[i][j + 1] * 0.125 + gray[i + 1][j - 1] * 0.0625 +
                              gray[i + 1][j] * 0.125 + gray[i + 1][j + 1] * 0.0625;
            }
        }
    }

    return gray2;
}

int **sobelX(int **gray, int N, int M)
{

    int **Yx = 0;
    Yx = new int *[N];

    for (int i = 0; i < N; i++)
    {
        Yx[i] = new int[M];
        for (int j = 0; j < M; j++)
        {
            if ((i == 0) || (i == N - 1) || (j == 0) || (j == M - 1))
            {
                Yx[i][j] = 0;
            }
            else
            {
                //sobel Yx
                Yx[i][j] = gray[i - 1][j - 1] * 1 + gray[i - 1][j] * 0 + gray[i - 1][j + 1] * -1 + gray[i][j - 1] * 2 + gray[i][j] * 0 + gray[i][j + 1] * -2 + gray[i + 1][j - 1] * 1 + gray[i + 1][j] * 0 + gray[i + 1][j + 1] * -1;
            }
        }
    }

    return Yx;
}

int **sobelY(int **gray, int N, int M)
{

    int **Yy = 0;
    Yy = new int *[N];

    for (int i = 0; i < N; i++)
    {
        Yy[i] = new int[M];
        for (int j = 0; j < M; j++)
        {
            if ((i == 0) || (i == N - 1) || (j == 0) || (j == M - 1))
            {
                Yy[i][j] = 0;
            }
            else
            {
                //sobel Yy
                Yy[i][j] = gray[i - 1][j - 1] * 1 + gray[i - 1][j] * 2 + gray[i - 1][j + 1] * 1 + gray[i][j - 1] * 0 + gray[i][j] * 0 + gray[i][j + 1] * 0 + gray[i + 1][j - 1] * -1 + gray[i + 1][j] * -2 + gray[i + 1][j + 1] * -1;
                //std::cout << Yy[i][j] << " ";
            }
        }
    }

    return Yy;
}

float **theta(int **Yx, int **Yy, int N, int M)
{

    float **th = 0;
    th = new float *[N];

    for (int i = 0; i < N; i++)
    {
        th[i] = new float[M];

        for (int j = 0; j < M; j++)
        {
            th[i][j] = ((atan2(Yx[i][j], Yy[i][j])) * 180 / pi);
            if (th[i][j] < 0)
                th[i][j] += 180;
        }
    }

    return th;
}

int **gradient(int **Yx, int **Yy, int N, int M)
{

    int **G = 0;
    G = new int *[N];

    for (int i = 0; i < N; i++)
    {
        G[i] = new int[M];

        for (int j = 0; j < M; j++)
        {
            G[i][j] = sqrt(pow(Yx[i][j], 2) + pow(Yy[i][j], 2)) / 2;
            //std::cout << G[i][j] << " ";
        }
    }
    return G;
}

int **nonMaxSuppresion(float **th, int **G, int N, int M)
{
    int **Z = 0;
    Z = new int *[N];
    int max = -1;
    for (int i = 0; i < N; i++)
    {
        Z[i] = new int[M];

        for (int j = 0; j < M; j++)
        {
            if ((i == 0) || (i == N - 1) || (j == 0) || (j == M - 1))
            {
                Z[i][j] = G[i][j];
            }
            else
            {
                int q = 255;
                int r = 255;

                if ((th[i][j] < 22.5 && th[i][j] >= 0) || (th[i][j] <= 180 && th[i][j] >= 157.5))
                {
                    q = G[i][j + 1];
                    r = G[i][j - 1];
                }
                else if ((th[i][j] < 67.5 && th[i][j] >= 22.5))
                {
                    q = G[i + 1][j - 1];
                    r = G[i - 1][j + 1];
                }
                else if ((th[i][j] < 112.5 && th[i][j] >= 67.5))
                {
                    q = G[i + 1][j];
                    r = G[i - 1][j];
                }
                else if ((th[i][j] < 157.5 && th[i][j] >= 112.5))
                {
                    q = G[i - 1][j - 1];
                    r = G[i + 1][j + 1];
                }

                if (G[i][j] >= q && G[i][j] >= r)
                {
                    Z[i][j] = G[i][j];
                    if (Z[i][j] > max)
                        max = Z[i][j];
                }
                else
                {
                    Z[i][j] = 0;
                    //std::cout << Z[i][j] << " ";
                }
            }
        }
    }

    return Z;
}

int **threshold(int **Z, int N, int M)
{
    int **final = 0;
    final = new int *[N];

    double highThresh = 55;
    double lowThresh = 18;

    for (int i = 0; i < N; i++)
    {
        final[i] = new int[M];
        for (int j = 0; j < M; j++)
        {
            if ((i == 0) || (i == N - 1) || (j == 0) || (j == M - 1))
            {
                final[i][j] = Z[i][j];
            }
            else
            {
                if (Z[i][j] >= highThresh)
                {
                    final[i][j] = 255;
                }

                else if (Z[i][j] <= highThresh && Z[i][j] >= lowThresh)
                {

                    if ((Z[i][j - 1] > highThresh) || (Z[i][j + 1] > highThresh) || (Z[i - 1][j] > highThresh) || (Z[i + 1][j] > highThresh) ||
                        (Z[i - 1][j - 1] > highThresh) || (Z[i + 1][j - 1] > highThresh) || (Z[i - 1][j + 1] > highThresh) || (Z[i + 1][j + 1] > highThresh))
                        final[i][j] = Z[i][j];
                }
                else
                    final[i][j] = 0;
            }
        }
    }

    return final;
}

void write(int **Z, int N, int M)
{
    FILE *frame_yuv;
    frame_yuv = fopen("result.yuv", "wb");

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = 0; j < M; j++)
        {

            fputc(Z[i][j], frame_yuv);
        }
    }

    fclose(frame_yuv);
}

int main()
{
    int height, width;

    int **gray = ReadBMP(height, width);

    std::cout << "height: " << height << "\n";
    std::cout << "width: " << width << "\n";

    int **noiseFree = noiseReduction(gray, height, width);

    int **xOperator = sobelX(noiseFree, height, width);
    int **yOperator = sobelY(noiseFree, height, width);

    int **Gradient = gradient(xOperator, yOperator, height, width);

    float **th = theta(xOperator, yOperator, height, width);

    int **nonMax = nonMaxSuppresion(th, Gradient, height, width);

    int **threshImg = threshold(nonMax, height, width);

    write(threshImg, height, width);

    return 0;
}
