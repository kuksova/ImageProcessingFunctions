//FingerprintSegmentation
//AOI 

#include "stdafx.h"
#include "FingerprintSegmentation.h"
#include <string.h>
#include <cstdlib>
#include <vector>
#include "queue.h"

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

typedef struct Position{
	unsigned int row;
	unsigned int col;
} Position;

//==================================================================================================
//
//==================================================================================================
unsigned int  EstimateMaxMin(BYTE * ImgAddress, ProcParam Param, unsigned int X_OFFSET, unsigned int Y_OFFSET, unsigned int x_size, unsigned int y_size, unsigned char *vmin, unsigned char *vmax)
{
	int gist[256];
	unsigned int j, i;
	unsigned int sum = 0;
	unsigned int SumInt = 0;
	(*vmin) = 255;
	(*vmax) = 0;

	memset(gist, 0, sizeof(gist));

	for (j = 0; j < (y_size); j++){
		for (i = 0; i < (x_size); i++)
		{
			(*vmin) = min(ImgAddress[(j + Y_OFFSET)*Param.Width + i + X_OFFSET], (*vmin));
			(*vmax) = max(ImgAddress[(j + Y_OFFSET)*Param.Width + i + X_OFFSET], (*vmax));
			gist[ImgAddress[(j + Y_OFFSET)*Param.Width + i + X_OFFSET]]++;
		}
	}
	return 0;
}
//==================================================================================================
//
//==================================================================================================
void CalculateGistAndMin(BYTE *InBuf, ProcParam Param, int *pnMin /*= NULL*/, int *pnGist /*= NULL*/)
{
	int nSum;
	int* nGist = new int[256];
	for (int t = 0; t < 256; t++)
		nGist[t] = 0;

	for (int i = 0; i < Param.Width*Param.Height; i++)
	{
		nGist[InBuf[i]]++;
	}
	if (pnMin)
	{
		int nMin = 0; nSum = 0;
		while (nSum < 0.0065*Param.Width*Param.Height && nMin < 255) // IT WAS 1000
		{
			nSum += nGist[nMin];
			nMin++;
		}
		if (nMin > 64) nMin = 64;
		*pnMin = nMin;
		//
	}
	for (int t = 0; t < 256; t++)
		pnGist[t] = nGist[t];
}

int Norm(BYTE *InBuf, ProcParam Param, int nMin)
{
	int i, nSum;
	const int nMax = 255;

	if (nMin < 0)
	{
		CalculateGistAndMin(InBuf, Param, &nMin, NULL);
	}
	if (nMax <= nMin)
		return 0;
	for (i = 0; i < Param.Width; i++)
	{
		for (int j = 0; j < Param.Height; j++)
		{
			nSum = (InBuf[i*Param.Width + j] - nMin) * nMax;
			nSum /= (nMax - nMin);
			if (nSum > 255)
			{
				nSum = 255;
			}
			else {
				if (nSum < 0) nSum = 0;
			}
			InBuf[i*Param.Width + j] = (BYTE)nSum;
		}

	}
	return 1;
}
//==================================================================================================
//
//==================================================================================================
void flood_fill(BYTE *Img, ProcParam Param, Position from, BYTE color)
{
	Queue q; /*used to conveniently Queue*/
	initQueue(&q);

	Position * _from = (Position*)malloc(sizeof(Position));
	_from->row = from.row; // yc
	_from->col = from.col; // xc

	Position * possible = (Position*)malloc(sizeof(Position));

	enqueue(&q, _from);
	int X, Y;
	while (!isEmpty(&q)){
		//while (n<1000){
		Position c = *(Position *)dequeue(&q);
		Img[c.row*Param.Width + c.col] = color;

		X = c.col;
		Y = c.row;
		int X_min = X; // found the leftmost boundary point
		while (Img[Y*Param.Width + X_min - 1] != 255) // 255 is white border
		{
			X_min--;
		}
		int X_max = X; // found the rightmost boundary point
		while (Img[Y*Param.Width + X_max + 1] != 255) 
		{
			X_max++;
		}
		// filling row contains the started point
		for (int i = X_min; i <= X_max; i++)
		{
			Img[Y*Param.Width + i] = color;
		}
		//processing row from the top
		bool flag = true;
		for (X = X_min; X <= X_max; X++)
		{
			if ((Img[(Y - 1)*Param.Width + X] != 255) && (Img[(Y - 1)*Param.Width + X] != color))
			{
				if (flag == true)
				{
					possible->row = Y - 1;
					possible->col = X;
					enqueue(&q, possible);
					flag = false;
				}
			}
			else flag = true;
		}
	}
	enqueue(&q, _from);

	while (!isEmpty(&q)){
		Position c = *(Position *)dequeue(&q);
		Img[c.row*Param.Width + c.col] = color;

		X = c.col;
		Y = c.row;
		int X_min = X;
		while (Img[Y*Param.Width + X_min - 1] != 255) // 255 is white border
		{
			X_min--;
		}
		int X_max = X;
		while (Img[Y*Param.Width + X_max + 1] != 255) // 255 is white border
		{
			X_max++;
		}

		for (int i = X_min; i <= X_max; i++)
		{
			Img[Y*Param.Width + i] = color;
		}

		//processing row from the bottom
		bool flag = true;
		for (int X = X_min; X <= X_max; X++)
		{
			if ((Img[(Y + 1)*Param.Width + X] != 255) && (Img[(Y + 1)*Param.Width + X] != color))
			{
				if (flag == true)
				{
					possible->row = Y + 1;
					possible->col = X;
					enqueue(&q, possible);
					flag = false;
				}
			}
			else flag = true;
		}

	}

	free(possible);
}
//==================================================================================================
//
//==================================================================================================
void del_pixels(BYTE *InBuf, int row, int col, int size, color_histo_t *h, ProcParam Param)
{
	int i;
	//rgb_t *pix;

	if (col < 0 || col >= Param.Width) return;
	for (i = row - size; i <= row + size && i < Param.Height; i++) {
		if (i < 0) continue;
		//pix = im->pix[i] + col;
		//h->r[InBuf[i] + col]--;
		h->r[InBuf[i*Param.Width + col]]--;
		h->n--;
	}
}

void add_pixels(BYTE *InBuf, int row, int col, int size, color_histo_t *h, ProcParam Param)
{
	int i;
	//rgb_t *pix;

	if (col < 0 || col >= Param.Width) return;
	for (i = row - size; i <= row + size && i < Param.Height; i++) {
		if (i < 0) continue;
		//		pix = im->pix[i] + col;
		h->r[InBuf[i*Param.Width + col]]++;
		h->n++;
	}
}

void init_histo(BYTE *InBuf, int row, int size, color_histo_t*h, ProcParam Param)
{
	int j;

	memset(h, 0, sizeof(color_histo_t));

	for (j = 0; j < size && j < Param.Width; j++)
		add_pixels(InBuf, row, j, size, h, Param);
}

int median(const int *x, int n)
{
	int i;
	for (n /= 2, i = 0; i < 256 && (n -= x[i]) > 0; i++);
	return i;
}

BYTE median_color(const color_histo_t *h)
{
	unsigned char pix;
	return pix = median(h->r, h->n);
}

void MedianFilter(BYTE *InBuf, BYTE *OutBuf, ProcParam Param, int size)
{
	color_histo_t h;
	BYTE pix = 0;
	for (int row = 0; row < Param.Height; row++)
	{
		for (int col = 0; col < Param.Width; col++)
		{
			if (!col) init_histo(InBuf, row, size, &h, Param);
			else {
				del_pixels(InBuf, row, col /*- 1 */ - size, size, &h, Param);
				add_pixels(InBuf, row, col /*- 1*/ + size, size, &h, Param);
			}
			pix = median_color(&h);
			OutBuf[row*Param.Width + col] = (BYTE)pix;
		}
	}
}
//==================================================================================================
//
//==================================================================================================
void dilate_mask(BYTE *InBuf, BYTE *OutBuf, BYTE *mask, ProcParam Param, ProcParam MaskPar)
{
	int c, d, k, l;
	int ind;
	//	int M = MaskPar.Width + 1;
	//	int ind;
	//	int offset1[9] = { -M, -M + 1, -M + 2, -1, 0, +1, M - 2, M - 1, M };
	//	int offset[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	BYTE *Im;
	Im = (BYTE *)malloc((Param.Width /*+ 10*/)*(Param.Height /*+ 10*/)*sizeof(BYTE));
	memset(&Im[0], 0, Param.Width*Param.Height);
	for (int j = MaskPar.Height / 2; j < Param.Height - MaskPar.Height / 2; j++)
	{
		for (int i = MaskPar.Width / 2; i < Param.Width - MaskPar.Width / 2; i++)
		{
			if (InBuf[j*Param.Width + i] == 255)
			{
				for (k = -MaskPar.Height / 2/*+2*/; k <= MaskPar.Height / 2; k++)
				{
					for (l = -MaskPar.Width / 2/*+2*/; l <= MaskPar.Width / 2; l++)
						//			for (int m = 0; m <= 8; m++)
					{
						ind = (k + MaskPar.Height / 2)* MaskPar.Width + (l + MaskPar.Width / 2);
						if (mask[ind] == 255)
							//	if (mask[offset[m]] == 255)
						{
							c = j + k;
							d = i + l;
							//	printf("%5d %5d", c, d);
							//Im[j*Param.Width + i+offset1[m]] = 255;
							Im[c*Param.Width + d] = 255;
						}
					}
				}
			}

		}
	}
	memcpy(OutBuf, Im, Param.Width*Param.Height);
	free(Im);
}
//==================================================================================================
//
//==================================================================================================
void erode_mask(BYTE *InBuf, BYTE *OutBuf, BYTE *mask, ProcParam Param, ProcParam MaskPar)
{
	int c = 0, d;
	BYTE *Im;
	Im = (BYTE *)malloc((Param.Width /*+ 10*/)*(Param.Height /*+ 10*/)*sizeof(BYTE));
	memset(&Im[0], 0, Param.Width*Param.Height);

	int nm = 0;
	for (int k1 = 0; k1 < MaskPar.Height; k1++)
	{
		for (int l1 = 0; l1 < MaskPar.Width; l1++)
		{
			if (mask[k1*MaskPar.Width + l1] == 255)
				nm++;
		}
	}
	/*erosion squre element*/
	for (int j = 0; j < Param.Height - MaskPar.Height; j++)
	{
		for (int i = 0; i < Param.Width - MaskPar.Width; i++)
		{
			if (InBuf[j*Param.Width + i] == 255)
			{
				for (int k = 0; k < MaskPar.Height; k++)
				{
					for (int l = 0; l < MaskPar.Width; l++)
					{
						if (mask[k*MaskPar.Width + l] == InBuf[(j + k)*Param.Width + (i + l)])
						{
							c++;
							//c = j - k;
							//d = i - l;
							//Im[(j /*+ MaskPar.Height / 2)*Param.Width + (i /*+ MaskPar.Width / 2)] = 255;
						}
					}
				}

				if (c == nm /*(MaskPar.Height*MaskPar.Width)*/)
				{
					//printf("\n Im = %5d", (j + MaskPar.Height / 2)*Param.Width + (i + MaskPar.Width / 2));
					Im[(j + MaskPar.Height / 2)*Param.Width + (i + MaskPar.Width / 2)] = 255;
				}
				c = 0;
			}

		}
	}
	memcpy(OutBuf, Im, Param.Width*Param.Height);
	free(Im);
}
//==================================================================================================
//
//==================================================================================================
void  imerode_(BYTE *Img_d, BYTE *rect_pol, int H, ProcParam Param)
{
	int l;
	BYTE T;
	BYTE *M;
	int y, x, j, i;
	/*M = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	for (y = H; y<Param.Height - H; y++){
	l = Param.Width*y;
	for (x = 0; x<Param.Width; x++){
	T = rect_pol[x + (y - H)*Param.Width];
	for (j = y - H + 1; j <= y + H; j++){
	if (rect_pol[x + j*Param.Width]<T) T = rect_pol[x + j*Param.Width];
	}
	M[x + l] = T;
	}
	for (x = H; x<Param.Width - H; x++){
	T = M[x - H + l];
	for (i = x - H + 1; i <= x + H; i++){
	if (M[i + l]<T) T = (BYTE)M[i + l];
	}
	Img_d[x + y*Param.Width] = T;
	}
	}
	free(M);*/

	//int K[7/*3 * H + 1*/] = { 1, 2, 3, 3, 3, 2, 1 };
	//int K[11/*3 * H + 1*/] = { 1, 2, 3, 4, 5, 5, 5, 4, 3, 2, 1 };
	//int K[11/*3 * H + 1*/] = { 3, 4, 5, 5, 5, 5, 5, 5, 5, 4, 3 };
	//int K[21/*3 * H + 1*/] = { 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 8, 7, 6, 5 };
	//int K[11] = { 3, 4, 5, 5, 5, 5, 5, 5, 5, 4, 3 };
	int K[13] = { 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 5, 4, 3 };
	int k = 0;
	M = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	for (y = H; y<Param.Height - H; y++){
		l = Param.Width*y;
		for (x = 0; x<Param.Width; x++){
			T = rect_pol[x + (y - H)*Param.Width];
			//for (j = y - H + 1; j <= y + H; j++){
			for (j = y - K[k] + 1; j <= y + K[k]; j++){
				//k = k + 1;
				if (rect_pol[x + j*Param.Width]<T) T = rect_pol[x + j*Param.Width];
			}
			k = k + 1;
			if (k == (H + H + 1))
				k = 0;
			k = 0;
			M[x + l] = T;
		}
		for (x = H; x<Param.Width - H; x++){
			T = M[x - H + l];
			for (i = x - H + 1; i <= x + H; i++){
				if (M[i + l]<T) T = (BYTE)M[i + l];
			}
			Img_d[x + y*Param.Width] = T;
		}
	}
	free(M);
	// Processing boards H
	for (y = 0; y<H; y++){
		l = Param.Width*y;
		for (x = 0; x<Param.Width; x++)
		{
			T = rect_pol[x + y*Param.Width];
			for (i = max(0, x - H); i <= min(x + H, Param.Width - 1); i++)
				for (j = max(0, y - H); j <= min(y + H, Param.Height - 1); j++)
					if (rect_pol[i + j*Param.Width]<T) T = rect_pol[i + j*Param.Width];
			Img_d[x + y*Param.Width] = T;
		}
	}
	for (y = Param.Height - H; y<Param.Height; y++){
		l = Param.Width*y;
		for (x = 0; x<Param.Width; x++)
		{
			T = rect_pol[x + y*Param.Width];
			for (i = max(0, x - H); i <= min(x + H, Param.Width - 1); i++)
				for (j = max(0, y - H); j <= min(y + H, Param.Height - 1); j++)
					if (rect_pol[i + j*Param.Width]<T) T = rect_pol[i + j*Param.Width];
			Img_d[x + y*Param.Width] = T;
		}
	}
	for (x = 0; x<H; x++)
	{
		for (y = 0; y<Param.Height; y++)
		{
			T = rect_pol[x + y*Param.Width];
			for (i = max(0, x - H); i <= min(x + H, Param.Width - 1); i++)
			{
				for (j = max(0, y - H); j <= min(y + H, Param.Height - 1); j++)
				{
					if (rect_pol[i + j*Param.Width]<T) T = rect_pol[i + j*Param.Width];
				}
			}
			Img_d[x + y*Param.Width] = T;
		}
	}
	for (x = Param.Width - H; x<Param.Width; x++)
	{
		for (y = 0; y<Param.Height; y++)
		{
			T = rect_pol[x + y*Param.Width];
			for (i = max(0, x - H); i <= min(x + H, Param.Width - 1); i++)
			{
				for (j = max(0, y - H); j <= min(y + H, Param.Height - 1); j++)
				{
					if (rect_pol[i + j*Param.Width]<T) T = rect_pol[i + j*Param.Width];
				}
			}
			Img_d[x + y*Param.Width] = T;
		}
	}
}
//==================================================================================================
//
//==================================================================================================
void imclose_(BYTE *rect_pol, BYTE *rect_e, int H, ProcParam Param)
{
	BYTE *rect;
	rect = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect1;
	rect1 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect2;
	rect2 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect3;
	rect3 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));

	ProcParam MPar1;
	MPar1.Height = 13;
	MPar1.Width = 13;
	static BYTE Mask1[] = {
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, };
	dilate_mask(rect_pol, rect, Mask1, Param, MPar1);
	dilate_mask(rect, rect1, Mask1, Param, MPar1);
	//	imdilate_(rect, rect_pol,H, Param);
	//imerode_(rect1, rect, H/*8*/, Param);
	imerode_(rect2, rect1, H/*8*/, Param); // djn
	//	imerode_(rect_e, rect2, H/*8*/, Param);
	erode_mask(rect1, rect_e, Mask1, Param, MPar1);
	//	erode_mask(rect1, rect_e, Mask1, Param, MPar1);
	//	imerode_(rect_e, rect1, 6/*8*/, Param);
	free(rect);
	free(rect1);
	free(rect2);
	free(rect3);
}
//==================================================================================================
//
//==================================================================================================
void imopen_(BYTE *rect_pol, BYTE *rect_e, int H, ProcParam Param)
{
	BYTE *rect;
	rect = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect1;
	rect1 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect2;
	rect2 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect3;
	rect3 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect4;
	rect4 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *rect5;
	rect5 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	imerode_(rect, rect_pol, H, Param);
	//	imerode_(rect1, rect, H, Param);
	//	imerode_(rect2, rect1, H, Param);

	//imdilate_(rect_e, rect1,H, Param);

	ProcParam MPar1;
	MPar1.Height = 13;
	MPar1.Width = 13;
	static BYTE Mask1[] = {
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0 };

	ProcParam MPar2;
	MPar2.Height = 7;
	MPar2.Width = 7;
	static BYTE Mask2[] = {
		0, 0, 255, 255, 255, 0, 0,
		0, 255, 255, 255, 255, 255, 0,
		255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255,
		0, 255, 255, 255, 255, 255, 0,
		0, 0, 255, 255, 255, 0, 0
	};

	//	erode_mask(rect_pol, rect, Mask1, Param, MPar1);
	erode_mask(rect, rect1, Mask1, Param, MPar1);
	//	imdilate_(rect2, rect1, H, Param);
	//	imdilate_(rect_e, rect2, H, Param);


	dilate_mask(rect1, rect3, Mask1, Param, MPar1);
	dilate_mask(rect3, rect_e, Mask1, Param, MPar1);

	//	dilate_mask(rect, rect3, Mask1, Param, MPar1);
	//	dilate_mask(rect3, rect4, Mask1, Param, MPar1);
	//	dilate_mask(rect4, rect_e, Mask2, Param, MPar2);
	free(rect);
	free(rect1);
	free(rect2);
	free(rect3);
	free(rect4);
	free(rect5);
}
//==================================================================================================
//
//==================================================================================================
void MorfOperations(BYTE *Img, BYTE *ImgD, int H, ProcParam Param)
{
	BYTE *temp;
	temp = (BYTE *)malloc(Param.Height*Param.Width*sizeof(BYTE));

	ProcParam MPar1;
	MPar1.Height = 13;
	MPar1.Width = 13;
	static BYTE Mask1[] = {
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
		0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
		0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0,
		0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, };

	//	dilate_mask(Img, temp, Mask1, Param, MPar1);
	//	erode_mask(temp, ImgD, Mask1, Param, MPar1);
	// Apllying Close various masks
	imclose_(Img, temp, H, Param);
	//	imdilate_(ImgD, Img, 6, Param);

	// Apllying Open various masks
	imopen_(temp, ImgD, H, Param);

	free(temp);
	return;
}
//==================================================================================================
//
//==================================================================================================
int Bwlabel(BYTE * mas, BYTE * Label, ProcParam Param)
{
	int Num = 0;
	std::vector<int> pos, pos1, ind, v_x, v_y;
	BYTE *mas1;
	mas1 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&mas1[0], 0, Param.Width*Param.Height);
	BYTE *mas2;
	mas2 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&mas2[0], 0, Param.Width*Param.Height);
	BYTE *mas3;
	mas3 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&mas3[0], 0, Param.Width*Param.Height);
	BYTE *Y;
	Y = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Y[0], 0, Param.Width*Param.Height);

	ProcParam MPar1;
	MPar1.Height = 3;
	MPar1.Width = 3;
	static BYTE Mask1[] = {
		255, 255, 255,
		255, 255, 255,
		255, 255, 255,
	};


	for (int i = 0; i < Param.Height; i++){
		for (int j = 0; j < Param.Width; j++){
			if (mas[i*Param.Width + j] == 255)
			{
				pos.push_back(i); //y
				pos1.push_back(j); //x
			}
		}
	}
	int x, y;
	int N1 = 0;
	bool equal = true;
	while ((!pos.empty()) && (!pos1.empty()))
	{
		Num++;
		x = pos1[0];
		y = pos[0];

		memset(&mas1[0], 0, Param.Height*Param.Width);
		mas1[y * Param.Width + x] = 255;
		dilate_mask(mas1, mas2, Mask1, Param, MPar1);
		//Dilation(mas2, mas1, 1, Param);

		for (int i = 0; i < Param.Height; i++)
		{
			for (int j = 0; j < Param.Width; j++)
			{

				if ((mas[i*Param.Width + j] == 255) &(mas2[i*Param.Width + j]) == 255)
				{
					Y[i*Param.Width + j] = 255;
				}
				else Y[i*Param.Width + j] = 0;
				if (mas1[i*Param.Width + j] == Y[i*Param.Width + j])
					N1++;
			}
		}

		if (N1 != Param.Width*Param.Height)
			equal = false;

		while (!equal)
		{
			memcpy(mas1, Y, Param.Width*Param.Height);
			dilate_mask(mas1, mas3, Mask1, Param, MPar1);
			//	Dilation(mas3, mas1, 1, Param);

			for (int i = 0; i < Param.Height; i++)
			{
				for (int j = 0; j < Param.Width; j++)
				{
					if ((mas[i*Param.Width + j] == 255) &(mas3[i*Param.Width + j]) == 255)
					{
						Y[i*Param.Width + j] = 255;
					}
				}
			}

			/////////////////////////////
			N1 = 0;
			for (int i = 0; i < Param.Height; i++)
			{
				for (int j = 0; j < Param.Width; j++)
				{
					if (mas1[i*Param.Width + j] == Y[i*Param.Width + j])
						N1++;
				}
			}
			if (N1 != Param.Width*Param.Height)
				equal = false;
			else
			{
				equal = true;
			}
		}

		for (int i = 0; i < Param.Height; i++){
			for (int j = 0; j < Param.Width; j++){
				if (Y[i*Param.Width + j] == 255)
				{
					v_y.push_back(i); //y
					v_x.push_back(j); //x

				}
			}
		}

		for (int i = 0; i < v_x.size(); i++){
			mas[v_y[i] * Param.Width + v_x[i]] = 0;
			Label[v_y[i] * Param.Width + v_x[i]] = Num;
		}

		pos.clear();
		pos1.clear();
		v_x.clear();
		v_y.clear();
		for (int i = 0; i < Param.Height; i++){
			for (int j = 0; j < Param.Width; j++){
				if (mas[i*Param.Width + j] == 255)
				{
					pos.push_back(i); //y
					pos1.push_back(j); //x
					//ind.push_back(i*Param.Width + j);
				}
			}
		}
	}
	free(mas1);
	free(mas2);
	free(mas3);
	free(Y);
	return Num;
}
//==================================================================================================
//
//==================================================================================================
void ind2sub(const int sub, const int cols, const int rows, int &row, int &col)
{
	row = sub / cols;
	col = sub%cols;
}
void Boundary(BYTE *Img, ProcParam Param, std::vector<int> &B, std::vector<int>&B1, std::vector<int> &scratch)
{
	int M = Param.Width;
	//int offset[8] = { -1, M - 1, M, M + 1, 1, -M + 1, -M, -M - 1 };
	int offset[8] = { -M, -M + 1, 1, M + 1, M, M - 1, -1, -M - 1 };
	int next_search_direction_lut[8] = { 8, 8, 2, 2, 4, 4, 6, 6 };
	int next_direction_lut[8] = { 2, 3, 4, 5, 6, 7, 8, 1 };
	int START = 254;//-1;
	int BOUNDARY = 253;//-2;

	std::vector<int> rr;
	std::vector<int> cc;
	int k2 = 0, k3 = 0;
	int idx, which, numPixels = 0, currentPixel, done, next_search_direction, direction, found_next_pixel;
	int neighbor;
	int row = 0;
	int col = 0;
	BYTE *Lp1;
	Lp1 = (BYTE *)malloc((Param.Width)*(Param.Height + 2)*sizeof(BYTE));
	memset(&Lp1[0], 0, (Param.Width)*(Param.Height + 2));
	//memcpy(Lp1, mas, Param.Width*Param.Height);

	BYTE *Lp2;
	Lp2 = (BYTE *)malloc((Param.Width)*(Param.Height)*sizeof(BYTE));
	BYTE *Lp3;
	Lp3 = (BYTE *)malloc((Param.Width)*(Param.Height)*sizeof(BYTE));
	memset(&Lp2[0], 0, (Param.Width)*(Param.Height));
	BYTE *Lp2_1;
	Lp2_1 = (BYTE *)malloc((Param.Width)*(Param.Height)*sizeof(BYTE));
	BYTE *Lp3_1;
	Lp3_1 = (BYTE *)malloc((Param.Width)*(Param.Height)*sizeof(BYTE));
	memset(&Lp2_1[0], 0, (Param.Width)*(Param.Height));
	BYTE *Out;
	Out = (BYTE *)malloc((Param.Width)*(Param.Height)*sizeof(BYTE));
	memset(&Out[0], 0, (Param.Width)*(Param.Height));

	////////////////////////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j < (Param.Height) - 1; j++)
	{
		for (int i = 0; i < (Param.Width); i++)
		{
			/*if (Lp1[(j+1)*Param.Width + (i)] >0)
			Lp2[(j)*Param.Width + (i)] = Lp1[(j+1)*Param.Width + (i)];*/
			Lp1[(j + 1)* Param.Width + i] = Img[(j)* Param.Width + i];
			//printf("%5d ", Lp1[(j) * Param.Width + i]);
		}
	}
	/////////////////////////////////////////////////////// это соотвествует А1
	for (int j = 0; j < (Param.Height); j++)
	{
		for (int i = 0; i < (Param.Width); i++)
		{
			Lp2[(j)* Param.Width + i] = Lp1[(j + 1)* Param.Width + i];
		}
	}
	for (int j = 0; j < (Param.Height); j++)
	{
		for (int i = 0; i < (Param.Width); i++)
		{
			if (Lp2[(j)* Param.Width + i]>0)
				Lp2_1[(j)* Param.Width + i] = Lp2[(j)* Param.Width + i];
		}
	}
	////////////////////////////////////////////////// это соответсвует А2
	for (int j = 0; j < (Param.Height); j++)
	{
		for (int i = 0; i < (Param.Width); i++)
		{
			Lp3[(j)* Param.Width + i] = Lp1[(j)* Param.Width + i];
		}
	}

	for (int j = 0; j < (Param.Height); j++)
	{
		for (int i = 0; i < (Param.Width); i++)
		{
			if (Lp3[(j)* Param.Width + i] == 0)
				Lp3_1[(j)* Param.Width + i] = 255;
			else
			{
				Lp3_1[(j)* Param.Width + i] = 0;
			}
			/*printf("%5d ", */Lp3_1[(j)* Param.Width + i]/*)*/;
		}
	}
	//////////////////////////////////////////////////
	/*Lp2_1 & Lp3_1*/
	for (int j = 0; j < Param.Height; j++)
	{
		for (int i = 0; i < Param.Width; i++)
		{

			if (Lp2_1[(j)*Param.Width + (i)] == Lp3_1[(j)*Param.Width + (i)])
				Out[(j)*Param.Width + (i)] = (BYTE)255;
			if ((Lp2_1[(j)*Param.Width + (i)] == (BYTE)255) && (Lp3_1[(j)*Param.Width + (i)] == (BYTE)0))
				Out[(j)*Param.Width + (i)] = (BYTE)0;
			if ((Lp2_1[(j)*Param.Width + (i)] == (BYTE)0) && (Lp3_1[(j)*Param.Width + (i)] == (BYTE)255))
				Out[(j)*Param.Width + (i)] = (BYTE)0;
			if ((Lp2_1[(j)*Param.Width + (i)] == (BYTE)0) && (Lp3_1[(j)*Param.Width + (i)] == (BYTE)0))
				Out[(j)*Param.Width + (i)] = (BYTE)0;


		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j < Param.Height; j++)
	{
		for (int i = 0; i < Param.Width; i++)
		{
			if (Out[j*Param.Width + i] == 255){
				rr.push_back(i);
				cc.push_back(j);
			}
			Out[j*Param.Width + i];
		}
	}
	k2 = rr.size();
	k3 = cc.size();

	for (int i = 0; i < cc.size(); i++)
	{
		cc[i] = cc[i] + 1;
	}

	int m1, m2, m3;
	int initial_departure_direction = 0;
	for (int k = 0; k < rr.size(); k++)
	{
		m1 = cc[k];
		m2 = rr[k];
		if ((Lp1[cc[k] * Param.Width + rr[k]] > 0) == (Lp1[(cc[k] - 1) * Param.Width + (rr[k])] == 0 == B.empty()))
		{
			idx = (cc[k])*Param.Width + rr[k]; // first point of contour
			which = Lp1[idx];
			scratch.push_back(idx); // put index point in scracth
			Lp1[idx] = (BYTE)START; // Noted as visited
			numPixels = 1;
			currentPixel = idx;
			initial_departure_direction = 0/*[]*/;
			done = 0;
			next_search_direction = 2; // Further looks okrestnosti points in certain areas
			while (!done){
				direction = next_search_direction;
				found_next_pixel = 0;
				for (int n5 = 0; n5 < 8; n5++)
				{
					m3 = offset[direction - 1];
					neighbor = currentPixel + offset[direction - 1];
					if (Lp1[neighbor] != (BYTE)0){
						if (((Lp1[currentPixel] == (BYTE)START) == true) && ((initial_departure_direction == 0) == true)){
							initial_departure_direction = direction;
						}
						//else
						//{
						//mnt = (Lp1[currentPixel] == (BYTE)START) == (initial_departure_direction == direction);

						//else if (((Lp1[currentPixel] == (BYTE)START) == true) && ((initial_departure_direction == (direction)) == true))
						else if (((Lp1[currentPixel] == (BYTE)START)) && ((initial_departure_direction == (direction))))
						{
							done = 1;
							found_next_pixel = 1;
							break;
						}
						//}
						next_search_direction = next_search_direction_lut[direction - 1];
						found_next_pixel = 1;
						numPixels = numPixels + 1;
						//if (numPixels > scratch.size())
						//	scratch[2 * scratch.size()] = 0;
						// scratch[numPixels] = neighbor;
						scratch.push_back(neighbor); // to define by the adjacent element data of 8-connected region
						if (Lp1[neighbor] != (BYTE)START)
							Lp1[neighbor] = (BYTE)BOUNDARY; // to verify that it is not starting, then note how the boundary 
						currentPixel = neighbor;
						break;
					}
					direction = next_direction_lut[direction - 1];
				}
				if (!found_next_pixel){
					numPixels = 2;
					scratch[2] = scratch[1];
					done = 1;
				}
			}

			/*convert indices into the column and row numbers*/
			/*put in В, В1*/
			for (int i1 = 0; i1 < scratch.size(); i1++)
			{
				//  for (int i2 = 0; i2 < 1; i2++)
				// {
				//printf("%5d",scratch[i1]);
				ind2sub(scratch[i1], Param.Width, Param.Height, row, col);

				B.push_back(row);
				B1.push_back(col);
				//printf("%5d %5d\n", row, col);
				// floor(scratch[i1*1 + numPixels] / (sizeof(Lp1) / sizeof(FLOAT))) + 1;// 
				//  modf(scratch[i1 * 1 + numPixels], (float*)(sizeof(Lp1) / sizeof(float)));
			}
		}
	}
	free(Lp1);
	free(Lp2);
	free(Lp2_1);
	free(Lp3);
	free(Lp3_1);
	free(Out);
}
//==================================================================================================
//
//==================================================================================================
int Processing(BYTE *InBuf, BYTE *OutBuf, ProcParam Param)
{
	BYTE *Im_test;
	Im_test = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_test[0], 0, Param.Width*Param.Height);
	BYTE *Im_range_filt;
	Im_range_filt = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *Im_aver;
	Im_aver = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *Im_bin;
	Im_bin = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_bin[0], 0, Param.Width*Param.Height);
	BYTE *Im_bin1;
	Im_bin1 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_bin1[0], 0, Param.Width*Param.Height);

	/*Step 1*/
	//returns the array Im_test, where each output pixel contains the range value (maximum value − minimum value) 
	//of the 3-by-3 neighborhood around the corresponding pixel in the input image
	const int vCONTRAST_CELL_X = 3;
	const int vCONTRAST_CELL_Y = 3;
	unsigned char vmin, vmax;
	for (int j = 0; j < Param.Height - vCONTRAST_CELL_X / 2; j++)
	{
		for (int i = 0; i < Param.Width - vCONTRAST_CELL_X / 2; i++)
		{
			EstimateMaxMin(InBuf, Param, i, j, vCONTRAST_CELL_X, vCONTRAST_CELL_Y, &vmin, &vmax);
			//us_min_matrix[x_cell_max*j_c + i_c] = vmin;
			//us_max_matrix[x_cell_max*j_c + i_c] = vmax;
			vmax;
			vmin;
			Im_test[(j + vCONTRAST_CELL_X / 2)*Param.Width + (i + vCONTRAST_CELL_X / 2)] = vmax - vmin; 
		}
	}

	/*Step 1.1 Convert to binary image*/
	memcpy(Im_range_filt, Im_test, Param.Width*Param.Height);

	// normalize image
	int t;
	t = Norm(Im_range_filt, Param, 30);
	if (t == 1) { return 0; } // error

	//Average image
	MedianFilter(Im_range_filt, Im_aver, Param, 16);

	//Difference images & Binarisation
	for (int j_c = 0; j_c < Param.Height; j_c++)
	{
		for (int i_c = 0; i_c < Param.Width; i_c++)
		{			
			if ((Im_aver[j_c*Param.Width + i_c] < Im_range_filt[j_c*Param.Width + i_c]))
				Im_bin[j_c*Param.Width + i_c] = 0;
			if ((Im_aver[j_c*Param.Width + i_c] - Im_range_filt[j_c*Param.Width + i_c]) > 0)
				Im_bin[j_c*Param.Width + i_c] = (BYTE)(Im_aver[j_c*Param.Width + i_c] - Im_range_filt[j_c*Param.Width + i_c]);
			if (Im_bin[j_c*Param.Width + i_c] > 0){
				Im_bin1[j_c*Param.Width + i_c] = 255;
			}
			else
				Im_bin1[j_c*Param.Width + i_c] = 0;
		}
	}

	/*Step 2 Morphological Processing*/
	BYTE *Im_bin5;
	Im_bin5 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	BYTE *Im_bin56;
	Im_bin56 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_bin56[0], 0, Param.Height*Param.Width);

	//Morfological operations (dilate, erode)
	// Appliig dilate, erode with mask which inside this funtion
	MorfOperations(Im_bin1, Im_bin5, 6, Param);

	/*Step 3 Label connected components in a binary image*/
	BYTE *Label;
	Label = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Label[0], 0, Param.Height*Param.Width);
	BYTE *Im_Label;
	Im_Label = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_Label[0], 0, Param.Height*Param.Width);
	BYTE *Im_bin5_copy;
	Im_bin5_copy = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_bin5_copy[0], 0, Param.Height*Param.Width);

	int Num;
	std::vector<int> v_ind, v_y;

	memcpy(Im_bin5_copy, Im_bin5, Param.Height*Param.Width);
	// Returns a matrix Label, of the same size as InBuf(bin), containing labels for the connected objects in InBuf(bin)
	Num = Bwlabel(Im_bin5_copy, Label, Param);

	int max = 0;
	int m = 0, id = 0;
	BYTE *Im_bin7;
	Im_bin7 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));
	memset(&Im_bin7[0], 0, Param.Height*Param.Width);
	// Note connected components
	if (Num > 0)
	{
		for (int k = 1; k <= Num; k++)
		{
			for (int j = 0; j < Param.Height; j++)
			{
				for (int i = 0; i < Param.Width; i++)
				{
					if (Label[j*Param.Width + i] == k)
					{
						v_ind.push_back(j*Param.Width + i);
					}
				}
			}
			m = v_ind.size();
			if (max < m)
			{
				id = k;
				max = m;
			}
			v_ind.clear();
		}
		for (int j = 0; j < Param.Height; j++)
		{
			for (int i = 0; i < Param.Width; i++)
			{
				if (Label[j*Param.Width + i] == id)
				{
					Im_Label[j*Param.Width + i] = 255;
				}
			}
		}
		/*Step 4 Define the boundary contour*/
		std::vector<int> B, B1;
		std::vector<int> scratch;

		// Scratch contains index of boundary contour
		Boundary(Im_Label, Param, B, B1, scratch); 

		for (int ib = 0; ib < scratch.size(); ib++)
		{
			Im_bin7[scratch[ib]] = 255;
		}

		// Define started point from with to start flood fill
		int s = 0, s1 = 0, t = 1;
		for (int j_c = 0; j_c < Param.Height; j_c++)
		{
			for (int i_c = 0; i_c < Param.Width; i_c++)
			{
				if (Im_Label[j_c*Param.Width + i_c] == 255)
				{
					s = s + j_c;
					s1 = s1 + i_c;
					t++;
				}
			}
		}
		int yc = s / t; // (xc, yc) started point
		int xc = s1 / t;

		//if the started point is found and inside the contour
		if ((xc != 0) || (yc != 0)){
			Position start_at;
			start_at.col = xc;
			start_at.row = yc;
			BYTE fill_color = 255;

			flood_fill(Im_bin7, Param, start_at, fill_color); /*the flood fill contour line-by-line method*/
		}
	}
	else
	{
		memset(&Im_bin7[0], 0, Param.Width*Param.Height);
	}

	// convert to binary
	BYTE *Im_bin8;
	Im_bin8 = (BYTE *)malloc(Param.Width*Param.Height*sizeof(BYTE));

	for (int j_c = 0; j_c < Param.Height; j_c++)
	{
		for (int i_c = 0; i_c < Param.Width; i_c++)
		{

			if (Im_bin7[j_c*Param.Width + i_c] == (BYTE)255){
				Im_bin8[(j_c)*Param.Width + i_c] = InBuf[(j_c)*Param.Width + i_c];
			}
			else
				Im_bin8[(j_c)*Param.Width + i_c] = (BYTE)0;
		}
	}

	// result image
	for (int j_c = 0; j_c < Param.Height; j_c++)
	{
		for (int i_c = 0; i_c < Param.Width; i_c++)
		{
			Im_bin56[j_c*Param.Width + i_c] = 255 - Im_bin8[j_c*Param.Width + i_c];
		}
	}

	// can to show any image
	memcpy(OutBuf, Im_Label /*Im_bin56*/, Param.Height * Param.Width);

	free(Im_range_filt);
	free(Im_aver);
	free(Im_bin);
	free(Im_bin1);
	free(Im_test);
	free(Im_bin5);
	free(Im_bin8);
	free(Im_bin56);
	free(Im_bin7);
	free(Label);
	free(Im_Label);
	free(Im_bin5_copy);

	return 0;
}