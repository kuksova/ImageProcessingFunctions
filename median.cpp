#include "stdafx.h"
#include <stdlib.h>
#include <string.h>

#define Filtr_size  3
typedef unsigned char       BYTE;

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define MASK1_WIDTH 3
#define MASK1_HEIGHT 3

#define Filtr_size  3
#define Filtr_iterations  6
typedef struct
{
	int Width;
	int Height;
}  ProcParam, *pProcParam;

unsigned char CalcMedian(BYTE *uc_array, int f_size)
{
	int i1, changes;
	unsigned char *M = uc_array;
	unsigned char Temp;
	do
	{
		changes = 0;
		for (i1 = 0; i1<(f_size - 1); i1++)
		{
			if (M[i1]>M[i1 + 1])
			{
				Temp = M[i1 + 1];
				M[i1 + 1] = M[i1];
				M[i1] = Temp;
				changes++;
			}
		}
	} while (changes);
	return M[f_size / 2];
}


void MedainFiltr(BYTE *InBuf, BYTE *OutBuf, ProcParam Param)
{
	// Create temporary array
	unsigned char *us_matrix = NULL;
	us_matrix = (BYTE*)malloc(sizeof(BYTE) * (Param.Height + 10)*(Param.Width + 10));

	memcpy(us_matrix, InBuf, Param.Height*Param.Width);

	unsigned char mask[Filtr_size*Filtr_size]; // маска
	int i_c, j_c;

	for ( j_c =0; j_c < Param.Height; j_c++)
	{
		us_matrix[Param.Height*(j_c)+(0)] = 0;
		us_matrix[Param.Height*(j_c)+(Param.Height - 1)] = 0;
		for (i_c = 0; i_c < Param.Width; i_c++)
		{
			int cnt = 0;
			for (int j = 0; j <Filtr_size; j++)
			{
				for (int i =0; i < Filtr_size; i++)
				{
					mask[cnt] = us_matrix[Param.Width*(j_c + j) + (i_c + i)];
					//printf("%3d", uc_array[cnt]);
					cnt++;
				}
			}
			OutBuf[Param.Width*(j_c)+(i_c)] = CalcMedian(&mask[0], cnt);
			//printf("%5d", OutBuf[j_c*(Param.Width) + i_c]);
		}
		//printf("\n");
	}

	free(us_matrix);
}
int main(int argc, char* argv[])
{
	ProcParam Param;
	Param.Height = 12;
	Param.Width = 7;

	unsigned char *IM = NULL;
	IM = (unsigned char*)malloc(sizeof(unsigned char) * (Param.Width)*(Param.Height));
	for (int i = 0; i < Param.Height; i++){
		for (int j = 0; j < Param.Width; j++){
			IM[i*(Param.Width) + j] = 0 + rand() % 100;
			printf("%5d", IM[i*(Param.Width) + j]);
		}
		printf("\n");
	}
	//Result image
	unsigned char *Im_aver;
	Im_aver = (unsigned char *)malloc(Param.Width*Param.Height*sizeof(unsigned char));

	MedainFiltr(IM, Im_aver, Param);
	printf("\n");
	for (int i = 0; i < Param.Height; i++){
		for (int j = 0; j < Param.Width; j++){
			printf("%5d", Im_aver[i*(Param.Width) + j]);
		}
		printf("\n");
	}

	free(IM);
	free(Im_aver);
	getchar();
	return 0;
}