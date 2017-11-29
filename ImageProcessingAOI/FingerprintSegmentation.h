//FingerprintSegmentation
//AOI 
// 26 october 2017
typedef unsigned char       BYTE;

typedef struct
{
	int Width;
	int Height;
} ProcParam, *pProcParam;

typedef struct {
	int r[256];
	int n;
} color_histo_t;


//==================================================================================================
//
//==================================================================================================

int Processing(BYTE *InBuf, BYTE *OutBuf, ProcParam Param);
void MedianFilter(BYTE *InBuf, BYTE *OutBuf, ProcParam Param, int size);
void MorfOperations(BYTE *Img, BYTE *ImgD, int H, ProcParam Param);
int Bwlabel(BYTE * mas, BYTE * Label, ProcParam Param)