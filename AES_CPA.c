#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define _FOLD_ "C:\\"
#define TraceFN "AES.traces"
#define AlignedTraceFN "AlignedAES.traces"
#define PlaintextFN "plaintext.txt"
#define CiphertextFN "ciphertext.txt"
#define startpoint 22551
#define endpoint 31050

static unsigned char SBOX[256] = {
0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

double cov(float* x, float* y, int size) {
	double Sxy = 0, Sx = 0, Sy = 0;
	int i;

	for (i = 0; i < size; i++) {
		Sxy += x[i] * y[i];	//E(XY)
		Sx += x[i];
		Sy += y[i];
	}
	//E(XY)-E(X)E(Y)
	return (Sxy - Sx * Sy / (double)size) / (double)size;
	//return Sxy / (double)size - (Sx / (double)size) * (Sy / (double)size);
}
double corr(float* x, float* y, int size) {
	double Sxy = 0, Sx = 0, Sy = 0, Sxx = 0, Syy = 0;	//var(X) = E(X^2) - E(X)^2
	int i;

	for (i = 0; i < size; i++) {
		Sxy += x[i] * y[i];	//E(XY)
		Sx += x[i];
		Sy += y[i];
		Sxx += x[i] * x[i];
		Syy += y[i] * y[i];
	}
	//상관계수
	//return ((Sxy - Sx * Sy / (double)size) / (double)size)/sqrt((Sxx/(double)size-Sx*Sx/(double)size/(double)size)* (Syy / (double)size - Sy * Sy / (double)size / (double)size));
	return ((double)size * Sxy - Sx * Sy) / sqrt(((double)size * Sxx - Sx * Sx) * ((double)size * Syy - Sy * Sy));
}

void subalign(float *data, float *data1, int windowsize, int stepsize, int threshold, int TraceLength) {
	//data 배열에 저장되어 있는 전력파형을 기준으로 data1 배열에 저장되어 있는 전력파형을 정렬

	int i,j,size,maxcovpos,k;
	float* x, * y;
	double covval,maxcov;
	

	for (i = 0; i < TraceLength-windowsize; i+=stepsize) { //stepsize만큼 띄우면서 정렬
		maxcovpos = 0; //파형을 흔들 때 마다 얼만큼 흔들지 사용하는 변수
		maxcov = 0;
		for (j = -threshold; j < threshold; j++) { //파형 j를 좌우로 흔들면서 정렬
			//파형을 좌측으로 정렬할 때
			if (j < 0) {
				x = data + i;
				y = data1 + i - j;
				size = windowsize +j;
			}
			//파형을 우측으로 정렬할 때
			else {
				x = data + i + j;
				y = data1 + i;
				size = windowsize - j;
			}

			covval=cov(x, y, size); //covariance 값 계산
			//covariance 값의 최대값을 찾음
			if (covval > maxcov) {

				maxcovpos = j; //covariance 값이 최대값이 되게끔하는 j 값 저장(j값을 기준으로 data1의 파형을 data 파형으로 맞춤)
				maxcov = covval; //covariance 최대 값 저장
			}
		}
		//파형을 좌측으로 정렬
		if (maxcovpos < 0) {
			for (k = i; k < (TraceLength + maxcovpos); k++) {
				data1[k] = data1[k - maxcovpos]; //data1+i-j
			}
		}
		//파형을 우측으로 정렬
		else {
			//우측으로 정렬할 때는 마지막부터 우측으로 하나씩 밀어냄
			for (k = (TraceLength - maxcovpos - 1); k >= i; k--) {
				data1[k + maxcovpos] = data1[k]; //data1+i -> data1+i+maxcovpos
			}
		}
	}
}

void Alignment() {
	int windowsize = 500;	//부분부분 정렬을 맞추고 싶은 파형의 길이
	int stepsize = 450;	//한 부분의 정렬을 맞춘 후에 몇 포인트씩 이동 후 다시 정렬을 맞출 것인지 결정
	int threshold = 100;	//좌우로 얼마나 흔들면서 covariance값을 계산해서 최대값과 이동할 포인트 수를 계산
	char buf[256];
	int err,TraceLength,TraceNum,i;
	FILE* rfp,*wfp;

	float* data, * data1;

	sprintf_s(buf, 256 * sizeof(char), "%s%s", _FOLD_, TraceFN); //해당경로안에 있는 문자열이 buf안에 들어옴
	if ((err = fopen_s(&rfp, buf, "rb")))
	{
		printf("File Open Error!!\n");
	}

	sprintf_s(buf, 256 * sizeof(char), "%s%s", _FOLD_, AlignedTraceFN); 
	if ((err = fopen_s(&wfp, buf, "wb")))
	{
		printf("File Open Error!!\n");
	}

	//첫번 째 파형은 동일한 데이터를 읽고 동일한 데이터를 씀
	fread(&TraceLength, sizeof(int), 1, rfp);
	fwrite(&TraceLength, sizeof(int), 1, wfp);
	fread(&TraceNum, sizeof(int), 1, rfp);
	fwrite(&TraceNum, sizeof(int), 1, wfp);

	data = (float*)calloc(TraceLength, sizeof(float));
	data1 = (float*)calloc(TraceLength, sizeof(float));
	fread(data, sizeof(float), TraceLength, rfp);
	fwrite(data, sizeof(float), TraceLength, wfp);

	//두번 째 파형부터 읽음
	for (i = 1; i < TraceNum; i++) {
		fread(data1, sizeof(float), TraceLength, rfp);
		subalign(data, data1, windowsize, stepsize, threshold,TraceLength); //데이터를 정렬하여 data1에 저장
		fwrite(data1, sizeof(float), TraceLength, wfp);
	}

	fclose(rfp);
	fclose(wfp);

	free(data);
	free(data1);
}
void CPA() {
	unsigned char** plaintext = NULL; //2000개의 16byte의 평문을 2000*16 크기의 배열에 저장
	float** data; //정렬된 파형을 한번에 메모리에 올려서 작업
	double* Sx, * Sxx, * Sxy, * corrT;
	double Sy, Syy,max;
	unsigned char temp[34],x,y,iv, hw_iv;
	char buf[256];
	int err,TraceLength,TraceNum,i,j,k,key,maxkey;
	FILE* rfp,* wfp;

	sprintf_s(buf, 256 * sizeof(char), "%s%s", _FOLD_, AlignedTraceFN);
	if ((err = fopen_s(&rfp, buf, "rb")))
	{
		printf("File Open Error1!!\n");
	}
	fread(&TraceLength, sizeof(int), 1, rfp);
	fread(&TraceNum, sizeof(int), 1, rfp);

	//2차원 포인터 배열에 파형을 다 읽어옴
	data = (float**)calloc(TraceNum, sizeof(float*));
	for (i = 0; i < TraceNum; i++) {
		data[i] = (float*)calloc(TraceLength, sizeof(float));
	}
	for (i = 0; i < TraceNum; i++) {
		fread(data[i], sizeof(float), TraceLength, rfp);
	}
	fclose(rfp);

	//평문 읽기(아스키코드로 되어있는 평문을 한 줄씩 읽으면 34bytes(32+2)에 해당)
	sprintf_s(buf, 256 * sizeof(char), "%s%s", _FOLD_, PlaintextFN);
	if ((err = fopen_s(&rfp, buf, "r")))
	{
		printf("File Open Error2!!\n");
	}
	plaintext = (unsigned char**)calloc(TraceNum, sizeof(unsigned char*));
	for (i = 0; i < TraceNum; i++) {
		plaintext[i] = (unsigned char*)calloc(16, sizeof(unsigned char));
	}
	for (i = 0; i < TraceNum; i++) {
		fgets(temp, 34, rfp); //-->16bytes로 바꿔서 plaintext[i]에 저장 필요
		//printf("%02x ",temp);
		for (j = 0; j < 16; j++) {
			//순차적으로 문자열 처리 ex)x=15,y=16...
			x = temp[2 * j];
			y = temp[2 * j + 1];

			if (x >= 'A' && x <= 'Z')x = x - 'A'+10; //'0'~'9','A'~'F','a'~'f'
			else if(x >= 'a' && x <= 'z')x = x - 'a' + 10;
			else if (x >= '0' && x <= '9')x -= '0';
			if (y >= 'A' && y <= 'Z')y = y - 'A' + 10; //'0'~'9','A'~'F','a'~'f'
			else if (y >= 'a' && y <= 'z')y = y - 'a' + 10;
			else if (y >= '0' && y <= '9')y -= '0';
			plaintext[i][j] = x * 16 + y;

		}
	}
	fclose(rfp);

	//배열 선언
	Sx = (double*)calloc(TraceLength, sizeof(double));
	Sxx = (double*)calloc(TraceLength, sizeof(double));
	Sxy = (double*)calloc(TraceLength, sizeof(double));
	corrT = (double*)calloc(TraceLength, sizeof(double));


	for (i = 0; i < TraceNum; i++) {
		for (j = startpoint; j < endpoint; j++) {
			Sx[j] += data[i][j];
			Sxx[j] += data[i][j] * data[i][j];
		}
	}

	//키 추측
	for (i = 0; i < 16; i++) {
		max = 0;
		maxkey = 0;
		for (key = 0; key < 256; key++) {
			Sy = 0;
			Syy = 0;
			memset(Sxy, 0, sizeof(double) * TraceLength);
			for (j = 0; j < TraceNum; j++) { //평문의 개수만큼 반복
				iv = SBOX[plaintext[j][i] ^ key]; //추측한 키마다 SBOX의 중간값을 알아냄
				hw_iv = 0;
				//해밍 웨이트 값 계산(1의 개수를 계산)
				for (k = 0; k < 8; k++)hw_iv += ((iv >> k) & 1);
				Sy += hw_iv;
				Syy += hw_iv * hw_iv;
				for (k = startpoint; k < endpoint; k++) {
					Sxy[k] += hw_iv * data[j][k];
				} 
			}
			//상관계수 파형 계산
			for (k = startpoint; k < endpoint; k++) {
				corrT[k] = ((double)TraceNum * Sxy[k] - Sx[k]*Sy)/sqrt(((double)TraceNum*Sxx[k]-Sx[k]*Sx[k])* ((double)TraceNum * Syy - Sy * Sy));
				if (fabs(corrT[k]) > max) { //상관계수 최대값 구하기
					maxkey = key;
					max = fabs(corrT[k]);
				}
			
			}

			sprintf_s(buf, 256 * sizeof(char), "%scorrtrace\\%02dth_block_%02x.corrtrace",_FOLD_, i,key);
			if ((err = fopen_s(&wfp, buf, "wb")))
			{
				printf("File Open Error3!!\n");
			}
			fwrite(corrT, sizeof(double), TraceLength, wfp);
			fclose(wfp);
			printf(".");
		}
		printf("%02dth_block : maxkey(%02x),maxcorr(%lf)\n",i,maxkey,max);
	}
	free(Sx);
	free(Sxx);
	free(Sxy);
	free(corrT);
	free(data);
	free(plaintext);
}

int main() {
	//float X[10] = { 1,2,3,5,6,4,6,7,6,5 };
	//float Y[10] = { 2,3,6,8,10,7,11,12,11,9 };

	//printf("covariance: %lf, correlation coefficient: %lf\n", cov(X, Y, 10), corr(X, Y, 10));
	//Alignment();
	CPA();
	return 0;
}