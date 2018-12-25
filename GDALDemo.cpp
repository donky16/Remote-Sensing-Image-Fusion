// GDALDemo.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "gdal/gdal_priv.h"
#include<iostream>
#include<math.h>
#include<time.h>
#pragma comment(lib, "gdal_i.lib")

using namespace std;
// *****************************************�ָ���*************************************************
/*
	�޸�һ����ͼ�������ֵ
	RasterIO() readһ��һ�еĴ洢��bufftmp�е�
*/
void modifyPix()
{
	GDALDataset* poSrcDS;
	GDALDataset* poDstDS;
	int imgXlen, imgYlen;
	char* srcPath = "Data\\pictures\\Taeyeon.jpg";
	char* dstPath = "Data\\output\\result_modify.tif";
	GByte* bufftmp;
	int bandNum;
	// �޸�ͼ��block����
	int blackStartX = 300, blackStartY = 300;
	int whiteStartX = 500, whiteStartY = 500;
	int blackXlen = 100;
	int blackYlen = 50;
	int whiteXlen = 50;
	int whiteYlen = 100;

	GDALAllRegister();

	poSrcDS = (GDALDataset*)GDALOpenShared(srcPath, GA_ReadOnly);

	imgXlen = poSrcDS->GetRasterXSize();
	imgYlen = poSrcDS->GetRasterYSize();
	bandNum = poSrcDS->GetRasterCount();

	cout << "Image X Length:" << imgXlen << endl;
	cout << "Image Y Length:" << imgYlen << endl;
	cout << "Band Number:" << bandNum << endl;
	// �����ڴ�
	bufftmp = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));

	// ��������
	poDstDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(dstPath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	for (int i = 0; i < bandNum; i++)
	{
		// ��ȡͼ�����ݵ��ڴ�
		poSrcDS->GetRasterBand(i + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, bufftmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		// ���ڴ�������д��ͼ��
		poDstDS->GetRasterBand(i + 1)->RasterIO(GF_Write, 0, 0, imgXlen, imgYlen, bufftmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		// �޸��ڴ��в�������
		for (int k = 0; k < blackXlen; k++)
		{
			for (int j = 0; j < blackYlen; j++)
			{
				bufftmp[k*blackXlen + j] = (GByte)0;
			}
		}
		poDstDS->GetRasterBand(i + 1)->RasterIO(GF_Write, blackStartX, blackStartY, blackXlen, blackYlen, bufftmp, blackXlen, blackYlen, GDT_Byte, 0, 0);
		for (int n = 0; n < whiteXlen; n++)
		{
			for (int m = 0; m < whiteYlen; m++)
			{
				bufftmp[n*whiteXlen + m] = (GByte)255;

			}
		}
		poDstDS->GetRasterBand(i + 1)->RasterIO(GF_Write, whiteStartX, whiteStartY, whiteXlen, whiteYlen, bufftmp, whiteXlen, whiteYlen, GDT_Byte, 0, 0);
	}
	CPLFree(bufftmp);

	GDALClose(poSrcDS);
	GDALClose(poDstDS);
}
// *****************************************�ָ���*************************************************
/*
	ͼ��ϳɺ�������supermanͼ�����ɫ�����۵�����supermanǶ��spaceͼ����
	��ʽ��
		1.������ͼ�������ͨ���ֱ�洢��supermanBufftmp[3]��spaceBufftmp[3]��
		2.�Ƚ�superman������ֵ�����̱���ȥ������������ֵ��ֵ��spaceͼ������ͬ��λ��
		3.��spaceBufftemp����д����ͼ�����ɺϳ�ͼ��
*/
void imageSynthesis()
{
	GDALDataset* supermanSrcDS;
	GDALDataset* spaceSrcDS;
	GDALDataset* resultDS;

	char* supermanPath = "Data\\pictures\\superman.jpg";
	char* spacePath = "Data\\pictures\\space.jpg";
	char* supermanAndSpacePath = "Data\\output\\superman-on-the-space.tif";

	int imgXlen, imgYlen, bandNum;

	GByte RValue, GValue, BValue;

	GDALAllRegister();

	supermanSrcDS = (GDALDataset*)GDALOpenShared(supermanPath, GA_ReadOnly);
	spaceSrcDS = (GDALDataset*)GDALOpenShared(spacePath, GA_ReadOnly);

	imgXlen = spaceSrcDS->GetRasterXSize();
	imgYlen = spaceSrcDS->GetRasterYSize();
	bandNum = spaceSrcDS->GetRasterCount();

	GByte* supermanBufftmp[3];		// bandNum��Ϊ3
	GByte* spaceBufftmp[3];
	// ͼƬ����д��bufftmp
	for (int n = 0; n < bandNum; n++)
	{
		supermanBufftmp[n] = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
		supermanSrcDS->GetRasterBand(n + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, supermanBufftmp[n], imgXlen, imgYlen, GDT_Byte, 0, 0);
		spaceBufftmp[n] = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
		spaceSrcDS->GetRasterBand(n + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, spaceBufftmp[n], imgXlen, imgYlen, GDT_Byte, 0, 0);
	}
	// ��superman�г��̱���������д��space��
	int flag;
	for (int j = 0; j < imgYlen; j++)
	{
		for (int i = 0; i < imgXlen; i++)
		{
			flag = 1;
			RValue = (GByte)supermanBufftmp[0][i*imgYlen + j];
			GValue = (GByte)supermanBufftmp[1][i*imgYlen + j];
			BValue = (GByte)supermanBufftmp[2][i*imgYlen + j];
			// ���̱�����Χ����
			if (RValue > 50 && RValue < 160 && GValue >130 && GValue < 220 && BValue >50 && BValue < 150)
			{
				flag = 0;
			}
			if (flag)
			{
				spaceBufftmp[0][i*imgYlen + j] = RValue;
				spaceBufftmp[1][i*imgYlen + j] = GValue;
				spaceBufftmp[2][i*imgYlen + j] = BValue;
			}
		}
	}
	// ��������
	resultDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(supermanAndSpacePath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	// д��ϳ�ͼ���ͷ�
	for (int m = 0; m < 3; m++)
	{
		resultDS->GetRasterBand(m + 1)->RasterIO(GF_Write, 0, 0, imgXlen, imgYlen, spaceBufftmp[m], imgXlen, imgYlen, GDT_Byte, 0, 0);
		CPLFree(supermanBufftmp[m]);
		CPLFree(spaceBufftmp[m]);
	}
	GDALClose(supermanSrcDS);
	GDALClose(spaceSrcDS);
	GDALClose(resultDS);

}
// *****************************************�ָ���*************************************************
/*
	�����˲�ʵ��
*/
GByte* getConvolutionSum(GByte* bufftemp, int Xlen, int Ylen, double convolutionKernel[], int kernelLen, int param, int offset)
{
	/*
		��ȡ�������صľ���͵�ֵ��һ�����ε�����ֵ
		params: Դͼ�����ݣ�ԭͼ���X��Y���;���ˣ�����˵ĳ���(���)���������Ҫ����������ƫ����
		return: ��������ͼ����ͼ������GByte���͵�ָ��

		return ͼ����� * ����
	*/
	GByte* sumData = (GByte*)CPLMalloc(Xlen*Ylen * sizeof(GByte));
	double tmp = 0;		// �洢ÿһ�����صľ����
	for (int i = 0; i < Xlen; i++)
	{
		for (int j = 0; j < Ylen; j++)
		{
			tmp = 0;
			// ɸѡ���߽粢��0
			if ((i == 0) || (i == Xlen - 1) || (j == 0) || (j == Ylen - 1))
			{
				sumData[i*Ylen + j] = tmp;
			}
			else
			{
				for (int m = 0; m < kernelLen; m++)
				{
					for (int n = 0; n < kernelLen; n++)
					{
						/*
							���������������[(kernelLen - 1) / 2, (kernelLen - 1) / 2]
							���������һ����Ϊ(m,n)����(XDiff,YDiff)Ϊ�������λ��
							�������λ�ã��������ͼ������һ���ص���������ֵ
						*/
						int XDiff = n - (kernelLen - 1) / 2;
						int YDiff = m - (kernelLen - 1) / 2;
						tmp += convolutionKernel[m*kernelLen + n] * bufftemp[i*Ylen + j + XDiff*Ylen + YDiff];
					}
				}
				tmp = tmp / param + offset;
				// �ضϣ��������ֵ���ڣ�0,255��֮��
				if (tmp > 255)
				{
					sumData[i*Ylen + j] = 255;
				}
				else if (tmp < 0)
				{
					sumData[i*Ylen + j] = 0;
				}
				else
				{
					sumData[i*Ylen + j] = (GByte)round(tmp);
				}
			}

		}
	}
	return sumData;
}
void convolutionConvert(char* poSrcPath, char* poResultPath, double convolutionKernel[], int kernelLen, int param, int offset)
{
	/*
		����ԭͼ�񣬾���ˣ�����˵ĳ���(���)���������Ҫ����������ƫ�����õ�ת��֮���ͼ��
	*/
	GDALAllRegister();

	GDALDataset* poSrcDS = (GDALDataset*)GDALOpenShared(poSrcPath, GA_ReadOnly);

	//cout << "Loading lean.jpg..." << endl;

	int imgXlen = poSrcDS->GetRasterXSize();
	int imgYlen = poSrcDS->GetRasterYSize();
	int bandNum = poSrcDS->GetRasterCount();
	cout << "X len: " << imgXlen << endl;
	cout << "Y len: " << imgYlen << endl;
	cout << "Band num : " << bandNum << endl;

	GDALDataset* resultDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(poResultPath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	GByte* bufftmpSrc = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
	GByte* bufftmpResult = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));


	for (int b = 0; b < bandNum; b++)
	{
		poSrcDS->GetRasterBand(b + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, bufftmpSrc, imgXlen, imgYlen, GDT_Byte, 0, 0);
		bufftmpResult = getConvolutionSum(bufftmpSrc, imgXlen, imgYlen, convolutionKernel, kernelLen, param, offset);
		resultDS->GetRasterBand(b + 1)->RasterIO(GF_Write, 0, 0, imgXlen, imgYlen, bufftmpResult, imgXlen, imgYlen, GDT_Byte, 0, 0);
	}

	CPLFree(bufftmpSrc);
	CPLFree(bufftmpResult);
	GDALClose(poSrcDS);
	GDALClose(resultDS);
}

void convolutionMain()
{
	/*
		�����˲������������ݲ�ͬ����˻�ȡ��ͬ��ͼ��
	*/
	char* poSrcPath = "Data\\pictures\\lena.jpg";
	char* resultPath = "Data\\output\\convolution.jpg";
	int kernekLen = 0;

	/*  �����1  ��ֵģ��
		 0 , 1 ,  0
		 1 , 1 ,  1		/ 5  ��param=5��
		 0 , 1 ,  0
	*/
	resultPath = "Data\\output\\convolution-1.jpg";
	double convolutionKernel_1[9] = { 0, 1, 0, 1, 1, 1, 0, 1, 0 };
	kernekLen = sqrt(sizeof(convolutionKernel_1) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_1, kernekLen, 5, 0);

	/*  �����2  �˶�ģ��
		1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0,    / 5	(param=5)
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1,
	*/
	resultPath = "Data\\output\\convolution-2.jpg";
	double convolutionKernel_2[25] = {
		1,0,0,0,0,
		0,1,0,0,0,
		0,0,1,0,0,
		0,0,0,1,0,
		0,0,0,0,1 };
	kernekLen = sqrt(sizeof(convolutionKernel_2) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_2, kernekLen, 5, 0);

	/* �����3 ��Ե���
		-1,-1,-1,
		-1, 8,-1,
		-1,-1,-1,

	*/
	resultPath = "Data\\output\\convolution-3.jpg";
	double convolutionKernel_3[9] = { -1,-1,-1,-1,8,-1,-1,-1,-1 };
	kernekLen = sqrt(sizeof(convolutionKernel_3) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_3, kernekLen, 1, 0);

	/* �����4  ��
		-1,-1,-1,
		-1, 9,-1,
		-1,-1,-1,
	*/
	resultPath = "Data\\output\\convolution-4.jpg";
	double convolutionKernel_4[9] = { -1,-1,-1,-1,9,-1,-1,-1,-1 };
	kernekLen = sqrt(sizeof(convolutionKernel_4) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_4, kernekLen, 1, 0);
	/* �����5  ����
		-1,-1, 0,
		-1, 0, 1,		offset=128
		 0, 1, 1,
	*/
	resultPath = "Data\\output\\convolution-5.jpg";
	double convolutionKernel_5[9] = { -1,-1,0,-1,0,1,0,1,1 };
	kernekLen = sqrt(sizeof(convolutionKernel_5) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_5, kernekLen, 1, 128);
	/* �����6 ��˹ģ��
		0.0120,0.1253,0.2736,0.1253,0.0120,
		0.1253,1.3054,2.8514,1.3054,0.1253,
		0.2736,2.8514,6.2279,2.8514,0.2736,
		0.1253,1.3054,2.8514,1.3054,0.1253,
		0.0120,0.1253,0.2736,0.1253,0.0120,
	*/
	resultPath = "Data\\output\\convolution-6.jpg";
	double convolutionKernel_6[25] = {
		0.0120,0.1253,0.2736,0.1253,0.0120,
		0.1253,1.3054,2.8514,1.3054,0.1253,
		0.2736,2.8514,6.2279,2.8514,0.2736,
		0.1253,1.3054,2.8514,1.3054,0.1253,
		0.0120,0.1253,0.2736,0.1253,0.0120
	};
	kernekLen = sqrt(sizeof(convolutionKernel_6) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_6, kernekLen, 25, 0);
}
void IHSImageFusion(float **RGB, float *panBuffTmp, float **IHS, float **resultRGB, int imgXlen, int imgYlen)
{
	/*
		IHSͼ���ںϺ����㷨
		param: RGB-�����ͼ���RGBֵ panBuffTmp-ȫɫͼ�������ֵ IHS-�����ⶨ���IHS 
		       resultRGB-�洢�ںϺ��ͼ������ֵ imgXlen,imgYlen-ͼ��ĳ��ȺͿ��	
	*/
	float matrix[3][3] = { { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 },{ -sqrt(2.0) / 6.0, -sqrt(2.0) / 6.0, sqrt(2.0) / 3.0 },{ 1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0 } };
	float matrixReverse[3][3] = { { 1.0, -1.0 / sqrt(2.0), 1.0 / sqrt(2.0) },{ 1.0, -1.0 / sqrt(2.0), -1.0 / sqrt(2.0) },{ 1.0, sqrt(2.0), 0 } };
	// RGB2IHS
	for (int m = 0; m < 3; m++)
	{
		for (int n = 0; n < imgXlen*imgYlen; n++)
		{
			IHS[m][n] = 0;	// ��ʼ��Ϊ0��������Ҫ��ε���
			for (int k = 0; k < 3; k++)
			{
				IHS[m][n] += matrix[m][k] * RGB[k][n];
			}
		}
	}
	// ��panͼ�������û�I
	// IHS[0] = panBuffTmp; ָ�븳ֵ�ᵼ��IHS[0]ָ��panBuffTmp����һ�ε��ô˺����ͻḳֵʧ��
	for (int i = 0; i < imgXlen*imgYlen; i++)
	{
		IHS[0][i] = panBuffTmp[i];
	}
	// ��任
	for (int mm = 0; mm < 3; mm++)
	{
		for (int nn = 0; nn < imgXlen*imgYlen; nn++)
		{
			resultRGB[mm][nn] = 0;	// ��ʼ��Ϊ0��������Ҫ��ε���
			for (int kk = 0; kk < 3; kk++)
			{
				resultRGB[mm][nn] += matrixReverse[mm][kk] * IHS[kk][nn];
			}
		}
	}
}

void IHSImageFusionMain(char* srcMULPath, char* srcPANPath, char* resultPath)
{
	/*
		ң��ͼ���ں�������
		param: �����ͼ��·����ȫɫͼ��·�����ںϽ��ͼ��·��
	*/
	GDALAllRegister();

	GDALDataset* srcMULDS = (GDALDataset*)GDALOpenShared(srcMULPath, GA_ReadOnly);
	GDALDataset* srcPANDS = (GDALDataset*)GDALOpenShared(srcPANPath, GA_ReadOnly);

	int imgMULXlen = srcMULDS->GetRasterXSize();
	int imgMULYlen = srcMULDS->GetRasterYSize();
	int MULBandNum = srcMULDS->GetRasterCount();
	// ������������ͼ�񣬿��Ƿ���ֱ������GDT_Float32�Ĳ��У��ܻ�
	GDALDataset *resDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(resultPath, imgMULXlen, imgMULYlen, MULBandNum, GDT_Byte, NULL);
	
	cout << "MUL image X len: " << imgMULXlen << endl;
	cout << "MUL image Y len: " << imgMULYlen << endl;
	cout << "MUL image Band num : " << MULBandNum << endl;

	int imgPANXlen = srcPANDS->GetRasterXSize();
	int imgPANYlen = srcPANDS->GetRasterYSize();
	int PANBandNum = srcPANDS->GetRasterCount();
	// ��float32���ͣ�ԭͼ��ÿһ�����ص㲢��������
	float* panBuffTmp = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
	srcPANDS->GetRasterBand(PANBandNum)->RasterIO(GF_Read, 0, 0, imgMULXlen, imgMULYlen, panBuffTmp, imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);
	
	cout << "PAN image X len: " << imgPANXlen << endl;
	cout << "PAN image Y len: " << imgPANYlen << endl;
	cout << "PAN image Band num : " << PANBandNum << endl;
	
	// ��ά���飬ÿһ���������鶼����CPLMalloc���м��IHSҲ�ڴ�malloc������ʡmalloc��free��ʱ��
	float **RGB;
	float **resultRGB;
	float **IHS;
	RGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	resultRGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	IHS = (float**)CPLMalloc(MULBandNum * sizeof(float));
	
	// ��ȡ���ݵ�������
	for (int i = 0; i < 3; i++)
	{
		RGB[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
		srcMULDS->GetRasterBand(i + 1)->RasterIO(GF_Read, 0, 0, imgMULXlen, imgMULYlen, RGB[i], imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);	
		resultRGB[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
		IHS[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
	}

	// �����㷨ֱ�Ӽ��㣬����洢��resultRGB��
	IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, imgMULXlen, imgMULYlen);

	// ֱ�ӽ�float��������ֵ�Ž�ȥ�����ﴴ������GDT_Byte�ģ����ﾹȻ����ֱ�ӷŽ�ȥ����
	for (int j = 0; j < 3; j++)
	{
		resDS->GetRasterBand(j + 1)->RasterIO(GF_Write, 0, 0, imgMULXlen, imgMULYlen, resultRGB[j], imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);
	}

	// �ͷ��ڴ�
	CPLFree(panBuffTmp);
	CPLFree(RGB);
	CPLFree(IHS);
	CPLFree(resultRGB);
	GDALClose(srcMULDS);
	GDALClose(srcPANDS);
	GDALClose(resDS);

}

void IHSBigImageFusionByBlock(char *srcMULPath, char *srcPANPath, char * resPath)
{
	/*
		�ֿ�ͼ���ں�
		param: �����ͼ��·����ȫɫͼ��·�����ںϽ��ͼ��·��
	*/
	GDALAllRegister();
	// ��Ĵ�С
	int blockLen = 256;
	GDALDataset* srcMULDS = (GDALDataset*)GDALOpenShared(srcMULPath, GA_ReadOnly);
	GDALDataset* srcPANDS = (GDALDataset*)GDALOpenShared(srcPANPath, GA_ReadOnly);

	int imgMULXlen = srcMULDS->GetRasterXSize();
	int imgMULYlen = srcMULDS->GetRasterYSize();
	int MULBandNum = srcMULDS->GetRasterCount();

	GDALDataset *resDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(resPath, imgMULXlen, imgMULYlen, MULBandNum, GDT_Byte, NULL);

	cout << "MUL image X len: " << imgMULXlen << endl;
	cout << "MUL image Y len: " << imgMULYlen << endl;
	cout << "MUL image Band num : " << MULBandNum << endl;

	int imgPANXlen = srcPANDS->GetRasterXSize();
	int imgPANYlen = srcPANDS->GetRasterYSize();
	int PANBandNum = srcPANDS->GetRasterCount();
	cout << "PAN image X len: " << imgPANXlen << endl;
	cout << "PAN image Y len: " << imgPANYlen << endl;
	cout << "PAN image Band num : " << PANBandNum << endl;

	// ������ұ�ʣ������غ�����ʣ�������
	int blocksX = imgMULXlen / blockLen;
	int blocksY = imgMULYlen / blockLen;
	int restX = imgMULXlen - blocksX * blockLen;
	int restY = imgMULYlen - blocksY * blockLen;
	cout << "Blocks num: X-" << blocksX << "   Y-" << blocksY << endl;
	cout << "Rest pixels: X-" << restX << "   Y-" << restY << endl;

	float *panBuffTmp;
	float **RGB;
	float **resultRGB;
	float **IHS;
	RGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	resultRGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	IHS = (float**)CPLMalloc(MULBandNum * sizeof(float));
	panBuffTmp = (float*)CPLMalloc(blockLen*blockLen * sizeof(float));

	for (int k = 0; k < 3; k++)
	{
		RGB[k] = (float*)CPLMalloc(blockLen*blockLen * sizeof(float));
		IHS[k] = (float*)CPLMalloc(blockLen*blockLen * sizeof(float));
		resultRGB[k] = (float*)CPLMalloc(blockLen*blockLen * sizeof(float));
	}
	int X,Y;
	for (int m = 0; m <= blocksY; m++)
	{
		for (int n = 0; n <= blocksX; n++)
		{
			// ������ʵ�ʴ�С
			if (n == blocksX)
			{
				X = restX;
			}
			else if (m == blocksY)
			{
				Y = restY;
			}
			else
			{
				X = Y = blockLen;
			}
			//cout << "Working..." << endl;
			//cout << "Xoff:" << n*blockLen << "    Yoff:" << m*blockLen << endl;
			//cout << "Block size:" << X << "*" << Y << "=" << X*Y << endl;
			
			// ���ݿ�Ĵ�С����ͼ���ж�ȡ����
			srcPANDS->GetRasterBand(1)->RasterIO(GF_Read, n*blockLen, m*blockLen, X, Y, panBuffTmp, X, Y, GDT_Float32, 0, 0);
			for (int j = 0; j < 3; j++)
			{		
				srcMULDS->GetRasterBand(j + 1)->RasterIO(GF_Read, n*blockLen, m*blockLen, X, Y, RGB[j], X, Y, GDT_Float32, 0, 0);
			}

			// �����㷨�ں�
			IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, X, Y);

			// ����д��ͼ��
			for (int i = 0; i < 3; i++)  
			{
				resDS->GetRasterBand(i + 1)->RasterIO(GF_Write, n*blockLen, m*blockLen, X, Y, resultRGB[i], X, Y, GDT_Float32, 0, 0);
			}
		}
	}
	// �ͷ��ڴ�
	CPLFree(panBuffTmp);
	CPLFree(RGB);
	CPLFree(IHS);
	CPLFree(resultRGB);
	GDALClose(srcMULDS);
	GDALClose(srcPANDS);
	GDALClose(resDS);
}
void IHSBigImageFusionByLine(char *srcMULPath, char *srcPANPath, char * resPath)
{
	/*
		����ͼ���ں�
		param: �����ͼ��·����ȫɫͼ��·�����ںϽ��ͼ��·��
	*/
	GDALAllRegister();

	GDALDataset* srcMULDS = (GDALDataset*)GDALOpenShared(srcMULPath, GA_ReadOnly);
	GDALDataset* srcPANDS = (GDALDataset*)GDALOpenShared(srcPANPath, GA_ReadOnly);

	int imgMULXlen = srcMULDS->GetRasterXSize();
	int imgMULYlen = srcMULDS->GetRasterYSize();
	int MULBandNum = srcMULDS->GetRasterCount();
	int lineHeight = 256;		// �и߶�
	int lineLen = imgMULXlen;	// �г���

	GDALDataset *resDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(resPath, imgMULXlen, imgMULYlen, MULBandNum, GDT_Byte, NULL);

	cout << "MUL image X len: " << imgMULXlen << endl;
	cout << "MUL image Y len: " << imgMULYlen << endl;
	cout << "MUL image Band num : " << MULBandNum << endl;

	int imgPANXlen = srcPANDS->GetRasterXSize();
	int imgPANYlen = srcPANDS->GetRasterYSize();
	int PANBandNum = srcPANDS->GetRasterCount();
	cout << "PAN image X len: " << imgPANXlen << endl;
	cout << "PAN image Y len: " << imgPANYlen << endl;
	cout << "PAN image Band num : " << PANBandNum << endl;

	// ����ʣ�������涨�и߶ȵ��еĸ߶�
	int linesY = imgMULYlen / lineHeight;
	int restY = imgMULYlen - lineHeight * linesY;
	cout << "Line num:" << linesY << endl;
	cout << "Rest pixel:" << restY << endl;

	float *panBuffTmp;
	float **RGB;
	float **resultRGB;
	float **IHS;
	RGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	resultRGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	IHS = (float**)CPLMalloc(MULBandNum * sizeof(float));
	panBuffTmp = (float*)CPLMalloc(lineHeight*lineLen * sizeof(float));

	for (int k = 0; k < 3; k++)
	{
		RGB[k] = (float*)CPLMalloc(lineHeight*lineLen * sizeof(float));
		IHS[k] = (float*)CPLMalloc(lineHeight*lineLen * sizeof(float));
		resultRGB[k] = (float*)CPLMalloc(lineHeight*lineLen * sizeof(float));
	}
	int Y;
	for (int i = 0; i <= linesY; i++)
	{
		// �����е�ʵ�ʸ߶�
		if (i == linesY)
		{
			Y = restY;
		}
		else
		{
			Y = lineHeight;
		}
		//cout << "Working..." << endl;
		//cout << "Line " << i << "   Xoff:0  Yoff:" << i*lineHeight << endl;
		//cout << "Line size:" << Y << "*" << lineLen << "=" << lineLen*Y << endl;
		// ��ȡͼ������
		srcPANDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i*lineHeight, lineLen, Y, panBuffTmp, lineLen, Y, GDT_Float32, 0, 0);
		for (int j = 0; j < 3; j++)
		{
			srcMULDS->GetRasterBand(j + 1)->RasterIO(GF_Read, 0, i*lineHeight, lineLen, Y, RGB[j], lineLen, Y, GDT_Float32, 0, 0);
		}
		// �����㷨����
		IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, lineLen, Y);
		// ����д��ͼ��
		for (int k = 0; k < 3; k++)
		{
			resDS->GetRasterBand(k + 1)->RasterIO(GF_Write, 0, i*lineHeight, lineLen, Y, resultRGB[k], lineLen, Y, GDT_Float32, 0, 0);
		}
	}
	// �ͷ��ڴ�
	CPLFree(panBuffTmp);
	CPLFree(RGB);
	CPLFree(IHS);
	CPLFree(resultRGB);
	GDALClose(srcMULDS);
	GDALClose(srcPANDS);
	GDALClose(resDS);
}

int main()
{
	char* srcMULPath = "Data\\pictures\\Mul_large.tif";
	char* srcPANPath = "Data\\pictures\\Pan_large.tif";
	char* resPathByBlock= "Data\\output\\Big-image-fusion-block.tif";
	char* resPathByLine = "Data\\output\\Big-image-fusion-line.tif";
	int t1 = clock();
	IHSBigImageFusionByLine(srcMULPath, srcPANPath, resPathByLine);
	int t2 = clock();
	IHSBigImageFusionByBlock(srcMULPath, srcPANPath, resPathByBlock);
	int t3 = clock();

	cout << "Image fusion by line time:" << (t2 - t1) / CLOCKS_PER_SEC << endl;
	cout << "Image fusion by block time:" << (t3 - t2) / CLOCKS_PER_SEC << endl;
	system("pause");
	return 0;
}
