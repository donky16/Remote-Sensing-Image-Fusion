// GDALDemo.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "gdal/gdal_priv.h"
#include<iostream>
#include<math.h>
#include<time.h>
#pragma comment(lib, "gdal_i.lib")

using namespace std;
// *****************************************分割线*************************************************
/*
	修改一部分图像得像素值
	RasterIO() read一列一列的存储到bufftmp中的
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
	// 修改图像block数据
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
	// 分配内存
	bufftmp = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));

	// 加载驱动
	poDstDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(dstPath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	for (int i = 0; i < bandNum; i++)
	{
		// 读取图像数据到内存
		poSrcDS->GetRasterBand(i + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, bufftmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		// 将内存中数据写入图像
		poDstDS->GetRasterBand(i + 1)->RasterIO(GF_Write, 0, 0, imgXlen, imgYlen, bufftmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		// 修改内存中部分数据
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
// *****************************************分割线*************************************************
/*
	图像合成函数，将superman图像得绿色背景扣掉，把superman嵌入space图像上
	方式：
		1.将两个图像得三个通道分别存储在supermanBufftmp[3]，spaceBufftmp[3]中
		2.比较superman中像素值，将绿背景去掉，其他像素值赋值给space图像中相同得位置
		3.将spaceBufftemp数据写入结果图像，生成合成图像
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

	GByte* supermanBufftmp[3];		// bandNum数为3
	GByte* spaceBufftmp[3];
	// 图片数据写入bufftmp
	for (int n = 0; n < bandNum; n++)
	{
		supermanBufftmp[n] = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
		supermanSrcDS->GetRasterBand(n + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, supermanBufftmp[n], imgXlen, imgYlen, GDT_Byte, 0, 0);
		spaceBufftmp[n] = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
		spaceSrcDS->GetRasterBand(n + 1)->RasterIO(GF_Read, 0, 0, imgXlen, imgYlen, spaceBufftmp[n], imgXlen, imgYlen, GDT_Byte, 0, 0);
	}
	// 将superman中除绿背景得像素写到space上
	int flag;
	for (int j = 0; j < imgYlen; j++)
	{
		for (int i = 0; i < imgXlen; i++)
		{
			flag = 1;
			RValue = (GByte)supermanBufftmp[0][i*imgYlen + j];
			GValue = (GByte)supermanBufftmp[1][i*imgYlen + j];
			BValue = (GByte)supermanBufftmp[2][i*imgYlen + j];
			// 将绿背景范围除掉
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
	// 加载驱动
	resultDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(supermanAndSpacePath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	// 写入合成图像并释放
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
// *****************************************分割线*************************************************
/*
	线性滤波实现
*/
GByte* getConvolutionSum(GByte* bufftemp, int Xlen, int Ylen, double convolutionKernel[], int kernelLen, int param, int offset)
{
	/*
		获取所有像素的卷积和的值，一个波段的像素值
		params: 源图像数据，原图像的X，Y，和卷积核，卷积核的长度(宽度)，卷积和需要除的数，和偏移量
		return: 经过卷积和计算的图像数据GByte类型的指针

		return 图像矩阵 * 矩阵
	*/
	GByte* sumData = (GByte*)CPLMalloc(Xlen*Ylen * sizeof(GByte));
	double tmp = 0;		// 存储每一个像素的卷积和
	for (int i = 0; i < Xlen; i++)
	{
		for (int j = 0; j < Ylen; j++)
		{
			tmp = 0;
			// 筛选出边界并置0
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
							卷积核中心坐标是[(kernelLen - 1) / 2, (kernelLen - 1) / 2]
							卷积核中任一坐标为(m,n)，则(XDiff,YDiff)为相对中心位置
							根据相对位置，就能算出图像中任一像素的邻域像素值
						*/
						int XDiff = n - (kernelLen - 1) / 2;
						int YDiff = m - (kernelLen - 1) / 2;
						tmp += convolutionKernel[m*kernelLen + n] * bufftemp[i*Ylen + j + XDiff*Ylen + YDiff];
					}
				}
				tmp = tmp / param + offset;
				// 截断，如果像素值不在（0,255）之间
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
		根据原图像，卷积核，卷积核的长度(宽度)，卷积和需要除的数，和偏移量得到转化之后的图像
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
		线性滤波主函数，根据不同卷积核获取不同的图像
	*/
	char* poSrcPath = "Data\\pictures\\lena.jpg";
	char* resultPath = "Data\\output\\convolution.jpg";
	int kernekLen = 0;

	/*  卷积核1  均值模糊
		 0 , 1 ,  0
		 1 , 1 ,  1		/ 5  （param=5）
		 0 , 1 ,  0
	*/
	resultPath = "Data\\output\\convolution-1.jpg";
	double convolutionKernel_1[9] = { 0, 1, 0, 1, 1, 1, 0, 1, 0 };
	kernekLen = sqrt(sizeof(convolutionKernel_1) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_1, kernekLen, 5, 0);

	/*  卷积核2  运动模糊
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

	/* 卷积核3 边缘检测
		-1,-1,-1,
		-1, 8,-1,
		-1,-1,-1,

	*/
	resultPath = "Data\\output\\convolution-3.jpg";
	double convolutionKernel_3[9] = { -1,-1,-1,-1,8,-1,-1,-1,-1 };
	kernekLen = sqrt(sizeof(convolutionKernel_3) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_3, kernekLen, 1, 0);

	/* 卷积核4  锐化
		-1,-1,-1,
		-1, 9,-1,
		-1,-1,-1,
	*/
	resultPath = "Data\\output\\convolution-4.jpg";
	double convolutionKernel_4[9] = { -1,-1,-1,-1,9,-1,-1,-1,-1 };
	kernekLen = sqrt(sizeof(convolutionKernel_4) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_4, kernekLen, 1, 0);
	/* 卷积核5  浮雕
		-1,-1, 0,
		-1, 0, 1,		offset=128
		 0, 1, 1,
	*/
	resultPath = "Data\\output\\convolution-5.jpg";
	double convolutionKernel_5[9] = { -1,-1,0,-1,0,1,0,1,1 };
	kernekLen = sqrt(sizeof(convolutionKernel_5) / sizeof(double));
	convolutionConvert(poSrcPath, resultPath, convolutionKernel_5, kernekLen, 1, 128);
	/* 卷积核6 高斯模糊
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
		IHS图像融合核心算法
		param: RGB-多光谱图像的RGB值 panBuffTmp-全色图像的像素值 IHS-函数外定义的IHS 
		       resultRGB-存储融合后的图像像素值 imgXlen,imgYlen-图像的长度和宽度	
	*/
	float matrix[3][3] = { { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 },{ -sqrt(2.0) / 6.0, -sqrt(2.0) / 6.0, sqrt(2.0) / 3.0 },{ 1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0 } };
	float matrixReverse[3][3] = { { 1.0, -1.0 / sqrt(2.0), 1.0 / sqrt(2.0) },{ 1.0, -1.0 / sqrt(2.0), -1.0 / sqrt(2.0) },{ 1.0, sqrt(2.0), 0 } };
	// RGB2IHS
	for (int m = 0; m < 3; m++)
	{
		for (int n = 0; n < imgXlen*imgYlen; n++)
		{
			IHS[m][n] = 0;	// 初始化为0，函数需要多次调用
			for (int k = 0; k < 3; k++)
			{
				IHS[m][n] += matrix[m][k] * RGB[k][n];
			}
		}
	}
	// 将pan图像数据置换I
	// IHS[0] = panBuffTmp; 指针赋值会导致IHS[0]指向panBuffTmp，下一次调用此函数就会赋值失败
	for (int i = 0; i < imgXlen*imgYlen; i++)
	{
		IHS[0][i] = panBuffTmp[i];
	}
	// 逆变换
	for (int mm = 0; mm < 3; mm++)
	{
		for (int nn = 0; nn < imgXlen*imgYlen; nn++)
		{
			resultRGB[mm][nn] = 0;	// 初始化为0，函数需要多次调用
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
		遥感图像融合主函数
		param: 多光谱图像路径，全色图像路径，融合结果图像路径
	*/
	GDALAllRegister();

	GDALDataset* srcMULDS = (GDALDataset*)GDALOpenShared(srcMULPath, GA_ReadOnly);
	GDALDataset* srcPANDS = (GDALDataset*)GDALOpenShared(srcPANPath, GA_ReadOnly);

	int imgMULXlen = srcMULDS->GetRasterXSize();
	int imgMULYlen = srcMULDS->GetRasterYSize();
	int MULBandNum = srcMULDS->GetRasterCount();
	// 创建驱动生成图像，可是发现直接生成GDT_Float32的不行，很慌
	GDALDataset *resDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(resultPath, imgMULXlen, imgMULYlen, MULBandNum, GDT_Byte, NULL);
	
	cout << "MUL image X len: " << imgMULXlen << endl;
	cout << "MUL image Y len: " << imgMULYlen << endl;
	cout << "MUL image Band num : " << MULBandNum << endl;

	int imgPANXlen = srcPANDS->GetRasterXSize();
	int imgPANYlen = srcPANDS->GetRasterYSize();
	int PANBandNum = srcPANDS->GetRasterCount();
	// 用float32类型，原图像每一个像素点并不是整数
	float* panBuffTmp = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
	srcPANDS->GetRasterBand(PANBandNum)->RasterIO(GF_Read, 0, 0, imgMULXlen, imgMULYlen, panBuffTmp, imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);
	
	cout << "PAN image X len: " << imgPANXlen << endl;
	cout << "PAN image Y len: " << imgPANYlen << endl;
	cout << "PAN image Band num : " << PANBandNum << endl;
	
	// 二维数组，每一个开创数组都是用CPLMalloc，中间的IHS也在此malloc，这样省malloc和free的时间
	float **RGB;
	float **resultRGB;
	float **IHS;
	RGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	resultRGB = (float**)CPLMalloc(MULBandNum * sizeof(float));
	IHS = (float**)CPLMalloc(MULBandNum * sizeof(float));
	
	// 读取数据到数组中
	for (int i = 0; i < 3; i++)
	{
		RGB[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
		srcMULDS->GetRasterBand(i + 1)->RasterIO(GF_Read, 0, 0, imgMULXlen, imgMULYlen, RGB[i], imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);	
		resultRGB[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
		IHS[i] = (float*)CPLMalloc(imgMULXlen*imgMULYlen * sizeof(float));
	}

	// 核心算法直接计算，结果存储在resultRGB中
	IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, imgMULXlen, imgMULYlen);

	// 直接将float类型像素值放进去，这里创建的是GDT_Byte的，这里竟然可以直接放进去，迷
	for (int j = 0; j < 3; j++)
	{
		resDS->GetRasterBand(j + 1)->RasterIO(GF_Write, 0, 0, imgMULXlen, imgMULYlen, resultRGB[j], imgMULXlen, imgMULYlen, GDT_Float32, 0, 0);
	}

	// 释放内存
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
		分块图像融合
		param: 多光谱图像路径，全色图像路径，融合结果图像路径
	*/
	GDALAllRegister();
	// 块的大小
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

	// 计算出右边剩余的像素和下面剩余的像素
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
			// 计算块的实际大小
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
			
			// 根据块的大小来从图像中读取数据
			srcPANDS->GetRasterBand(1)->RasterIO(GF_Read, n*blockLen, m*blockLen, X, Y, panBuffTmp, X, Y, GDT_Float32, 0, 0);
			for (int j = 0; j < 3; j++)
			{		
				srcMULDS->GetRasterBand(j + 1)->RasterIO(GF_Read, n*blockLen, m*blockLen, X, Y, RGB[j], X, Y, GDT_Float32, 0, 0);
			}

			// 核心算法融合
			IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, X, Y);

			// 数据写入图像
			for (int i = 0; i < 3; i++)  
			{
				resDS->GetRasterBand(i + 1)->RasterIO(GF_Write, n*blockLen, m*blockLen, X, Y, resultRGB[i], X, Y, GDT_Float32, 0, 0);
			}
		}
	}
	// 释放内存
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
		分行图像融合
		param: 多光谱图像路径，全色图像路径，融合结果图像路径
	*/
	GDALAllRegister();

	GDALDataset* srcMULDS = (GDALDataset*)GDALOpenShared(srcMULPath, GA_ReadOnly);
	GDALDataset* srcPANDS = (GDALDataset*)GDALOpenShared(srcPANPath, GA_ReadOnly);

	int imgMULXlen = srcMULDS->GetRasterXSize();
	int imgMULYlen = srcMULDS->GetRasterYSize();
	int MULBandNum = srcMULDS->GetRasterCount();
	int lineHeight = 256;		// 行高度
	int lineLen = imgMULXlen;	// 行长度

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

	// 计算剩余最后不足规定行高度的行的高度
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
		// 计算行的实际高度
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
		// 读取图像数据
		srcPANDS->GetRasterBand(1)->RasterIO(GF_Read, 0, i*lineHeight, lineLen, Y, panBuffTmp, lineLen, Y, GDT_Float32, 0, 0);
		for (int j = 0; j < 3; j++)
		{
			srcMULDS->GetRasterBand(j + 1)->RasterIO(GF_Read, 0, i*lineHeight, lineLen, Y, RGB[j], lineLen, Y, GDT_Float32, 0, 0);
		}
		// 核心算法计算
		IHSImageFusion(RGB, panBuffTmp, IHS, resultRGB, lineLen, Y);
		// 数据写入图像
		for (int k = 0; k < 3; k++)
		{
			resDS->GetRasterBand(k + 1)->RasterIO(GF_Write, 0, i*lineHeight, lineLen, Y, resultRGB[k], lineLen, Y, GDT_Float32, 0, 0);
		}
	}
	// 释放内存
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
