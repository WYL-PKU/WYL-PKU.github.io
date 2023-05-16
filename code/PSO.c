/*cern ROOT*/
#include "TRandom3.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THStack.h"
#include<iostream>
#include<fstream>
using namespace std;

const int N = 10000;/*粒子群规模*/
const int D = 2;/*参数空间维度*/
const int loop = 100;/*迭代次数*/
const double wmax = 0.8;/*惯性指数最大值*/
const double wmin = 0.2;/*惯性指数最小值*/
double w = wmax;/*惯性指数*/
const double c1 = 1.6;
const double c2 = 1.8;
double r1, r2;
int k = 0;

/*h和d数据*/
double hs[5] = { 1000,828,800,600,300 };
double ds[5] = { 1500,1340,1328,1172,800 };
double sigma = 15;

double x0[N][D] = { 0 }; double x1[N][D] = { 0 };
double v0[N][D] = { 0 }; double v1[N][D] = { 0 };
double pp[N][D] = { 0 }; double pg[D] = { 0 };
double fp[N]; double fg;

TRandom3* rng = new TRandom3(0);

double chi2(double alpha, double beta) {/*目标函数，使之最小化*/
	double chisq = 0;
	for (int i = 0; i < 5; i++) {
		chisq += TMath::Power(ds[i] - alpha * TMath::Power(hs[i], beta), 2) / TMath::Power(sigma, 2);
	}
	return chisq;
}

void initializing() {/*初始化坐标速度和拟合值*/
	double temp = 0;
	for (int i = 0; i < N; i++) {
		temp = rng->Rndm(); x1[i][0] = temp * 20+40;
		temp = rng->Rndm(); v1[i][0] = temp * 10-5;
		temp = rng->Rndm(); x1[i][1] = temp;
		temp = rng->Rndm(); v1[i][1] = temp * 0.4 - 0.2;
		fp[i] = 3000;
	}
	fg = 3000;
	return;
}

void swapping() {/*将上次迭代得到的坐标放入x0，准备下一次更新*/
	double temp;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < D; j++) {
			x0[i][j] = x1[i][j];
			v0[i][j] = v1[i][j];
		}
	}
	return;
}

void recursing() {/*进行一次迭代*/
	swapping();
	double temp;
	double rounding = 100000;
	/*计算拟合值，更新位置优劣信息*/
	for (int i = 0; i < N; i++) {
		temp = chi2(x1[i][0], x1[i][1]);
		rounding = temp < rounding ? temp : rounding;
		if (temp < fp[i]) {
			fp[i] = temp;
			pp[i][0] = x1[i][0];
			pp[i][1] = x1[i][1];
		}
		if (temp < fg) {
			fg = temp;
			pg[0] = x1[i][0];
			pg[1] = x1[i][1];
		}
	}
	cout << "this round best is " << rounding << endl;
	w = wmax - (wmax - wmin) * k / loop;/*更新惯性权重*/
	/*更新坐标和速度*/
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < D; j++) {
			r1 = rng->Rndm(); r2 = rng->Rndm();
			v1[i][j] = w * v0[i][j] + c1 * r1 * (pp[i][j] - x0[i][j]) + c2 * r2 * (pg[j] - x0[i][j]);
			x1[i][j] = x0[i][j] + v1[i][j];
		}
	}
	cout << "chi2 = " << fg << "\talpha = " << pg[0] << "\tbeta = " << pg[1] << endl;
	return;
}

void particalswarm() {
	while (k < loop) {
		recursing();
		k++;
	}
	return;
}

void PSO() {
	cout << "Hello, world!" << endl;
	initializing();
	particalswarm();
	return;
}