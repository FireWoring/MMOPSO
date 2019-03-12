#ifndef __MMOPSO_H__
#define __MMOPSO_H__

#include "canshu.h"
#include "Coordinate.h"
#include "Particle.h"
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <math.h>
#include <iterator>

#define MUTATIONNUM    5
#define NUM            20
#define PARTICLENUM    NUM    //粒子群粒子数
#define NODOMINATENUM  NUM
#define ObjectFuncNum  2    //两个目标函数 
#define DECISIONLENGTH 6    //决策变量的长度

class MMOPSO
{
private:
	static int Cmax, Cmin;
	static double Umax, Umin;
	static double AlphaHmax, AlphaHmin;
	static double AlphaCmax, AlphaCmin;
	static double BetaHmax, BetaHmin;
	static double BetaCmax, BetaCmin;
	static double Vcmax, Vcmin;
	static double Vumax, Vumin;
	static double VAlphaHmax, VAlphaHmin;
	static double VAlphaCmax, VAlphaCmin;
	static double VBetaHmax, VBetaHmin;
	static double VBetaCmax, VBetaCmin;

	static double pc, pm;//交叉率、突变率
	static double gc, gm;//SBX与PM的分布指数
	static double delta;
	static double theta; //分解方法中的参数
	static int    max_ev;//最大评估次数
public:
	vector<Particle *> pop;//粒子群
	vector<Particle *> rep;//非支配粒子集合
	vector<vector<double> > weight_vec;//权重向量
	vector<double> z;        //参考点位置
	vector<double> distanceAll; //总的距离度量
	vector<Particle *> pbest;
	Particle * gbest;
	vector<Particle *> E; //进化策略中的E，即ref的一半
	pair<Particle*, Particle*> C; //交叉之后生成的结果
	vector<Particle *> S; //进化算法生成的S集合

	double c1, c2;//学习因子
	double w;
public:
	MMOPSO(){};
	void InitWeightVector(); //初始化权重向量
	void InitZ(); //初始化z*
	void InitPopulation();//初始化粒子群的位置和速度
	double GetG(Particle *p, vector<double>& weight); //计算所以粒子的g值
	int CheckDominance(Particle*& x, Particle*& y); // 如果函数返回1，则表示x支配y;否则，当y占优或与x相等时，函数返回-1
	void UpdateArchive(vector<Particle *>& pop); //更新rep
	void CrowdingDistanceAssignment(); //计算rep中的拥挤距离值
	void DeleteMostCrowdedOne(); //删除最拥挤的一个粒子
	void SelectPbest(); //选取Xpbest
	void SelectGbest(); //选取Xgbest
	void UpdateVelocity(); //更新粒子速度
	void UpdatePosition(); // 更新粒子位置
	void UpdateReferencePoint(); //更新z*
	void EvolutionarySearchStrategy(); //粒子进化策略
	void GetHalfRef(); //获取一半ref，且保证包含具有更大拥挤距离值的非支配解
	void SBX(Particle *x,Particle *y); //模拟二进制交叉
	void MP(Particle *x); //变异操作
	void CompleteMMOPSO(); //完整的MMOPSO算法

	void SetW();
	void SetC();
};

struct SortObject
{
	int id;
	Particle *p;
	SortObject(int i, Particle *pt) :id(i), p(pt){}
};

struct SortDistance
{
	int id;
	double distance;
	SortDistance(int i, double dis) :id(i), distance(dis){}
};

#endif