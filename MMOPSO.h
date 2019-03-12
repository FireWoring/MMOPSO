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
#define PARTICLENUM    NUM    //����Ⱥ������
#define NODOMINATENUM  NUM
#define ObjectFuncNum  2    //����Ŀ�꺯�� 
#define DECISIONLENGTH 6    //���߱����ĳ���

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

	static double pc, pm;//�����ʡ�ͻ����
	static double gc, gm;//SBX��PM�ķֲ�ָ��
	static double delta;
	static double theta; //�ֽⷽ���еĲ���
	static int    max_ev;//�����������
public:
	vector<Particle *> pop;//����Ⱥ
	vector<Particle *> rep;//��֧�����Ӽ���
	vector<vector<double> > weight_vec;//Ȩ������
	vector<double> z;        //�ο���λ��
	vector<double> distanceAll; //�ܵľ������
	vector<Particle *> pbest;
	Particle * gbest;
	vector<Particle *> E; //���������е�E����ref��һ��
	pair<Particle*, Particle*> C; //����֮�����ɵĽ��
	vector<Particle *> S; //�����㷨���ɵ�S����

	double c1, c2;//ѧϰ����
	double w;
public:
	MMOPSO(){};
	void InitWeightVector(); //��ʼ��Ȩ������
	void InitZ(); //��ʼ��z*
	void InitPopulation();//��ʼ������Ⱥ��λ�ú��ٶ�
	double GetG(Particle *p, vector<double>& weight); //�����������ӵ�gֵ
	int CheckDominance(Particle*& x, Particle*& y); // �����������1�����ʾx֧��y;���򣬵�yռ�Ż���x���ʱ����������-1
	void UpdateArchive(vector<Particle *>& pop); //����rep
	void CrowdingDistanceAssignment(); //����rep�е�ӵ������ֵ
	void DeleteMostCrowdedOne(); //ɾ����ӵ����һ������
	void SelectPbest(); //ѡȡXpbest
	void SelectGbest(); //ѡȡXgbest
	void UpdateVelocity(); //���������ٶ�
	void UpdatePosition(); // ��������λ��
	void UpdateReferencePoint(); //����z*
	void EvolutionarySearchStrategy(); //���ӽ�������
	void GetHalfRef(); //��ȡһ��ref���ұ�֤�������и���ӵ������ֵ�ķ�֧���
	void SBX(Particle *x,Particle *y); //ģ������ƽ���
	void MP(Particle *x); //�������
	void CompleteMMOPSO(); //������MMOPSO�㷨

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