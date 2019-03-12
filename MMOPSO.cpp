#define _CRT_SECURE_NO_WARNINGS
#include "MMOPDO.h"

int MMOPSO::Cmin = 1;
int MMOPSO::Cmax = 15;
double MMOPSO::Umin = 2.0;
double MMOPSO::Umax = 3.0;
double MMOPSO::AlphaHmin = 3.0;
double MMOPSO::AlphaHmax = 4.5;
double MMOPSO::AlphaCmin = 1.0;
double MMOPSO::AlphaCmax = 3.0;
double MMOPSO::BetaHmin = 5.0;
double MMOPSO::BetaHmax = 8.0;
double MMOPSO::BetaCmin = 3.0;
double MMOPSO::BetaCmax = 4.5;
double MMOPSO::Vcmax = MMOPSO::Cmax - MMOPSO::Cmin;
double MMOPSO::Vcmin = 0 - MMOPSO::Vcmax;
double MMOPSO::Vumax = MMOPSO::Umax - MMOPSO::Umin;
double MMOPSO::Vumin = 0 - MMOPSO::Vumax;
double MMOPSO::VAlphaHmax = MMOPSO::AlphaHmax - MMOPSO::AlphaHmin;
double MMOPSO::VAlphaHmin = 0 - MMOPSO::VAlphaHmax;
double MMOPSO::VAlphaCmax = MMOPSO::AlphaCmax - MMOPSO::AlphaCmin;
double MMOPSO::VAlphaCmin = 0 - MMOPSO::VAlphaCmax;
double MMOPSO::VBetaHmax = MMOPSO::BetaHmax - MMOPSO::BetaHmin;
double MMOPSO::VBetaHmin = 0 - MMOPSO::VBetaHmax;
double MMOPSO::VBetaCmax = MMOPSO::BetaCmax - MMOPSO::BetaCmin;
double MMOPSO::VBetaCmin = 0 - MMOPSO::VBetaCmax;

double MMOPSO::pc = 0.9;
double MMOPSO::pm = 1.0 / NUM;
double MMOPSO::gc = 20;
double MMOPSO::gm = 20;
double MMOPSO::delta = 0.9;
double MMOPSO::theta = 0.5;
int    MMOPSO::max_ev = 10;

void MMOPSO::InitWeightVector()
{
	double weight;
	for (int i = 0; i < NUM; i++)
	{
		weight = rand() / (double)RAND_MAX;
		vector<double> tmp = { weight, 1 - weight };
		weight_vec.push_back(tmp);
	}
}

void MMOPSO::InitZ()
{
	z.push_back(pop[0]->f->F(pop[0]->answer));
	z.push_back(pop[0]->f->TNom(pop[0]->answer));
	for (int i = 1; i < NUM; i++)
	{
		if (pop[i]->f->F(pop[i]->answer) < z[0])
			z[0] = pop[i]->f->F(pop[i]->answer);
		if (pop[i]->f->TNom(pop[i]->answer) < z[1])
			z[1] = pop[i]->f->TNom(pop[i]->answer);
	}
}

void MMOPSO::InitPopulation()
{
	int randC;
	double randU, randAlphaH, randAlphaC, randBetaH, randBetaC;
	srand((unsigned)time(NULL));
	for (int i = 0; i < PARTICLENUM; i++)
	{
		randC = static_cast<int>((rand() / (double)RAND_MAX)*(Cmax - Cmin)) + Cmin;
		randU = ((rand() / (double)RAND_MAX)*(Umax - Umin)) + Umin;
		randAlphaH = ((rand() / (double)RAND_MAX)*(AlphaHmax - AlphaHmin)) + AlphaHmin;
		randAlphaC = ((rand() / (double)RAND_MAX)*(AlphaCmax - AlphaCmin)) + AlphaCmin;
		randBetaH = ((rand() / (double)RAND_MAX)*(BetaHmax - BetaHmin)) + BetaHmin;
		randBetaC = ((rand() / (double)RAND_MAX)*(BetaCmax - BetaCmin)) + BetaCmin;

		pop.push_back(new Particle(randC, randU, randAlphaH, randAlphaC, randBetaH, randBetaC));

		/*cout << "c:" << pop[i]->getCorrdinate().c << " u:" << pop[i]->getCorrdinate().u <<
			" AlphaH:" << pop[i]->getCorrdinate().alphaH << " AlphaC:" << pop[i]->getCorrdinate().alphaC <<
			" BetaH:" << pop[i]->getCorrdinate().betaH << " BetaC:" << pop[i]->getCorrdinate().betaC << endl;*/
		cout << "F:" << pop[i]->f->F(pop[i]->answer) << " Tnom:" << pop[i]->f->TNom(pop[i]->answer) << endl;
	}

	//初始化pBest数组
	pbest.assign(NNUM, 0);
}

double MMOPSO::GetG(Particle *p,vector<double>& weight)
{
	vector<double> F;
	F.push_back(p->f->F(p->answer));
	F.push_back(p->f->TNom(p->answer));

	vector<double> tmp;
	tmp.push_back((z[0] - F[0])*weight[0]);
	tmp.push_back((z[1] - F[1])*weight[1]);
	double d1 = sqrt(pow(tmp[0], 2) + pow(tmp[1], 2)) / sqrt(pow(weight[0], 2) + pow(weight[1], 2));

	tmp[0] = F[0] - (z[0] - d1*weight[0]);
	tmp[1] = F[1] - (z[1] - d1*weight[1]);
	double d2 = sqrt(pow(tmp[0], 2) + pow(tmp[1], 2));

	return (d1 + theta * d2);
}

int MMOPSO::CheckDominance(Particle*& x, Particle*& y)
{
	//x支配y返回true，否则返回false
	if (x->f->F(x->answer) <= y->f->F(y->answer) && x->f->TNom(x->answer) < y->f->TNom(y->answer))
		return 1;
	else if (x->f->F(x->answer) < y->f->F(y->answer) && x->f->TNom(x->answer) <= y->f->TNom(y->answer))
		return 1;
	else if (x->f->F(x->answer) >= y->f->F(y->answer) && x->f->TNom(x->answer) > y->f->TNom(y->answer))
		return -1;
	else if (x->f->F(x->answer) > y->f->F(y->answer) && x->f->TNom(x->answer) >= y->f->TNom(y->answer))
		return -1;
	else
		return 0;
}

void MMOPSO::UpdateArchive(vector<Particle *>& popSet)
{
	int state = 1;
	//vector<int> 
	for (int i = 0; i < popSet.size(); i++)
	{
		vector<int> del;
		for (int j = 0; j < rep.size(); j++)
		{
			state = CheckDominance(popSet[i], rep[j]);
			if (state == 1)
			{
				del.push_back(j);
			}
			else if (state == -1)
			{
				break;
			}
		}
		for (int j = del.size()-1; j >=0; j--)
		{
			Particle *tmp = rep[del[j]];
			rep.erase(rep.begin() + del[j]);
			delete tmp;
		}
		if (rep.size() == 0 || state != -1)
		{
			//检测该粒子是否已经在rep当中
			int has = false;
			for (int k = 0; k < rep.size(); k++)
			{
				if (popSet[i]->getCorrdinate() == rep[k]->getCorrdinate())
					has = true;
				if (has == true)
					break;
			}
			//如果不在rep当中，就加入进去
			if (has == false)
			{
				Particle* newOne = new Particle(popSet[i]->getCorrdinate().c, popSet[i]->getCorrdinate().u, popSet[i]->getCorrdinate().alphaH, popSet[i]->getCorrdinate().alphaC, popSet[i]->getCorrdinate().betaH, popSet[i]->getCorrdinate().betaC);
				rep.push_back(newOne);
				if (rep.size() > NUM)
				{
					CrowdingDistanceAssignment();
					DeleteMostCrowdedOne();
				}
			}
		}


		cout << endl;
		for (int i = 0; i < rep.size(); i++)
		{
			/*cout << "c:" << rep[i]->getCorrdinate().c << " u:" << rep[i]->getCorrdinate().u <<
			" AlphaH:" << rep[i]->getCorrdinate().alphaH << " AlphaC:" << rep[i]->getCorrdinate().alphaC <<
			" BetaH:" << rep[i]->getCorrdinate().betaH << " BetaC:" << rep[i]->getCorrdinate().betaC << endl;*/
			cout << "F:" << rep[i]->f->F(rep[i]->answer) << " Tnom:" << rep[i]->f->TNom(rep[i]->answer) << endl;
		}
	}
}

void MMOPSO::CrowdingDistanceAssignment()
{
	vector<SortObject> sortObject;
	vector<double> distance(rep.size(), -1);
	for (int i = 0; i < rep.size(); i++)
	{
		sortObject.push_back(SortObject(i, rep[i]));
	}

	//求第一个函数的距离度量
	sort(sortObject.begin(), sortObject.end(),
		[&](SortObject& x, SortObject& y)
		{
			return x.p->f->F(x.p->answer) < y.p->f->F(y.p->answer);
		});
	double max_distance = sortObject[rep.size()-1].p->f->F(sortObject[rep.size()-1].p->answer) - sortObject[0].p->f->F(sortObject[0].p->answer);
	for (int i = 1; i < rep.size() - 1; i++)
	{
		distance[sortObject[i].id] = (sortObject[i + 1].p->f->F(sortObject[i + 1].p->answer) - sortObject[i - 1].p->f->F(sortObject[i - 1].p->answer))/max_distance;
	}

	//求第二个函数的距离度量
	sort(sortObject.begin(), sortObject.end(),
		[&](SortObject& x, SortObject& y)
	{
		return x.p->f->TNom(x.p->answer) < y.p->f->TNom(y.p->answer);
	});
	max_distance = sortObject[rep.size() - 1].p->f->F(sortObject[rep.size() - 1].p->answer) - sortObject[0].p->f->F(sortObject[0].p->answer);
	distance[sortObject[rep.size() - 1].id] = -1;
	distance[sortObject[0].id] = -1;
	for (int i = 1; i < rep.size() - 1 && distance[sortObject[0].id]!= -1; i++)
	{
		distance[sortObject[i].id] += ((sortObject[i + 1].p->f->F(sortObject[i + 1].p->answer) - sortObject[i - 1].p->f->F(sortObject[i - 1].p->answer))/max_distance);
	}

	//设置总的距离度量
	distanceAll.clear();
	distanceAll.assign(distance.begin(), distance.end());
}

void MMOPSO::DeleteMostCrowdedOne()
{
	double min_distance = DBL_MAX;
	int index = 0, min_index = 0;
	for_each(distanceAll.begin(), distanceAll.end(), 
		[&](double& x)
		{ 
			if (x > 0 && x < min_distance)
			{
				min_distance = x;
				min_index = index;
			}
			index++;
		});

	Particle *tmp = rep[min_index];
	rep.erase(rep.begin() + min_index);
	delete tmp;
}

void MMOPSO::SelectPbest()
{
	for (int i = 0; i < NUM; i++)
	{
		if (!rep.empty())
			pbest[i] = rep[0];
		else
			break;

		for (int j = 1; j < rep.size(); j++)
		{
			if (GetG(pbest[i], weight_vec[i])>GetG(rep[j], weight_vec[i]))
				pbest[i] = rep[j];
		}
	}
}

void MMOPSO::SelectGbest()
{
	int index = static_cast<int>((double)rand() / RAND_MAX*rep.size());
	Particle tmp(*(rep[index]));
	gbest = &tmp;
}

void MMOPSO::UpdateVelocity()
{
	double r = (double)rand() / RAND_MAX;
	if (r < delta)
	{
		for (int i = 0; i < pop.size(); i++)
		{
			pop[i]->Vc = static_cast<int>(w*pop[i]->Vc + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().c - pop[i]->getCorrdinate().c));
			pop[i]->Vu = w*pop[i]->Vu + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().u - pop[i]->getCorrdinate().u);
			pop[i]->ValphaH = w*pop[i]->ValphaH + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().alphaH - pop[i]->getCorrdinate().alphaH);
			pop[i]->ValphaC = w*pop[i]->ValphaC + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().alphaC - pop[i]->getCorrdinate().alphaC);
			pop[i]->VbetaH = w*pop[i]->VbetaH + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().betaH - pop[i]->getCorrdinate().betaH);
			pop[i]->VbetaC = w*pop[i]->VbetaC + c1*(double)rand() / RAND_MAX*(pbest[i]->getCorrdinate().betaC - pop[i]->getCorrdinate().betaC);
		}
	}
	else
	{
		SelectGbest();
		for (int i = 0; i < pop.size(); i++)
		{
			pop[i]->Vc = static_cast<int>(w*pop[i]->Vc + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().c - pop[i]->getCorrdinate().c));
			pop[i]->Vu = w*pop[i]->Vu + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().u - pop[i]->getCorrdinate().u);
			pop[i]->ValphaH = w*pop[i]->ValphaH + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().alphaH - pop[i]->getCorrdinate().alphaH);
			pop[i]->ValphaC = w*pop[i]->ValphaC + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().alphaC - pop[i]->getCorrdinate().alphaC);
			pop[i]->VbetaH = w*pop[i]->VbetaH + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().betaH - pop[i]->getCorrdinate().betaH);
			pop[i]->VbetaC = w*pop[i]->VbetaC + c1*(double)rand() / RAND_MAX*(gbest->getCorrdinate().betaC - pop[i]->getCorrdinate().betaC);
		}
	}
}

void MMOPSO::UpdatePosition()
{
	for (int i = 0; i < pop.size(); i++)
	{
		int newC = min(max((int)(pop[i]->getCorrdinate().c + pop[i]->Vc), Cmin), Cmax);
		double newU = pop[i]->getCorrdinate().u + pop[i]->Vu > Umin ? pop[i]->getCorrdinate().u + pop[i]->Vu : Umin;
		newU = newU < Umax ? newU : Umax;
		double newAlphaH = pop[i]->getCorrdinate().alphaH + pop[i]->ValphaH > AlphaHmin ? pop[i]->getCorrdinate().alphaH + pop[i]->ValphaH : AlphaHmin;
		newAlphaH = newAlphaH < AlphaHmax ? newAlphaH : AlphaHmax;
		double newAlphaC = pop[i]->getCorrdinate().alphaC + pop[i]->ValphaC > AlphaCmin ? pop[i]->getCorrdinate().alphaC + pop[i]->ValphaC : AlphaCmin;
		newAlphaC = newAlphaC < AlphaCmax ? newAlphaC : AlphaCmax;
		double newBetaH = pop[i]->getCorrdinate().betaH + pop[i]->VbetaH > BetaHmin ? pop[i]->getCorrdinate().betaH + pop[i]->VbetaH : BetaHmin;
		newBetaH = newBetaH < BetaHmax ? newBetaH : BetaHmax;
		double newBetaC = pop[i]->getCorrdinate().betaC + pop[i]->VbetaC > BetaCmin ? pop[i]->getCorrdinate().betaC + pop[i]->VbetaC : BetaCmin;
		newBetaC = newBetaC < BetaCmax ? newBetaC : BetaCmax;
		pop[i]->setCoordinate(newC, newU, newAlphaH, newAlphaC, newBetaH, newBetaC);
	}
}

void MMOPSO::UpdateReferencePoint()
{
	for (int i = 0; i < NUM; i++)
	{
		if (pop[i]->f->F(pop[i]->answer) < z[0])
			z[0] = pop[i]->f->F(pop[i]->answer);
		if (pop[i]->f->TNom(pop[i]->answer) < z[1])
			z[1] = pop[i]->f->TNom(pop[i]->answer);
	}
}

void MMOPSO::EvolutionarySearchStrategy()
{
	GetHalfRef();
	S.clear();

	for (int i = 0; i < rep.size() && E.size() > 0; i++)
	{
		int j = static_cast<int>(rand() / (double)RAND_MAX*E.size());
		if (rand() / (double)RAND_MAX < pc)
		{
			SBX(rep[i], E[j]);
			if (rand() / (double)RAND_MAX < pm)
			{
				int k = rand() / (double)RAND_MAX < 0.5 ? 1 : 2;
				if (k == 1)
				{
					MP(C.first);
				}
				else
				{
					MP(C.second);
				}
			}
			delete C.first;
			delete C.second;
		}
	}

	E.clear();
}

void MMOPSO::GetHalfRef()
{
	//获取拥挤距离值
	CrowdingDistanceAssignment();

	//构造拥挤距离值和id之间相互对应的vector
	vector<SortDistance> distance_vec;
	int id = 0;
	for_each(distanceAll.begin(),distanceAll.end(),
		[&](double x){
			SortDistance tmp(id, x);
			distance_vec.push_back(tmp);
			++id;
		});

	//根据拥挤距离值排序
	sort(distance_vec.begin(), distance_vec.end(),
		[&](SortDistance& x, SortDistance& y){
			return x.distance < y.distance;
		});

	//构造E数组
	for (int i = 0; i < rep.size() / 2 &&distance_vec.size()>0; i++)
	{
		E.push_back(rep[distance_vec[i].id]);
	}
}

void MMOPSO::SBX(Particle *x, Particle *y)
{
	vector<double> betaRetio;
	for (int i = 0; i < DECISIONLENGTH; i++)
	{
		double u = rand() / (double)RAND_MAX;
		if (u <= 0.5)
		{
			betaRetio.push_back(pow((2 * u), 1.0 / (gc + 1)));
		}
		else
		{
			betaRetio.push_back(pow(1.0 / (2.0 * (1 - u)), 1.0 / (gc + 1)));
		}
	}	

	//计算交叉后的新个体
	int newC1 = static_cast<int>(0.5*((1 + betaRetio[0])*x->getCorrdinate().c + (1 - betaRetio[0])*y->getCorrdinate().c));
	int newC2 = static_cast<int>(0.5*((1 - betaRetio[0])*x->getCorrdinate().c + (1 + betaRetio[0])*y->getCorrdinate().c));
	double newU1 = 0.5*((1 + betaRetio[1])*x->getCorrdinate().u + (1 - betaRetio[1])*y->getCorrdinate().u);
	double newU2 = 0.5*((1 - betaRetio[1])*x->getCorrdinate().u + (1 + betaRetio[1])*y->getCorrdinate().u);
	double newAlphaH1 = 0.5*((1 + betaRetio[2])*x->getCorrdinate().alphaH + (1 - betaRetio[2])*y->getCorrdinate().alphaH);
	double newAlphaH2 = 0.5*((1 - betaRetio[2])*x->getCorrdinate().alphaH + (1 + betaRetio[2])*y->getCorrdinate().alphaH);
	double newAlphaC1 = 0.5*((1 + betaRetio[3])*x->getCorrdinate().alphaC + (1 - betaRetio[3])*y->getCorrdinate().alphaC);
	double newAlphaC2 = 0.5*((1 - betaRetio[3])*x->getCorrdinate().alphaC + (1 + betaRetio[3])*y->getCorrdinate().alphaC);
	double newBetaH1 = 0.5*((1 + betaRetio[4])*x->getCorrdinate().betaH + (1 - betaRetio[4])*y->getCorrdinate().betaH);
	double newBetaH2 = 0.5*((1 - betaRetio[4])*x->getCorrdinate().betaH + (1 + betaRetio[4])*y->getCorrdinate().betaH);
	double newBetaC1 = 0.5*((1 + betaRetio[5])*x->getCorrdinate().betaC + (1 - betaRetio[4])*y->getCorrdinate().betaC);
	double newBetaC2 = 0.5*((1 - betaRetio[5])*x->getCorrdinate().betaC + (1 + betaRetio[4])*y->getCorrdinate().betaC);

	C.first = new Particle(newC1, newU1, newAlphaH1, newAlphaC1, newBetaH1, newBetaC1);
	C.second = new Particle(newC2, newU2, newAlphaH2, newAlphaC2, newBetaH2, newBetaC2);
}

void MMOPSO::MP(Particle *x)
{
	vector<double> deltaRetio;
	vector<double> deltaMax{ double(Cmax), Umax, AlphaHmax, AlphaCmax, BetaHmax, BetaCmax };
	vector<double> deltaMin{ double(Cmin), Umin, AlphaHmin, AlphaCmin, BetaHmin, BetaCmin };
	vector<double> xi{ double(x->getCorrdinate().c), x->getCorrdinate().u, x->getCorrdinate().alphaH, x->getCorrdinate().alphaC, x->getCorrdinate().betaH, x->getCorrdinate().betaC };
	for (int i = 0; i < DECISIONLENGTH; i++)
	{
		double u = rand() / (double)RAND_MAX;
		double tmp = deltaMax[i] - xi[i]>xi[i] - deltaMin[i] ? deltaMax[i] - xi[i] : xi[i] - deltaMin[i];
		if (u < 0.5)
		{
			deltaRetio.push_back(pow((2 * u +(1-2*u)*pow(tmp/(deltaMax[i]-deltaMin[i]),gm+1)), 1.0 / (gm + 1))-1);
		}
		else
		{
			deltaRetio.push_back(1 - pow(2 * (1 - u) + 2 * (u - 0.5)*pow(tmp / (deltaMax[i] - deltaMin[i]), gm + 1), 1.0 / (gm + 1)));
		}
		xi[i] += deltaRetio[i]*(deltaMax[i] - deltaMin[i]);
	}
	S.push_back(new Particle(static_cast<int>(xi[0]), xi[1], xi[2], xi[3], xi[4], xi[5]));
}

void MMOPSO::SetW()
{
	w = (rand() / (double)RAND_MAX)*0.4 + 0.1;
}

void MMOPSO::SetC()
{
	c1 = (rand() / (double)RAND_MAX)*0.5 + 1.5;
	c2 = (rand() / (double)RAND_MAX)*0.5 + 1.5;
}

void MMOPSO::CompleteMMOPSO()
{
	cout << "step 1:" << endl;
	//初始化权重向量
	InitWeightVector();

	//初始化种群
	InitPopulation();

	//初始化z值，即参考点
	InitZ();

	//初始化外部档案
	UpdateArchive(pop);

	int ev = 0;
	while (ev < max_ev)
	{
		cout << "step 2:" << endl;
		cout << "ev:" << ev << endl;
		//选择xpbest
		SelectPbest();

		for (int i = 0; i < MUTATIONNUM; i++)
		{
			cout << "step 3:" << endl;
			//更新速度
			UpdateVelocity();

			//更新位置
			UpdatePosition();
			
			//更新参考点
			UpdateReferencePoint();
		}

		/*for (int i = 0; i < pop.size(); i++)
		cout << "F:" << pop[i]->f->F(pop[i]->answer) << " Tnom:" << pop[i]->f->TNom(pop[i]->answer) << endl;*/

		//更新外部存档
		UpdateArchive(pop);

		//实行进化策略
		EvolutionarySearchStrategy();

		//更新参考点
		UpdateReferencePoint();

		cout << "end" << endl;
		//更新外部存档
		UpdateArchive(S);

		//将得到的粒子写入文件
		char str[10] = { 0 };
		_itoa(ev, str, 10);
		char filename[20] = "./repx_";
		strcat(filename, str);
		strcat(filename, ".txt");
		ofstream out(filename);
		for (int k = 0; k < rep.size(); k++)
		{
			out << rep[k]->f->F(rep[k]->answer) << " " << rep[k]->f->TNom(rep[k]->answer) << endl;
		}
		out.close();

		ev++;
	}
}

void mainxxx()
{
	MMOPSO pso;
	for (int i = 0; i < 100; i++)
		cout << rand() / (double)RAND_MAX << endl;
	//pso.InitWeightVector();
	vector<double> deltaMaxAndMin{1,2,3,4,5};
	for_each(deltaMaxAndMin.begin(), deltaMaxAndMin.end(), [](double& x){cout << x << endl; });
	cin.get();
}