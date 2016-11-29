﻿#ifndef NASH
#define NASH
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

double NashEE(CGraph *G,CGraph *GOR,vector<demand> & req,double OPEN){
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr constraint(env);
		for(int d = 0; d <  num; d++)
			constraint += req[d].flow*x[d][i];
		model.add( constraint <= G->Link[i]->capacity );  
	}

	//优化目标 节能 OPEN+X^2
	IloExpr cost(env);
	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		IloIntVarArray tmp(env,num,0,1);
		for(int d = 0; d < num; d++){
			load += req[d].flow*x[d][i];
			tmp[d] = x[d][i];
		}
		//cost += ( IloPower(load,2) + IloMax(tmp)*OPEN );
		cost += ( load*load + IloMax(tmp)*OPEN );
	}
	model.add(IloMinimize(env,cost));

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	if(EEsolver.solve()){
		obj = EEsolver.getObjValue(); //energy

		for(int i = 0;i < G->m; i++){  
			double loadc = 0;
			for(int d = 0;d < num; d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->bw = G->Link[i]->capacity -loadc; //residual bw
		}

		for(int m = 0;m < GOR->m; m++){	
			bool flag = false;
			double resbw = INF;
			for(int d = 0;d < num;d ++){
				if ( GOR->Link[m]->tail == req[d].org && GOR->Link[m]->head == req[d].des && (req[d].flow>0)){
					if(flag==false){
						flag=true;
						for(int i=0;i<G->m;i++)
							if(EEsolver.getValue(x[d][i])>0.5)
								resbw  = min( resbw , EEsolver.getValue(x[d][i])*G->Link[i]->bw);
					}
				}
				if(flag==true){
					GOR->Link[m]->bw = resbw ;	
					break;
				}
			} 
		} 		
	}
	else{
		cout << "nashEE unfeasible"<<endl;
	}
	for(size_t i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

#endif