#ifndef EE_H
#define EE_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

double LBdictor(CGraph *G,vector<demand> & req,int ornum){
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	 //优化目标
	IloNumVar z(env, 0, 1);	
    model.add(IloMinimize(env, z));

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
		IloExpr load(env);
		for(int d = 0; d <  num; d++)
			load += req[d].flow*x[d][i];
		model.add( load <= G->Link[i]->capacity );  
		model.add( load <= z*G->Link[i]->capacity);	
	}

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	if(EEsolver.solve()){
		obj = EEsolver.getObjValue(); //mlu
		//thoughtput
		double output = 0;
		for(int i = 0;i < G->m; i++){  
			double loadc = 0;
			for(int d = 0;d < num; d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->bw = G->Link[i]->capacity -loadc; //residual bw
		}
		for(int d = 0;d < ornum; d++){
			double res = INF;
			for(int i = 0;i < G->m; i++){
				if(EEsolver.getValue(x[d][i]) > 0.5)
					res = min( res, G->Link[i]->bw );
			}
			output += res*req[d].flow;
		}
		G->throughput = output;
		cout << "LB\t利用率 "<< obj <<"\t吞吐率 "<<output<<endl;
	}
	else{
		cout << "LB unfeasible"<<endl;
	}
	for(size_t i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

//数据流 吞吐率 P18
double throughput(CGraph *G,vector<demand> req,int ornum){    
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);

	int totalnum = req.size();
	IloArray<IloIntVarArray> x(env, totalnum); 
	for(int d = 0; d < totalnum; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); // num * G->m 的二维矩阵

	//优化目标
	IloExpr goal(env);
	IloExprArray Load(env, G->m);
	IloIntVarArray Res(env,totalnum,0,INF); //Res[d]表示业务d可用的剩余带宽
	// IloInfinity无解 ？？？   INF

	for(int i=0;i<G->m;i++)
		Load[i] = IloExpr(env);


	for(int i = 0; i < G->m; i++){
		for(int d = 0; d < totalnum; d++)
			Load[i] += x[d][i]*req[d].flow;
	}

	for(int d = 0; d < ornum; d++){
		for(int i = 0; i < G->m; i++)
			model.add(Res[d] <= ((1-x[d][i])*INF + (G->Link[i]->capacity - Load[i])));
		goal += Res[d]*req[d].flow; //why req[d].flow
		//goal+=Y[d]*reqL[d]->flow/(int)g->cost_best[reqL[d]->id];
	}

	model.add(IloMaximize(env,goal));

	//流量守恒约束
	for(int d = 0; d < totalnum; d++)
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		//带宽约束
		for(int i = 0; i < G->m; i++){
			IloExpr constraint(env);
			for(int d = 0; d <  totalnum; d++)
				constraint += req[d].flow*x[d][i];
			model.add(constraint <= G->Link[i]->capacity);  
		}

		ORsolver.setOut(env.getNullStream());
		double obj = INF;
		ORsolver.solve();
		if(ORsolver.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "throughput unfeasible" << endl;
		else{
			obj = ORsolver.getObjValue();
			// utilization
			double util = 0;
			for(int i = 0; i < G->m; i++){  
				double load = 0;
				for(int d = 0; d < totalnum; d++)
					load += ORsolver.getValue(x[d][i])*req[d].flow;
				util = max(util,load/G->Link[i]->capacity);
			}
			G->mlu = util;
			cout << "OR\t利用率 "<<util<< "\t吞吐率 "<<obj<<endl;
		}
		for(size_t i = 0; i < req.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}

//nashbw
double NashBW(CGraph *g,vector<demand> req){    
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);

	int num = req.size();
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, g->m, 0, 1); 

	//优化目标
	IloExpr goal(env);
	
	IloExprArray cost(env, num);
	for(int i = 0;i < num; i++)
		cost[i] = IloExpr(env);

	for(int d = 0; d < num; d++){
		for(int i = 0; i < g->m; i++)
			model.add(cost[d] <= ( (1-x[d][i])*INF + g->Link[i]->bw) ); //g->Link[i]->bw : residual bandwidth
		goal += cost[d];	
	}
	model.add(IloMaximize(env,goal));

	//流量守恒约束
	for(int d = 0; d < num; d++)
		for(int i = 0; i < g->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < g->adjL[i].size(); k++) // 出度边
				constraint += x[d][g->adjL[i][k]->id];
			for(unsigned int k = 0; k < g->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][g->adjRL[i][k]->id];
			// 出 - 入
			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		ORsolver.setOut(env.getNullStream());
		double obj = SMALL;
		ORsolver.solve();
		if(ORsolver.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "NashBW unfeasible" << endl;
		else{
			obj = ORsolver.getObjValue();
			//use
			for(int i = 0; i < g->m; i++){  
				double load = 0;
				for(int d = 0; d < num; d++){
					load += ORsolver.getValue(x[d][i])*req[d].flow;
				}
				g->Link[i]->use = load;
			}
		}
		for(size_t i = 0; i < req.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}

#endif