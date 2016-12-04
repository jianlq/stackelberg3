#include"evolutionbit.h"
#include"solver.h"
#include"nash.h"

// mlu和吞吐率
int main(){
	srand((unsigned)time(NULL));

	int Time = 2;
	int CASEnum= 20;	

	vector<double>CONSIDER;
	int CON_VALUE = 25;
	double c = 0;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}
	int LOOP  = 100;

	for(int i =0;i<Time;i++){

		int conN = CONSIDER.size();
		vector<double> smlu(conN,0) ;
		vector<double> sbw(conN,0);

		double mlu ;
		double mlubw ;

		double bwmlu;
		double bw;

		double nashmlu ;
		double nashbw;

		vector<int> successCase (conN, 0) ;

		int sucCasemlu = 0,sucCaseBW = 0,sucCaseNash = 0;

		for(int casenum = 0; casenum < CASEnum; casenum++){

			genGraph(9,54,"inputFile//graph.txt");
			genGraphOR(9,4,9,"inputFile//graphOR.txt");
			CGraph *G = new CGraph("inputFile//graph.txt");
			CGraph *GOR = new CGraph("inputFile//graphOR.txt");

			vector<demand> eqOR;
			vector<demand> eqTE;
			vector<demand> eqbase;

			//eqbase.clear();//background流 
			for(int i = 0; i < BGNUM; i++){
				int s = rand()%G->n, t;
				do{
					t = rand()%G->n;
				}while( s == t || G->canNotReach(s,t));
				eqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
			}

			////Overlay  产生demand
			//eqOR.clear(); 
			for(int i =0 ;i<GOR->m;i++)
				if(G->canNotReach(GOR->Link[i]->tail, GOR->Link[i]->head))
					continue;
				else
					eqOR.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,rand()%(MAXDEMAND)+1));

			//eqTE.clear(); 
			for(unsigned int i=0;i<eqOR.size();i++)
				eqTE.push_back(eqOR[i]);
			int ornum = eqTE.size();
			for(unsigned int i=0;i<eqbase.size();i++)
				eqTE.push_back(eqbase[i]);

			double mludic = 0, ordic = 0;

			G->clearOcc();
			mludic = LBdictor(G,eqTE,ornum);

			if( mludic + 1e-5 <= INF ){
				sucCasemlu++;
				mlu += mludic;
				mlubw += G->throughput;
			}
			else
				break;

			G->clearOcc();
			ordic = throughput(G,eqTE,ornum);

			if( ordic - 1e-5 >=SMALL ){
				sucCaseBW++;
				bw += ordic;
				bwmlu += G->mlu;
			}
			else
				break;

			G->clearOcc();
			if(!G->GAinit(eqTE)){
				cout << "*****GA init failed***** "<<endl;
				break;
			}

			//// nash	
			FILE *nash = fopen("outputFile//nash.csv", "a");
			int nacase = 0;
			double loopnashmlu=0,loopnashbw=0;
			fprintf(nash,"\n\n nash \n");
			for(int i =0;i<LOOP;i++){
				G->clearOcc();
				GOR->clearOcc();
				double curmlu = NashLB(G,GOR,eqTE);
				if(curmlu + 1e-5 >= INF){
					fprintf(nash,"NashEE unfeasible\n");
					break;
				}
				double curbw = NashBW(GOR,eqOR);
				if( curbw - 1e-5 <= SMALL){
					fprintf(nash,"NashOR unfeasible\n");
					break;
				}
				eqTE.clear();
				for(int i=0;i<GOR->m;i++){
					if(GOR->Link[i]->use>0)
						eqTE.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,GOR->Link[i]->use));
				}
				for(unsigned int i=0;i<eqbase.size();i++)
					eqTE.push_back(eqbase[i]);

				nacase++;
				loopnashmlu += curmlu;
				loopnashbw += curbw;
				fprintf(nash,"%f,%f\n",curmlu,curbw);
				cout<<curmlu<<"\t"<<curbw<<endl;
			}
			fclose(nash);

			if(nacase){
				nashmlu += (loopnashmlu/nacase);
				nashbw += (loopnashbw/nacase);
				sucCaseNash++;
			}

			FILE *cur = fopen("outputFile//mlubw.csv", "a");
			fprintf(cur,"\n\n\n,%d\n",casenum);
			fprintf(cur,",LB,,%f,%f\n",mludic,G->throughput);
			fprintf(cur,",OR,,%f,%f\n",G->mlu,ordic);
			fprintf(cur,",Nash,,%f,%f\n\n",loopnashmlu/nacase,loopnashbw/nacase);
			fclose(cur);

			for(unsigned int con = 0;con < CONSIDER.size();con++){

				int n = 150;//种群个体数目
				int m = eqTE.size();
				evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,mludic,ordic,CONSIDER[con]);
				(popubit).evolution();
				cout<<"S\t"<<popubit.hero.mlu<<"\t"<<popubit.hero.throughput <<endl;

				if( (popubit.hero.mlu +1e-5) >= INF ||  (popubit.hero.throughput  - 1e-5) < SMALL ){
					break;
				}

				cur = fopen("outputFile//mlubw.csv", "a");	
				fprintf(cur,",S,%f,%f,%f\n",CONSIDER[con],popubit.hero.mlu,popubit.hero.throughput);
				fclose(cur);

				successCase[con] += 1;
				smlu[con] += popubit.hero.mlu;
				sbw[con] += popubit.hero.throughput;

			} // end of CONSIDER for

			delete G;
			delete GOR;

		} // end of CASENum for

		FILE *res = fopen("outputFile//result.csv", "a");		
		fprintf(res,"\n\n case average ,%d\n",CASEnum);
		fprintf(res,",,CONSIDER,successCase,Energy Efficiency,throughput\n");
		for(unsigned int con = 0;con < CONSIDER.size();con++){
			fprintf(res, ",S,%f,%d,%f,%f\n",CONSIDER[con],successCase[con],smlu[con]/successCase[con],sbw[con]/successCase[con]); 
		}
		fprintf(res,",,,successCase,Load Balance,throughput\n");
		fprintf(res, "\n\n,LB,,%d,%f,%f\n",sucCasemlu,mlu/sucCasemlu,mlubw/sucCasemlu);
		fprintf(res, ",OR,,%d,%f,%f\n",sucCaseBW,bw/sucCaseBW,bwmlu/sucCaseBW);
		fprintf(res, ",Nash,,%d,%f,%f\n",sucCaseNash,nashmlu/sucCaseNash,nashbw/sucCaseNash);
		fclose(res);

	} // end of Time for
	system("pause");
	return 0;	
}