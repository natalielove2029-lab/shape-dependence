/***********************************************************/
/* Rendezetlen kvantum spinrendszerek renormalizacioja     */
/*                                                         */
/* Alapeset: RTFIC lanc, vagy letra                        */
/*                                                         */
/* Kovacs Istvan                                           */
/***********************************************************/
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <set>
#include <iostream>
#include <time.h>
#include <vector>
#include "binheap.h" 
#include "MersenneTwister.h"
#include <map>
#include <iomanip>
/****************************************************************
*Itt kezdodik a program						 *
****************************************************************/

using namespace std;
//TRandomMersenne randgen(0);
//**************************************
typedef set<int, std::greater<int> > halmaz;
typedef set<int, std::less<int> > felhalmaz;
typedef std::map<int, int> intmap;
typedef std::map<int, double> doublemap;

struct edge_type 
	{int                        start, end, h_id;
	double                     weight;
	};

typedef std::map<int, edge_type>::iterator     edge_iterator;
// Egy pont tulajdonsagai
struct node_type 
	{int                        	id, mu,gap, h_id;
	double                     	weight;
	std::map<int, edge_iterator>   edges;
	felhalmaz   his;
	}; 
	
typedef BinaryHeap<double, int, std::less_equal<double> >    min_heap; 
min_heap                   m_h;
min_heap::pair_type        m_he,m_hh;

int minsize=0, hsize, midway,realizacio, voltej,perk=0, duplae, Lmin, 
	Lmax, Edges, Nodes, graphform, rep, muin, seed, d, ido, s, kiiratas, NODE_START_ID=0, EDGE_START_ID,
	pontokszama, elekszama;

double minweight=0,lastj,INF=-1000000.0,EDGE_FACTOR=1,hmin, hmax, Jmin, Jmax,  omega, magne;

typedef std::map<int, edge_iterator>::iterator iteredge;
halmaz jmaxok, vent;
std::map<int, node_type>                                 nodes;
std::map<int, edge_type>                                 edges;
	

int ipow(int alap, int kitevo)
    {//pozitiv egesz kitevok kiszamitasara szolgal
    //Ha a kitevo nulla, akkor 1-et ad vissza, negativ kitevonel pedig 0-t!
    if (kitevo>=0)
        {int eredmeny=1;
        while (kitevo>=1)
            {if ((kitevo%2)==0)
                {alap*=alap;
                kitevo/=2;
                }
            else
                {eredmeny*=alap;
                --kitevo;
                }
            }
        return (eredmeny);
        }
    else
        {return 0;
        }
    };

void entropiakiiras(int xi, ofstream *fne)
	{int ai;
	halmaz::iterator hsi;
	if (nodes[xi].his.size()>1)	
		{//ha meg nem irtuk ki:(az elso eleme jol jellemzi a klasztert!)
		if (vent.find(*(nodes[xi].his.begin()))==vent.end())
			{//el kell tenni, ha mar duplaztunk:
			if (duplae==1)
				{vent.insert(*(nodes[xi].his.begin()));
				}
			//ki kell irni a klasztert:
			hsi=nodes[xi].his.begin();	
			ai=*hsi;
			if (nodes[xi].weight>INF)
				{*fne << ai;
				hsi++;
				for(;hsi!=nodes[xi].his.end();hsi++)
					{*fne << "\11"<<(*hsi)-ai;
					ai=*hsi;
					}
				*fne << "\n";
				}
			else
				{*fne << ai;
				hsi++;
				for(;hsi!=nodes[xi].his.end();hsi++)
					{*fne << "\11"<<(*hsi)-ai;
					ai=*hsi;
					}
				*fne << "\n";
				}
			}
		else
			{vent.erase(*nodes[xi].his.begin());
			}
		}
	nodes[xi].his.clear();
	//cout << "_e_";
	}
	

	
void lastminute(ofstream *fneu)
	{//pozitiv egesz kitevok kiszamitasara szolgal
	//Ha a kitevo nulla, akkor 1-et ad vissza, negativ kitevonel pedig 0-t!
	//cout << "_l";
	int x;
	halmaz::iterator hs;
	iteredge masx;
	int ehelyett=0;
	while (!m_h.isEmpty())
		{m_h.pop();
		}	
	for (std::map<int, node_type>::iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
		{x=iter->first;
	if (kiiratas>=1) {entropiakiiras(x,fneu);}
		if (nodes[x].weight<=INF)
			{ehelyett+=nodes[x].mu/muin;
			}
		else
			{if (nodes[x].weight<minweight)
				{minsize=nodes[x].mu/muin;
				//if (nodes[x].weight>INF)
				minweight=nodes[x].weight;
				}
			}
		//pontokszama--;
		/*for (masx=nodes[x].edges.begin(); masx!=nodes[x].edges.end(); masx++)
			{edges.erase(masx->second->first);
			nodes[masx->first].edges.erase(x);
			//cout << "ED_"<< masx->second->id << "\11";
			elekszama--;
			}
		nodes[x].edges.clear();	*/	
		}
	nodes.clear();
	edges.clear();
	if (ehelyett>0)
		{minsize=ehelyett;
		}
	//hsize=0;	
	//cout << "_L";
	};	
	
int lef(int alap)
	{//duplazott rendszerben az alteregot adja meg 
	int eredmeny=alap;
	if (duplae==1)
		{if (2*alap>=Nodes)	
			{eredmeny=alap-div(Nodes,2).quot;
			}
		}
	return (eredmeny);
	};	

void in(int p, int q)
	{//p:honnan, q: hova
	//edges[Edges].id=Edges;
	edges[Edges].start=p;
	edges[Edges].end=p+q;
	nodes[p].edges[p+q]=edges.find(Edges);
	nodes[p+q].edges[p]=edges.find(Edges);
	Edges++;
	}	

void hinsert(int itt, ofstream *fneuq)
	{if (nodes[itt].weight<minweight)
		{minsize=nodes[itt].mu/muin;
		if (nodes[itt].weight>INF)
			{minweight=nodes[itt].weight;
			}
		}
	nodes[itt].gap=-1;
	if (kiiratas>=1) {entropiakiiras(itt,fneuq);}
	//a heapbol ki kell szedni az eleit!!, illetve a midwayhez menok hosszat csak frissiteni kell!
	for (iteredge eli=nodes[itt].edges.begin(); eli!=nodes[itt].edges.end() ; eli++)
		{if (nodes[eli->first].gap>-1)	
			{//modositani kell a heapben a hosszat
			m_h.move(eli->second->second.h_id,-2*eli->second->second.weight+nodes[itt].weight);
			//cout << "mo:" <<eli->second->first<<"\11";	
			//ha ide mutatott a gapje, akkor ujra kellene szamolni, ha megvaltozna az az el!
			/*if ((nodes[eli->first].gap-1)==eli->second->first)
				{gapupdate(eli->first);
				}*/
			}
		else
			{//ha lokmaxhoz ment az el, ki kell szedni a heapbol:
			//cout << "ki:" <<eli->second->first<<"\11";	
			m_h.pop(eli->second->second.h_id);
			}
		}
	};			

	
void egypontelmaxe(int ip)
	{// Valtozas: a midway pont osszes olyan midwayekhez meno elet el kell tenni a heapbe, ami nem lokmax, es a pont sulyat is mindig el kell tenni!
	if (nodes[ip].gap>0)
		{int egyik, volte=-1;
		//a pont sulyat a pontbesoroloperem-ben mar eltettuk!
		//nodes[ip].h_id=m_h.push(-nodes[ip].weight,-ip-1);	
		if (edges[nodes[ip].gap-1].start==ip)
			{egyik=edges[nodes[ip].gap-1].end;
			}
		else
			{egyik=edges[nodes[ip].gap-1].start;
			}
		if (nodes[egyik].gap==nodes[ip].gap)
			{//csak az egyik vegenel kell eltenni a lokmax elt:
			if (egyik>ip)
				{jmaxok.insert(nodes[ip].gap-1);
				volte=egyik;
				}
			}
		
		for (iteredge eli=nodes[ip].edges.begin(); eli!=nodes[ip].edges.end() ; eli++)
			{if (nodes[eli->first].gap==-1)	
				{//az elnel nem is kell eltarolni, hogy ut-e, a vegei lokmaxsaga elarulja
				eli->second->second.h_id=m_h.push(-2*eli->second->second.weight+nodes[eli->first].weight,eli->second->first);
				//cout << "Be:" <<eli->second->first<<"ip"<<ip<< "ef"<< eli->first<<"\11";	
				}	
			else
				{//csak akkor kell eltenni, ha nem ez az el a lokmax:
				//egy midwayhez meno elt csak egyszer szabad megnezni:
				if (eli->first>ip)
					{if (eli->first!=volte)
						{eli->second->second.h_id=m_h.push(-eli->second->second.weight,eli->second->first);
						//cout << "bm:" <<eli->second->first<<"ip"<<ip<< "ef"<< eli->first<< "\11";	
						}
					}
				}
			}
		}
	};

void pontbesoroloperem(double jlim)
	{//A legelso besorolast futtatja le a grafon.
	hsize=0;
	int vane=0, ize=0;
	double kappa;
	for(int i=0; i<Nodes; i++)
		{//Ha nem perempont:	
		ize=0;
		if (nodes[i].gap==-3) {ize=1;};
			{if (nodes[i].weight<=jlim)	
				{vane=-1;
				kappa=nodes[i].weight;
				for (iteredge eli=nodes[i].edges.begin(); eli!=nodes[i].edges.end() ; eli++)
					{if (eli->second->second.weight>=kappa)
						{vane=eli->second->first;
						kappa=eli->second->second.weight;		
						//Ha van perem szomszedja, akkor a vane-jet -2-re allitjuk:
						//de csak ha meg minusz 1
						}
					}
				if (vane>-1)
					{nodes[i].gap=vane+1;
					if (ize==0) {midway++;}
					//el kell tenni a heapbe a kulso teret:
					//cout << i<< "\11";
					nodes[i].h_id=m_h.push(-nodes[i].weight,-i-1);
					}
				else
					{//a pont maximalis
					nodes[i].gap=-1;
					if (ize==1) {midway--;}
					}
				}
			else
				{//a pont maximalis
				nodes[i].gap=-1;
				if (ize==1) {midway--;}
				}
			}
		}
	//elek besorolasa
	//cout <<"_e";
	for(int i=0; i<Nodes; i++)
		{//cout <<"_"<<i;
		egypontelmaxe(i);
		}			
	//cout <<"_f";	
	};


void grafparameterek(int bemeret, int sz)
	{//A grafok pontjainak illetve eleinek szamat allitja be
	switch (graphform)
		{case 0: 
			{//linearis lanc
			Edges=bemeret;
			Nodes=bemeret;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;
			break;
			}
		case 1: 
			{//L es d: ket merete, d-vel vonaltol sikig valtozik, d vegig fix, csak L novekszik a merettel
			Nodes=bemeret*sz;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;//(Igazabol =2 lenne a konzisztens megoldas!)
			break;
			}
		case 2: 
			{//Negyzetracstol kockaig
			//Egyelore csak sima negyzetracs
			Edges=bemeret*bemeret*2;
			Nodes=bemeret*bemeret;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;	
			break;
			}
		case 3: 
			{//Kocka:
			Nodes=bemeret*bemeret*bemeret;
			Edges=Nodes*3;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;	
			break;
			}	
		case 4: 
			{//4d Kocka:
			Nodes=bemeret*bemeret*bemeret*bemeret;
			Edges=Nodes*4;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;	
			break;
			}				
		case 9: 
			{//Teljes graf
			Nodes=bemeret;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=(double (Nodes-1))/2;
			//cout << EDGE_FACTOR << '\11';
			break;
			}
		case 10: 
			{//Vegtelen dimenzios Bethe-racs
			//Nodes=1;
			/*for (int i=0; i<=bemeret; i++)
			{Nodes*=sz;
			}*/
			Nodes=ipow(sz, bemeret);
			Nodes--;
			Nodes/=sz-1;
			EDGE_START_ID=Nodes;
			EDGE_FACTOR=1;
			break;
			}
		}
	};



void pajek(int pillido)
	{char paj[255],clu[255];	
	//int pillido=(int) time(NULL);
	sprintf(paj,"RG_%d.net",pillido);
	ofstream lpaj(paj, ios::out);
	if (!lpaj)
		{cout << "Kimeneti pajek fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
		sprintf(clu,"RG_%d.clu",ido);
	ofstream lclu(clu, ios::out);
	if (!lclu)
		{cout << "Kimeneti pajek clu fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	lpaj << "*Vertices" << "  " << Nodes << '\n';	
	lclu << "*Vertices" << "  " << Nodes << '\n';	
	int i=1;
	std::map<int, int> neve;
	for (std::map<int, node_type>::iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
		{//lpaj << i << "  " << "\""<<iter->first<<"\"" << '\n'; 
		lpaj << i << "  " << "\""<<iter->second.weight<<"\"" << '\n';
		neve[iter->first]=i;
		//lclu << iter->second.hova << '\n';
		i++;
		}
	lpaj << "*Edges" << '\n';		
	for (std::map<int, edge_type>::iterator iter=edges.begin(); iter!=edges.end(); iter++)
		{lpaj << neve[iter->second.start] << "  " << neve[iter->second.end] << "  " << -iter->second.weight << '\n';
		}		
lpaj.close();		
lclu.close();
neve.clear();
};

	
void elekbehuzasapbc(int bemeret, int sz)
	{//elek behuzasa
	if(graphform==0)
		{edges.clear();
		nodes.clear();
		//linearis lanc
		int i;
		for(i=0; i<Edges-1; i++)
			{//edges[i].id=i;
			edges[i].start=i;
			edges[i].end=i+1;
			}
		//pbc
		if (Nodes>2)
		{
		//edges[Edges-1].id=Edges-1;
		edges[Edges-1].start=Edges-1;
		edges[Edges-1].end = 0;
            nodes[0].gap=-3;
            nodes[Edges-1].gap=-3;
            midway+=2;
		}
		//Torolni kell a korabbi ellistat es utana kell beilleszteni az ujat
		for(i=1; i<Nodes-1; i++)
			{nodes[i].id=i;
			nodes[i].edges.clear();
			nodes[i].edges[i+1]=edges.find(i);
			nodes[i].edges[i-1]=edges.find(i-1);
			}
		nodes[0].id=0;	
		nodes[0].edges.clear();
		if (Nodes>2)
		{
		nodes[0].edges[Edges-1]=edges.find(Edges-1);
		}
		nodes[0].edges[1]=edges.find(0);
		nodes[Nodes-1].id=Nodes-1;	
		nodes[Nodes-1].edges.clear();
		if (Nodes>2)
		{
		nodes[Nodes-1].edges[0]=edges.find(Edges-1);
		}
		nodes[Nodes-1].edges[Nodes-2]=edges.find(Edges-2);
		if (Nodes==2)
		{Edges--;
		}
		}
	else
		{if(graphform==1)
			{//L es d: ket merete, d-vel vonaltol sikig valtozik
			
			}	
		else
			{if(graphform==2)
				{//Negyzetracs	
				edges.clear();
			nodes.clear();
			Edges=0;
			int i;
			for(i=0; i<Nodes; i++)
				{nodes[i].id=i;
				nodes[i].edges.clear();
				}
			//x irany:
			//akkor van a jobb szelen, ha div(x,bemeret).rem=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,bemeret).rem==(bemeret-1))
					{//atkotjuk a bal szelevel
					in(i,-(bemeret-1));
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)].gap=-3;
					//midway++;
					}
				else
					{//jobbra kotjuk
					in(i,1);
					}
				}	
			midway+=2*bemeret;
			//y irany:
			//akkor van az also szelen, ha div(i,bemeret).quot=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,bemeret).quot==(bemeret-1))
					{//atkotjuk a tetejevel
					in(i,-(bemeret-1)*bemeret);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*bemeret].gap=-3;
					//midway++;
					}
				else
					{//le kotjuk
					in(i,bemeret);
					}
				}
			midway+=2*(bemeret-2);				
				}
		else
		{if(graphform==3)
			{//Kocka	
			edges.clear();
			nodes.clear();
			Edges=0;
			int i;
			//cout << "Nodes:" << Nodes << "\n";
			for(i=0; i<Nodes; i++)
				{nodes[i].id=i;
				//nodes[i].m=-1;
				//if (i==274) {cout << "!";}
				nodes[i].gap=0;
				//nodes[i].q=-1;
				nodes[i].edges.clear();
				}
			//x irany:
			//akkor van a jobb szelen, ha div(x,bemeret).rem=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,bemeret).rem==(bemeret-1))
					{//atkotjuk a bal szelevel
					in(i,-(bemeret-1));
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)].gap=-3;
					//midway++;
					}
				else
					{//jobbra kotjuk
					in(i,1);
					}
				}
			midway+=2*bemeret*bemeret;	
			//y irany:
			//akkor van az also szelen, ha div(i,bemeret).quot=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(div(i,bemeret*bemeret).rem,bemeret).quot==(bemeret-1))
					{//atkotjuk a tetejevel
					//	cout << "!" << i;
					in(i,-(bemeret-1)*bemeret);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*bemeret].gap=-3;
					//midway++;
					}
				else
					{//le kotjuk
					in(i,bemeret);
					//cout << "?" << i;
					}
				}	
				midway+=2*bemeret*(bemeret-2);
			//z irany:
			//akkor van az felso szelen, ha div(i,bemeret*bemeret).quot=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,bemeret*bemeret).quot==(bemeret-1))
					{//atkotjuk az aljaval
					in(i,-(bemeret-1)*bemeret*bemeret);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*bemeret*bemeret].gap=-3;
					//midway++;
					}
				else
					{//fel kotjuk
					in(i,bemeret*bemeret);
					}
				}
			midway+=2*(bemeret-2)*(bemeret-2);				
			//formalisan analog modon lehet a tovabbi dimenziokat bevezetni.
			}
		else
		{if(graphform==4)
			{//4d kocka
			edges.clear();
			nodes.clear();
			Edges=0;
			int i, k;
			//cout << "Nodes:" << Nodes << "\n";
			for(i=0; i<Nodes; i++)
				{nodes[i].id=i;
				//nodes[i].m=-1;
				//if (i==274) {cout << "!";}
				nodes[i].gap=0;
				//nodes[i].q=-1;
				nodes[i].edges.clear();
				}
			//x irany:
			//akkor van a jobb szelen, ha div(x,bemeret).rem=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,bemeret).rem==(bemeret-1))
					{//atkotjuk a bal szelevel
					in(i,-(bemeret-1));
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)].gap=-3;
					//midway++;
					}
				else
					{//jobbra kotjuk
					in(i,1);
					}
				}
			midway+=2*bemeret*bemeret*bemeret;	
			//y irany:
			//akkor van az also szelen, ha div(i,bemeret).quot=bemeret-1;
			k=bemeret*bemeret;
			for(i=0; i<Nodes; i++)
				{//if (div(div(i,k).rem,bemeret).quot==(bemeret-1))
				if (div(div(div(i,k*bemeret).rem,k).rem,bemeret).quot==(bemeret-1))
					{//atkotjuk a tetejevel
					//	cout << "!" << i;
					in(i,-(bemeret-1)*bemeret);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*bemeret].gap=-3;
					//midway++;
					}
				else
					{//le kotjuk
					in(i,bemeret);
					//cout << "?" << i;
					}
				}	
				midway+=2*k*(bemeret-2);
			//z irany:
				
			//akkor van az felso szelen, ha div(i,bemeret*bemeret).quot=bemeret-1;
			for(i=0; i<Nodes; i++)
				{//if (div(i,k).quot==(bemeret-1))
				if (div(div(i,k*bemeret).rem,k).quot==(bemeret-1))
					{//atkotjuk az aljaval
					in(i,-(bemeret-1)*k);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*k].gap=-3;
					//midway++;
					}
				else
					{//fel kotjuk
					in(i,k);
					}
				}
			midway+=2*(bemeret-2)*(bemeret-2)*bemeret;			
			//u irany:
			k*=bemeret;
			//akkor van az felso szelen, ha div(i,bemeret*bemeret*bemeret).quot=bemeret-1;
			for(i=0; i<Nodes; i++)
				{if (div(i,k).quot==(bemeret-1))
					{//atkotjuk az aljaval
					in(i,-(bemeret-1)*k);
					nodes[i].gap=-3;
					//midway++;
					nodes[i-(bemeret-1)*k].gap=-3;
					//midway++;
					}
				else
					{//fel kotjuk
					in(i,k);
					}
				}
			midway+=2*(bemeret-2)*(bemeret-2)*(bemeret-2)*bemeret;	
			}		
		}		
		//
		}
		}
		}
	};
	
void elekbehuzasateljes(int bemeret)
	{//elek behuzasa
	if(graphform==4)
		{if (Nodes>1)
			{edges.clear();
			nodes.clear();
			Edges=0;
			for(int i=0; i<Nodes; i++)
				{nodes[i].id=i;
				nodes[i].edges.clear();
				}
			for(int i=0; i<Nodes; i++)
				{for(int j=i+1; j<Nodes; j++)
					{nodes[i].edges[j]=edges.find(Edges);
					nodes[j].edges[i]=nodes[i].edges[j];
					//edges[Edges].id=Edges;
					edges[Edges].start=i;
					edges[Edges].end=j;
					Edges++;
					//cout << Edges << '\n';
					}
				}
			}
		}
	};

void elekbehuzasaBethe(int bemeret, int sz)
	{//elek behuzasa
	if(graphform==5)
		{//int kezdo=0;
		//int kovkezdo=0;
		if (Nodes>1)
			{edges.clear();
			nodes.clear();
			Edges=0;
			for(int i=0; i<Nodes; i++)
				{nodes[i].id=i;
				nodes[i].edges.clear();
				}
			int	kezdo=0;
			int	kovkezdo=0;
			for(int i=0; i<bemeret; i++)
				{kezdo+=ipow(sz, i-1);
				kovkezdo+=ipow(sz, i);
				for(int j=0; j<=sz; j++)
					{for(int k=0; k<=sz; j++)
						{nodes[kezdo+i].edges[kovkezdo+j*sz+k]=edges.find(Edges);
						nodes[kovkezdo+j*sz+k].edges[kezdo+i]=nodes[kezdo+i].edges[kovkezdo+j*sz+k];
						//edges[Edges].id=Edges;
						edges[Edges].start=kezdo+i;
						edges[Edges].end=kovkezdo+j*sz+k;
						Edges++;
						//cout << Edges << '\n';
						}
					}
				}
			}
		}
	};
	
void elekbehuzasafree(int bemeret, int sz)
	{//elek behuzasa
	if(graphform==0)
		{Edges--;
		edges.clear();
		nodes.clear();
		//linearis lanc
		for(int i=0; i<Edges; i++)
			{//edges[i].id=i;
			edges[i].start=i;
			edges[i].end=i+1;
			}
		//edges[Edges-1].end = 0;
		//Torolni kell a korabbi ellistat es utana kell beilleszteni az ujat
		for(int i=1; i<Nodes-1; i++)
			{nodes[i].id=i;
			nodes[i].edges.clear();
			nodes[i].edges[i+1]=edges.find(i);
			nodes[i].edges[i-1]=edges.find(i-1);
			}
		nodes[0].id=0;	
		nodes[0].edges.clear();
		//nodes[0].edges[Edges-1]=&edges[Edges-1];
		nodes[0].edges[1]=edges.find(0);
		nodes[Nodes-1].id=Nodes-1;	
		nodes[Nodes-1].edges.clear();
		//nodes[Nodes-1].edges[0]=&edges[Edges-1];
		nodes[Nodes-1].edges[Nodes-2]=edges.find(Edges-1);
		}
	else
		{if(graphform==1)
			{//L es d: ket merete, d-vel vonaltol sikig valtozik
			edges.clear();
			nodes.clear();
			Edges=0;
			for(int i=0; i<Nodes; i++)
				{nodes[i].id=i;
				nodes[i].edges.clear();
				}
			for(int i=0; i<sz-1; i++)
				{for(int j=0; j<bemeret-1; j++)
					{nodes[i*bemeret+j].edges[i*bemeret+j+1]=edges.find(Edges);
					nodes[i*bemeret+j+1].edges[i*bemeret+j]=nodes[i*bemeret+j].edges[i*bemeret+j+1];
					//edges[Edges].id=Edges;
					edges[Edges].start=i*bemeret+j;
					edges[Edges].end=i*bemeret+j+1;
					Edges++;
					nodes[i*bemeret+j].edges[i*bemeret+j+bemeret]=edges.find(Edges);
					nodes[i*bemeret+j+bemeret].edges[i*bemeret+j]=nodes[i*bemeret+j].edges[i*bemeret+j+bemeret];
					//edges[Edges].id=Edges;
					edges[Edges].start=i*bemeret+j;
					edges[Edges].end=i*bemeret+j+bemeret;
					Edges++;
					}
				}
			//free
			for(int i=0; i<sz-1; i++)
				{//utolso oszlop
				/*if (bemeret>2)
					{
					nodes[(i+1)*bemeret-1].edges[i*bemeret]=&edges[Edges];
					nodes[i*bemeret].edges[(i+1)*bemeret-1]=&edges[Edges];
					edges[Edges].id=Edges;
					edges[Edges].start=(i+1)*bemeret-1;
					edges[Edges].end=i*bemeret;
					Edges++;
					}*/
				nodes[(i+1)*bemeret-1].edges[(i+2)*bemeret-1]=edges.find(Edges);
				nodes[(i+2)*bemeret-1].edges[(i+1)*bemeret-1]=nodes[(i+1)*bemeret-1].edges[(i+2)*bemeret-1];
				//edges[Edges].id=Edges;
				edges[Edges].start=(i+1)*bemeret-1;
				edges[Edges].end=(i+2)*bemeret-1;
				Edges++;
			}
			for(int i=0; i<bemeret-1; i++)
				{//utolso sor
				nodes[(sz-1)*bemeret+i].edges[(sz-1)*bemeret+i+1]=edges.find(Edges);
				nodes[(sz-1)*bemeret+i+1].edges[(sz-1)*bemeret+i]=nodes[(sz-1)*bemeret+i].edges[(sz-1)*bemeret+i+1];
				//edges[Edges].id=Edges;
				edges[Edges].start=(sz-1)*bemeret+i+1;
				edges[Edges].end=(sz-1)*bemeret+i;
				Edges++;
				/*if (sz>2)
					{
					nodes[(sz-1)*bemeret+i].edges[i]=&edges[Edges];
					nodes[i].edges[(sz-1)*bemeret+i]=&edges[Edges];
					edges[Edges].id=Edges;
					edges[Edges].start=i;
					edges[Edges].end=(sz-1)*bemeret+i;
					Edges++;
					}*/
				}
			
			}	
		else
			{if(graphform==2)
				{//Negyzetracs	
				}
			}
		}
	};	

void elekbehuzasa(int bemerete, int sze)
	{if (graphform<=1)
		{if (s==0)
			{elekbehuzasapbc(bemerete,sze);
			}
		else
			{elekbehuzasafree(bemerete, sze);
			}
		}
	else
		{switch (graphform)
		{case 2:
			{elekbehuzasapbc(bemerete,sze);
			break;
			}
		case 3: 
				{//cout << "A" << '\n';
				elekbehuzasapbc(bemerete,sze);
				break;
				}
		case 4: 
				{//cout << "A" << '\n';
				elekbehuzasapbc(bemerete,sze);
				break;
				}				
			case 9: 
				{//cout << "A" << '\n';
				elekbehuzasateljes(bemerete);
				break;
				}
			case 10: 
				{elekbehuzasaBethe(bemerete, sze);
				break;
				}
			}
		}
	};
	
void duplazas()
	{//vegigmegyunk a pontokon es Nodes-zal odebb megkeltjuk oket es az eleiket
	std::map<int, edge_type>::iterator mask=edges.end();
	mask--;
	int utolso=mask->first;
	for (std::map<int, edge_type>::iterator masx=edges.begin(); masx->first<=utolso; masx++)
		{//edges[Edges+masx->first].id=Edges+masx->first;
		edges[Edges+masx->first].weight=masx->second.weight;
		edges[Edges+masx->first].start=Nodes+masx->second.start;
		edges[Edges+masx->first].end=Nodes+masx->second.end;
		}
	for (std::map<int, node_type>::iterator iter=nodes.begin(); iter->first<Nodes; iter++)
		{nodes[iter->first+Nodes].id=iter->first+Nodes;
		nodes[iter->first+Nodes].mu=iter->second.mu;
		nodes[iter->first+Nodes].weight=iter->second.weight;
		//A .his-et is at kell masolni:
		if (iter->second.his.size()>0)
			{nodes[iter->first+Nodes].his.insert(iter->second.his.begin(),iter->second.his.end());
			}
		if (iter->second.gap<0)
			{if (iter->second.gap==-3)
				{midway++;
				}
			nodes[iter->first+Nodes].gap=iter->second.gap;
			}
		else
			{//a maximalis elet attettuk Edges-szel odebb:
			nodes[iter->first+Nodes].gap=iter->second.gap+Edges;
			//midway++;
			}
		for (iteredge masi=iter->second.edges.begin(); masi!=iter->second.edges.end(); masi++)
			{nodes[iter->first+Nodes].edges[masi->first+Nodes]=edges.find(masi->second->first+Edges);
			}
		}
	};	
	
	

void ch(int p, int q)
	{//Duplazas utan a pbc peremek atkotese
	//hosszu el atkotese	
	edges[nodes[p].edges[p+q]->first].start=Nodes+p+q;
	edges[nodes[p].edges[p+q]->first].end=p;
	//ugyanezt frissiteni kell a pontoknal is:
	nodes[p].edges[Nodes+p+q]=edges.find(nodes[p].edges[p+q]->first);
	//nodes[Nodes+p+q].edges[p]=edges.find(nodes[p].edges[p+q]->first);
	nodes[Nodes+p+q].edges[p]=edges.find(nodes[p].edges[p+q]->first);
	//rovid el atkotese
	edges[nodes[Nodes+p].edges[Nodes+p+q]->first].start=p+q;
	edges[nodes[Nodes+p].edges[Nodes+p+q]->first].end=Nodes+p;
	//ugyanezt frissiteni kell a pontoknal is:
	nodes[p+q].edges[Nodes+p]=edges.find(nodes[Nodes+p].edges[Nodes+p+q]->first);
	nodes[Nodes+p].edges[p+q]=edges.find(nodes[Nodes+p].edges[Nodes+p+q]->first);
	//regi elek torlese
	nodes[p].edges.erase(p+q);
	nodes[p+q].edges.erase(p);
	nodes[Nodes+p].edges.erase(Nodes+p+q);
	nodes[Nodes+p+q].edges.erase(Nodes+p);
	//cout << "." <<"\11";
	};
	

void cross1dperem(int bemerete)
	{//Minden hatarnal atkotjuk, es mindent decimalhatova teszunk
	//x irany
	ch(0,bemerete-1);
	Nodes*=2;
	pontokszama*=2;
	elekszama*=2;
	Edges*=2;
	};	
	
void cross2dperem(int bemerete, int sze)
	{//Minden hatarnal atkotjuk, es mindent decimalhatova teszunk
	int j;
	//x irany
	for(j=0; j<sze; j++)
		{ch(j*bemerete,bemerete-1);
		}
	//y irany	
	for(j=0; j<bemerete; j++)
		{ch(j,(sze-1)*bemerete);
		}	
	Nodes*=2;
	pontokszama*=2;
	elekszama*=2;
	Edges*=2;
	};	

void cross3dperem(int bemerete)
	{//Minden hatarnal atkotjuk, es mindent decimalhatova teszunk
	int j,k;
	//cout << "Q" << midway << "\11";
	for(k=0; k<bemerete; k++)
		{//x irany
		for(j=0; j<bemerete; j++)
			{ch(k*bemerete*bemerete+j*bemerete,bemerete-1);
			}	
		//y irany	
		for(j=0; j<bemerete; j++)
			{ch(k*bemerete*bemerete+j,(bemerete-1)*bemerete);
			}
		}		
	//cout << "W" << midway << "\11";
	//z irany:
	for(k=0; k<bemerete*bemerete; k++)
		{ch(k,(bemerete-1)*bemerete*bemerete);
		}
	Nodes*=2;
	pontokszama*=2;
	elekszama*=2;
	Edges*=2;
	};	

void cross4dperem(int bemerete)
	{//Minden hatarnal atkotjuk, es mindent decimalhatova teszunk
	int j,k,m,l,p;
	l=bemerete*bemerete;
	p=l*bemerete;
	//cout << "Q" << midway << "\11";
	for(m=0; m<bemerete; m++)
		{
	for(k=0; k<bemerete; k++)
		{//x irany
		for(j=0; j<bemerete; j++)
			{ch(k*l+j*bemerete+m*p,bemerete-1);
			}	
		//y irany	
		for(j=0; j<bemerete; j++)
			{ch(k*l+j+m*p,(bemerete-1)*bemerete);
			}
		}		
	//cout << "W" << midway << "\11";
	//z irany:
	for(k=0; k<l; k++)
		{ch(k+m*p,(bemerete-1)*l);
		}
	}
	//u irany:
	l*=bemerete;	
	for(k=0; k<l; k++)
		{ch(k,(bemerete-1)*l);
		}	
		
	Nodes*=2;
	pontokszama*=2;
	elekszama*=2;
	Edges*=2;
	};			
	
int Jkilld(int index, ofstream *fne, double neww)
	{//lokalis maximum J decimalasa es a kornyezet rendbetetele
	int hogyan=0;
	int kovetkezo=-1;
	int tempo=-1;
	double kappa=0;
	omega=neww;
	midway--;
	int remained, decimated;
	if (omega<lastj) 
		{if (omega>INF) 
			{lastj=omega;
			}
		}
	remained=edges[index].start;
	decimated=edges[index].end;
	int am,bm;
	am=nodes[decimated].his.size();
	bm=nodes[remained].his.size();
	if (am==0) {nodes[decimated].his.insert(lef(decimated)); am++;}
	if (bm==0) {nodes[remained].his.insert(lef(remained)); bm++;}
	if (am>bm)
		{//a kisebbet decimaljuk
		decimated=edges[index].start;
		remained=edges[index].end;
		}
	m_h.pop(nodes[decimated].h_id);
	//Tortenetek egyesitese:
	nodes[remained].his.insert(nodes[decimated].his.begin(),nodes[decimated].his.end());
	nodes[decimated].his.clear();
	if (nodes[remained].his.size()<(am+bm))
		{hogyan=1;
		//cout << "T";
		perk=1;
		}
	iteredge iter;
	//if (kiiratas>=3) {tesztkiiras(remained, decimated, lki);}
	pontokszama--;
	edges.erase(nodes[decimated].edges[remained]->first);
	nodes[decimated].edges.erase(remained);
	nodes[remained].edges.erase(decimated);
	elekszama--;
	nodes[remained].mu+=nodes[decimated].mu;
	//a fenti bonyolult helyett:	
	if (hogyan==1)
	{nodes[remained].weight=INF;
	}
	else
	{if ((nodes[decimated].weight==INF) || (nodes[remained].weight==INF))
		{nodes[remained].weight=INF;
		}
	else
		{nodes[remained].weight+=nodes[decimated].weight-omega;
		}
	}
	///////
	for (iter=nodes[decimated].edges.begin(); iter!=nodes[decimated].edges.end(); iter++)
		{//iter->first	az a masik vegpont!
		if (nodes[remained].edges.find(iter->first) !=nodes[remained].edges.end())
			{//r megfelelo elenek sulya modosul	altalaban lehetne maxot is venni
			m_h.pop(iter->second->second.h_id);
			//cout << "_2Ki:" <<iter->second->first<<"\11";	
			if (iter->second->second.weight>nodes[remained].edges[iter->first]->second.weight)
			{nodes[remained].edges[iter->first]->second.weight=iter->second->second.weight;
			//cout << "_2Mo:" <<iter->second->first<<"\11";	
			if (nodes[iter->first].gap==-1)
				{m_h.move(nodes[remained].edges[iter->first]->second.h_id,-2*nodes[remained].edges[iter->first]->second.weight+nodes[iter->first].weight);
				}
			else
				{m_h.move(nodes[remained].edges[iter->first]->second.h_id,-nodes[remained].edges[iter->first]->second.weight);
				}				
			}
			if (nodes[iter->first].gap==iter->second->first+1)
				{nodes[iter->first].gap=nodes[remained].edges[iter->first]->first+1;
				}
			//a decimated ele torlodik, a remained sulya megvaltozik
			nodes[iter->first].edges.erase(decimated);
			edges.erase(iter->second->first);
			elekszama--;
			}
		else 
			{//A decimated szomszedjahoz vezeto elt atrakom es a decimatedet remainedre cserelem
			//A decimated szomszedjanak ezt az elet kitorlom a listabol es berakom ujra mint ami a remainedhez vezet
			//az el atkerul a decimatedrol a remainedre
			nodes[iter->first].edges.erase(decimated);
			if (iter->second->second.start==decimated)
				{iter->second->second.start=remained;
				}
			else
				{iter->second->second.end=remained;
				}
			nodes[remained].edges[iter->first]=iter->second;
			nodes[iter->first].edges[remained]=iter->second;
			}
		}	
	nodes[decimated].edges.clear();
	//megkeressuk a pont maximalis elet, ha van nagyobb
	kappa=nodes[remained].weight;
	tempo=-1;
	for (iter=nodes[remained].edges.begin(); iter!=nodes[remained].edges.end(); iter++)
		{if (iter->second->second.weight>=kappa)
			{tempo=iter->first;
			kappa=iter->second->second.weight;
			}
		}
	if (tempo>-1)	
		{//beallitjuk a pont gapjet
		//cout << "van nagyobb ujel:"<< tempo<<"\n";
		nodes[remained].gap=nodes[remained].edges[tempo]->first+1;
		//talaltunk nagyobb elet, ha a masik vegpontjanak is ez a maximuma
		if (nodes[tempo].gap==nodes[remained].gap)
			{kovetkezo=nodes[remained].gap-1;	
			//nodes[remained].q=1+nodes[remained].edges[tempo]->first;
			//ki kell szedni a heapbol:
			m_h.pop(edges[kovetkezo].h_id);
			//cout << "_3Ki:" <<kovetkezo<<"\11";	
			}
		//frissiteni kell a heapben 
		m_h.move(nodes[remained].h_id,-nodes[remained].weight);
			
		}
	else
		{//a pont lett maximalis
		midway--;
		//ki kell szedni a heapbol a pont kulso teret:
		m_h.pop(nodes[remained].h_id);
		hinsert(remained,fne);
		}
	nodes.erase(decimated);
	return (kovetkezo);	
	};	
	
void Jdec(int melyikel, ofstream *fneuv)
	{//addig hivja meg Jkillt, amig csak lehet J-t decimalni lokalisan
	int kelle=melyikel;
	//jmaxok.erase(melyikel);
	//cout << "-";	
	while ((kelle>=0) && (midway>1))
		{melyikel=kelle;
		kelle=Jkilld(melyikel,fneuv,edges[melyikel].weight);
		}
		
	};	

void Jdecd(int melyikel, ofstream *fneuv, double newj)
	{//addig hivja meg Jkillt, amig csak lehet J-t decimalni lokalisan
	//int kelle=melyikel;
	//jmaxok.erase(melyikel);
	//cout << "<";	
	voltej=1;	
	int kelle=Jkilld(melyikel,fneuv, newj);
	
	while ((kelle>=0) && (midway>1))
		{melyikel=kelle;
		//cout << "W";
		kelle=Jkilld(melyikel,fneuv,edges[melyikel].weight);
		}
		
	};	

void Jspread(int index, double neww)
	{//midway pont dijkstras terjesztese (dspread helyett):
	//neww: a heapben ilyen sullyal volt az el
	//ilyenkor biztos, hogy decimated lokmax, mert kulonben Jkilld-t hivtuk volna meg!	
	//int tempo=-1;
	double kappa=0;
	//omega=neww;
	//midway--;
	int remained, decimated;
	if (nodes[edges[index].start].gap>-1)
		{remained=edges[index].start;
		decimated=edges[index].end;
		}
	else
		{decimated=edges[index].start;
		remained=edges[index].end;
		}		
	//Tortenetek egyesitese:
	iteredge iter;
	pontokszama--;
	edges.erase(nodes[decimated].edges[remained]->first);
	nodes[decimated].edges.erase(remained);
	nodes[remained].edges.erase(decimated);
	elekszama--;
	///////
	for (iter=nodes[decimated].edges.begin(); iter!=nodes[decimated].edges.end(); iter++)
		{//iter->first	az a masik vegpont!
		//Ki kell szamolni az uj leendo el sulyat vagy hosszat, es csak akkor kell valtoztatni, ha az lenne a jobb:
		//kappa:
		kappa=neww-2*iter->second->second.weight+nodes[decimated].weight;		
		if (nodes[iter->first].gap==-1)		
			{kappa+=nodes[iter->first].weight;
			}
		else
			{//ha midway iter->first, akkor el volt eddig teve a heapbe a decimatedhez jovo ele, ezt most ki kell venni onnan:
			m_h.pop(iter->second->second.h_id);	
			//cout << "_4Ki:" <<iter->second->first<<"\11";	
			}
		if (nodes[remained].edges.find(iter->first) !=nodes[remained].edges.end())
			{//r megfelelo elenek sulya modosul	altalaban lehetne maxot is venni
			//m_h.pop(iter->second->second.h_id);
			if (nodes[iter->first].gap==-1)
				{//ha az uj el jobb lenne:
				if (kappa<(-2*nodes[remained].edges[iter->first]->second.weight+nodes[iter->first].weight))
					{m_h.move(nodes[remained].edges[iter->first]->second.h_id,kappa);
					nodes[remained].edges[iter->first]->second.weight=(nodes[iter->first].weight-kappa)/2.0;
					//cout << "_5Mo:" <<nodes[remained].edges[iter->first]->first<<"\11";
					}
				}
			else
				{//Ha midway az iter->first: kappa=-2*potencialissuly
				if (kappa<-2*nodes[remained].edges[iter->first]->second.weight)
					{m_h.move(nodes[remained].edges[iter->first]->second.h_id,kappa/2.0);
					//cout << "_6Mo:" <<nodes[remained].edges[iter->first]->first<<"\11";
					nodes[remained].edges[iter->first]->second.weight=-kappa/2.0;
					}
				}	
				
			if (nodes[iter->first].gap==iter->second->first+1)
				{nodes[iter->first].gap=nodes[remained].edges[iter->first]->first+1;
				}
			//a decimated ele torlodik, a remained sulya megvaltozik
			nodes[iter->first].edges.erase(decimated);
			edges.erase(iter->second->first);
			elekszama--;
			}
		else 
			{//A decimated szomszedjahoz vezeto elt atrakom es a decimatedet remainedre cserelem
			//A decimated szomszedjanak ezt az elet kitorlom a listabol es berakom ujra mint ami a remainedhez vezet
			//az el atkerul a decimatedrol a remainedre
			nodes[iter->first].edges.erase(decimated);
			if (iter->second->second.start==decimated)
				{iter->second->second.start=remained;
				}
			else
				{iter->second->second.end=remained;
				}
			nodes[remained].edges[iter->first]=iter->second;
			nodes[iter->first].edges[remained]=iter->second;
			/*nodes[remained].edges[iter->first]->second.h_id=m_h.push(kappa,nodes[remained].edges[iter->first]->first);
			*/
			//cout << "_6Be:" <<nodes[remained].edges[iter->first]->first<<"\11";
			if (nodes[iter->first].gap==-1)
				{//el kell tenni a heapbe es frissiteni kell a sulyat!
				//nodes[remained].edges[iter->first]->second.h_id=m_h.push(kappa,nodes[remained].edges[iter->first]->first);
				nodes[remained].edges[iter->first]->second.weight=(nodes[iter->first].weight-kappa)/2.0;
				nodes[remained].edges[iter->first]->second.h_id=m_h.push(kappa,nodes[remained].edges[iter->first]->first);
				}
			else
				{//ha iter->first midway: kappa=-2*potencialis suly
				//el kell tenni a heapbe es frissiteni kell a sulyat!
				//nodes[remained].edges[iter->first]->second.h_id=m_h.push(kappa,nodes[remained].edges[iter->first]->first);
				nodes[remained].edges[iter->first]->second.weight=-kappa/2.0;
				nodes[remained].edges[iter->first]->second.h_id=m_h.push(kappa/2.0,nodes[remained].edges[iter->first]->first);
				}
			}
		}	
	nodes[decimated].edges.clear();
	nodes.erase(decimated);
	};	
		
void hdecSK(ofstream *fneu)
	{//ez a resz hivatott eldonteni, hogy milyen modon decimaljunk
	int x; 
	double cella;
	//cout<< "DEC"<< midway<< "\11";
	while ((!m_h.isEmpty()) && (midway>1)) 
		{m_he = m_h.pop();
		cella=m_he.first;
		x=m_he.second;
		//cout << "!"<< m_he.first<<" ";
		if (x<0)
			{//a maximum egy h kulso ter, amibol lok. maxot kell csinalni + entropiaszamolas:
			//cout << "H"<< -x-1<<"\11"<< cella<< "\11";
			midway--;
			/*if (nodes.find(-x-1)==nodes.end()) {cout<< "Herror\n";}
			else
			{
			if (nodes[-x-1].gap<0) {cout<< "Error-H\n";}
			}*/
			hinsert(-x-1,fneu);
			}
		else
			{//elt szedtunk ki:
			//cout << "_HKi:" <<x<<"\11";
			if (edges.find(x)!=edges.end())
				{//vagy csak az egyik vege midway, vagy mindketto:
				if ((nodes[edges[x].start].gap>-1) && (nodes[edges[x].end].gap>-1))
					{//jdec:
					//cout << "J"<< x<< "s"<< edges[x].start<< "e"<<edges[x].end<<"\11";						
					Jdecd(x,fneu,-cella);
					}
				else
					{//jspread:
					//cout << "S"<< x<< "s"<< edges[x].start<< "e"<<edges[x].end<<"\11";
					Jspread(x,cella);	
					}				
				}
			}
		//egyelore nem keletkezhet az RG soran uj jmax!
		/*if (jmaxok.size()>0)
			{for (halmaz::iterator iti=jmaxok.begin();iti!=jmaxok.end();iti++)
				{//cout << "_";
				Jdec(*iti,lkiv,flogv,fneu);
				}
			jmaxok.clear();
			}*/
		}
	};	
	
int main(int argc, char *argv[])
	{char f[255];
	int nov=8;
	ido=(int) time(NULL);
	if(argc<nov)
		{//Hibas inditas eseten hibauzenet
     	};
	double magn=0;
	d = atoi(argv[1]);
	hmin = atof(argv[2]);
	double hmaxmax = atof(argv[3]);
	//Jmin = atof(argv[4]);
        Jmin=0;
        Jmax=1;
	//Jmax = atof(argv[5]);
	Lmin = atoi(argv[4]);
    Lmax=Lmin;
	//Lmax = atoi(argv[7]);
	graphform = atoi(argv[5]);
	//graphform (pbc chain=0) (ladder...plain=1) (plain...cubic=2) (clique=3)
	rep = atoi(argv[6]);
	//betha = atof(argv[10]);
	//method = atoi(argv[11]);
	//Határfeltételt szabályozza: 0=pbc 1=free 
	//s = atoi(argv[12]);
	//a pontossag hanyad resze legyen a beirt szorasnak: javaslat: 25-od resze, de a beirt szoras legyen a velhetoen igazinak a fele
	//igy az erdekes tartomanyt 100 reszre osztottuk fel
	//int mennyibe=atoi(argv[13]);
	//kmethod=atoi(argv[15]);
	kiiratas=atoi(argv[7]);
	
	muin=1;
	int fend=0;
	char random_seed; 
	
	if(argc>nov)
		{fend = atoi(argv[nov]);
	//cout << "OUT";
		}
	else
		{FILE *random_file = fopen("/dev/urandom","r");
	//NEEDS ERROR CHECKING
	if (!random_file)
		{cout << "Zaj fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	random_seed = getc(random_file);
	fend=random_seed;
		}
		MTRand mtrand(fend*ido);
	

	char sta[255];	
	sprintf(sta,"stat_%d_%d_%d_%d_%d.dat",d,Lmin,rep,fend,ido);
	ofstream lsta(sta, ios::out);
	if (!lsta)
		{//cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};	
	char mag[255];	
	sprintf(mag,"mag_%d_%d_%d_%d_%d.dat",d,Lmin,rep,fend,ido);
	ofstream lmag(mag, ios::out);
	if (!lmag)
		{cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	char neu[255];	
	sprintf(neu,"vN_%d_%d_%d_%d_%d.dat",d,Lmin,rep,fend,ido);
	ofstream lneu(neu, ios::out);
	//Ebbe mehet az utolso klaszter merete:
	if (!lneu)
		{cout << "Kimeneti entropia fajl megnyitasa sikertelen!\n";
		exit(-1);
		};
	char kiiro[255];		
	sprintf(kiiro,"KIO_%d_%d_%d_%d_%d.dat",d,Lmin,rep,fend,ido);
	ofstream lkiiro(kiiro, ios::out);
	if (!lkiiro)
		{cout << "Kimeneti fajl megnyitasa sikertelen!\n";
		exit(-1);
		};	
	// A pontokat 0-tol kezdjuk el szamozni
	NODE_START_ID=0;
	Nodes=0; Edges=0;
	//logo << (int) time(NULL)<< "\n";
int kezdi=(int) time(NULL);
	double szoras=0,varhato=0;
	for(int meret=Lmin; meret<=Lmax;/*meret++*/)
		{lsta << meret << '\11';
		varhato=0;
		szoras=0;
		for(realizacio=0; realizacio<rep; realizacio++)
			{
//cout <<"R"<<realizacio <<"\11";
			while (!m_h.isEmpty()) 
				{//TT
				m_h.pop();	
				}
			vent.clear();
			duplae=0;
			lastj=log(Jmax);
			voltej=0;
			magn=0;
			//grafparameterek(meret, d);		
			perk=0;
			//cella=0;
			//omega=exp(hmaxmax);
			if (graphform>=2)
				{d=meret;
				}
			grafparameterek(meret, d);		
			minweight=-INF;
			minsize=0;
			midway=0;
			elekbehuzasa(meret, d);
                //cout << "\n";
			for(int i=0; i<Edges; i++)
				{edges[i].weight = log(Jmin+(Jmax-Jmin)*(1-mtrand.rand()));
                  //  cout << edges[i].weight << "\11";
				}
               // cout << "\n";
			for(int i=0; i<Nodes; i++)
				{nodes[i].mu=muin;
				if (hmin==0)
					{nodes[i].weight = hmaxmax+log(1-mtrand.rand());
                      //  cout << nodes[i].weight << "\11";
					}
				else
					{if (hmaxmax>hmin)
						{nodes[i].weight = log(hmin+(exp(hmaxmax)-hmin)*(1-mtrand.rand()));
						}
					else
						{nodes[i].weight = hmaxmax;
						}
					}
				}
               // cout << "\n";
			omega=max(log(Jmax), hmaxmax);
			EDGE_START_ID=Nodes;
			pontokszama=Nodes;
			elekszama=Edges;
			duplazas();
			if (graphform==0)
				{cross1dperem(meret);
				}
			if (graphform==2)
				{cross2dperem(meret,d);
				}
			if (graphform==3)
				{cross3dperem(meret);
				}	
			if (graphform==4)
				{cross4dperem(meret);
				}					
			duplae=1;
			pontbesoroloperem(log(Jmax));	
			//cout << "." << midway << "\11";	
			//pajek(0);
			if (jmaxok.size()>0)
				{//cout << ".j"<< jmaxok.size()<< "\11";
				voltej=1;
				for (halmaz::iterator iti=jmaxok.begin();iti!=jmaxok.end();iti++)
					{
					Jdec(*iti,&lneu);
					}
				jmaxok.clear();
				}
			hdecSK(&lneu);
			lastminute(&lneu);
			magn=0;	
		if (perk>0)
			{minsize=div(minsize,2).quot;	
			magn=minsize;
			}
		if ((lastj<minweight) && (voltej>0))
			{minweight=lastj;
			//if (kiiratas>=2) {lkiir << minweight<< '\11';}
			}
		else
			{//if (kiiratas>=2) {lkiirp << minweight<< '\11';}
			}
		if (kiiratas>=3) {lkiiro << minweight<< '\11';}
		if (kiiratas>=2)
			{if (perk>0)
				{lmag << minsize << '\11';
				}
			else
				{lmag << minsize << '\11';
				}
			}
		varhato+=magn;
		szoras+=minsize;
		}
		if (kiiratas>=2) {lmag << '\n';
           // lperk << '\n';
            
        }
		if (kiiratas>=3) {//lkiirp << '\n';lkiir << '\n';
            lkiiro << '\n';}
		varhato/=rep;
		szoras/=rep;
		//A szoras helyett a legnagyobb klaszter meretenek atlagat irjuk ki, hogy lassuk a kulonbseget:
		lsta << setprecision(12) << varhato << '\11' << szoras << '\11'<< rep << '\n';
		meret*=2;
		}
	//logo << (int) time(NULL)<< "\11" << Edges <<"\n";	
//logo << Lmin << "\11"<<((double) ((int) time(NULL)-kezdi))/((double) rep) <<"\n";
	nodes.clear();
	edges.clear();
	//logo.close();
	lneu.close();
	lsta.close();
	lmag.close();
	return 0;
	}
