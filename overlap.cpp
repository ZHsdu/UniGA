#include "main.h"
#include "overlap.h"
#include "consense.h"
#include "spoa.hpp"
#define threshold 10
#define THRESHOLD 656
#define deep 40
#define exact 200
#define rth 1
ofstream BH("seed1.fasta");
extern int HOR_LEN;
int BHNUM=0;
bool erase_flag;
vector<int> end_read;
class read_info{
public:
	int read_flag0;
	bool if_seed0;
	int ppath0 ;
	int pdirect0[2];
	bool cand_flag0;
	read_info(){};
	read_info(int a,bool b,int c,int d[],bool e)
	{
		read_flag0=a;
		if_seed0=b;
		ppath0=c;
		pdirect0[0]=d[0];
		pdirect0[1]=d[1];
		cand_flag0=e;
	}
};
double **P_;
extern int real_read_size;
extern int CAND_SIZE;
int filter;
string AA;
string consensus_total;
int mxseed;
int consensus_len;
int Nc;
int mxseed_locat;
int err_repeat=0;
int sig_end=1;
int end_id=-1;
int correct_start=0;
extern bool is_re;
extern int _deep;
pair<int,int> **num_;
int *n_rare_;
int *palign;
string _buf;
string _buf0;
bool correct_flag=false;
bool *if_seed;
bool seed_err=false;
bool seed_is_correct;
	int* ppath;
	int** pdirect;
vector<int> readtem;
vector<int> readtem2;
vector<int> readchange;
vector<int> readchange2;
vector<int> creadtem;
vector<int> _P;
vector<inff> lastcand;
vector<pair<int,int>> _lastcand ;
vector<int> lbra;
vector<int> rbra;
vector<thread> threads;
bool unrepeat=true;
bool thread_exit=false;
bool thread_snp_exit=false;
bool correct_suc=false;
int target=0;
int coverage[2]={0};
int extendsize[2]={0};
ofstream fP("cand.info");
bool end_flag_1,end_flag_2;
int cal(char c, int i,vector<char> &fla_a )
{
	while ( i < fla_a.size() && fla_a[i] == c )
		i++;
	return i;
}
int maximal(vector<rdid> x)
{
	int max = 0;
	for (const auto& vv : x)
	{
		if (vv.self_locat > max) max = vv.self_locat;
	}
	return max;
}
int minimal(vector<rdid> x)
{

	int min = 0xffffff;
	for (const auto& vv : x)
	{
		if (vv.self_locat < min) min = vv.self_locat;
	}
	return min;
}
class Path_
{
	public:
		int id;
		int len;
		vector<int> cand;
	void push(vector<inff> &Cand,inff &nextseed,int Right)
	{
		for (const auto& v:Cand)
			if (v.is_rextend||v.is_insert||v.is_lextend)
				cand.push_back(v.id);
		id=nextseed.id;
		len=nextseed.extend[Right];
	}
};
/*
void num_update(pair<int,int> num_[])
{
	double x;
	double l=200;
	double m=200;
	for (int i=0;i<50;i++)
		num_[i]=make_pair(0,0);
	for (int i=13;i<=25;i++)
	{
		x=pow(0.85,i-1)*0.15/(3.0*i)*l*m;
		num_[i].first=ceil(x);
	}
	for (int i=13;i<=25;i++)
		num_[i].second=num_[i].first+(pow(0.85,i)*40);
}
*/
void num_update(pair<int,int> num_[],int M,int kk,int d)
{
	double x;
	double p;
	double N=M,N0=kk*d;
	int t=2;
	for (int i=0;i<50;i++)
		num_[i]=make_pair(0,0);
	for (int i=k_len;i<=k_lenmax;i++)
	{
		double kp=0.12;
		double pr=1-kp;
		double kb=kp*i;
		double kd=sqrt(kp*(1-kp)*i);
		double yp=pow(pr,i);
		double pu=0;
//		for (int j=1;j<4;j++)
//			pu+=pow(pr,i-j)*pow(kp,j)/(pow(3,j));
		for (int j=1;j<4;j++)
			pu+=pow(pr,i-j)*pow(kp,j)/((pow(4,j)-1));
		P_[i][0]=yp;
		P_[i][1]=pu;
		n_rare_[i]=pu*N;
		num_[i].first=n_rare_[i]-sqrt(pu*(1-pu)*N);
		num_[i].second=yp*N0;
	}
}
/*void num_update(pair<int,int> num_[],int M,int kk,int d)
{
	double x;
	double p;
	double N=M,N0=kk*d;
	int t=2;
	for (int i=0;i<50;i++)
		num_[i]=make_pair(0,0);
	for (int i=k_len;i<=27;i++)
	{
		double kp=0.12;
		double pr=1-kp;
		double km=kp*i;
		double ku=km,kb=km;
		double kd=sqrt(kp*(1-kp)*i);
//		double ku=km+t*kd;
//		double kb=max(double(1),km-t*kd);
		double pu=pow(pr,i-kb)*pow(kp,kb)/pow(3.0,kb);
		double pb=pow(pr,i-ku)*pow(kp,ku)/pow(3.0,ku);
		double xu=pu*N+t*sqrt(pu*(1-pu)*N);
		n_rare_[i]=pu*N;
		double xb=pb*N-t*sqrt(pb*(1-pb)*N);
		double yp=pow(pr,i);
		double yu=yp*N0+t*sqrt(yp*(1-yp)*N0);
		double yb=yp*N0-t*sqrt(yp*(1-yp)*N0);
		num_[i].first=max(2,(int)(xb+yb));
		num_[i].second=ceil(xu+yu);
	}
}*/
void num_update(pair<int,int> **num_,int M,int max_rare,int d)
{
	for (int i=0;i<max_rare;i++)
		num_update(num_[i],M,i+1,d);
}
void num_init(int max_rare)
{
	num_=new pair<int,int> *[max_rare];
	n_rare_=new int[50];
	P_=new double*[50];
	for (int i=0;i<50;i++)
	{
		n_rare_[i]=0;
		P_[i]=new double[2];
	}
	for (int i=0;i<max_rare;i++)
	{
		num_[i]=new pair<int,int>[50];
	}
	for (int i=0;i<max_rare;i++)
		for (int j=0;j<50;j++)
		{
			num_[i][j].first=0;
			num_[i][j].second=0;
		}
}
void num_delete(int max_rare)
{
	for (int i=0;i<max_rare;i++)
		delete[] num_[i];
	delete[] num_;	
	delete[] n_rare_;
}
int n_rare(int l,int n,int max_rare)
{
        double score=99999,c,d;
        int j=0;
        for (int i=0;i<max_rare;i++)
        {
                double  a=num_[i][l].second,
                        b=num_[i+1][l].second;
                if (n<=a) return i+1;
                if (n>b) continue;
                if (n<=b)
                {
                        return i+1;

                }
        }
}
/*
int n_rare(int l,int n,int max_rare)
{
	double score=99999,c,d;
	int j=0;
	for (int i=0;i<max_rare;i++)
	{
		double 	a=num_[i][l].first,
			b=num_[i][l].second;
		if (n<=a) continue;	
		if (n>b) {continue;}
		if ((c=((abs(n-a)+abs(n-b))/2))<(b-a)*score)
		{
			score=c/(b-a);
			j=i+1;	
		}
	}
	return j;
}*/
void out_info(int read_flag[],bool cand_flag[])
{
	ofstream gg("info_read",ios::out | ios::binary);
	for(int i=0;i<real_read_size;i++)
	{
		read_info s(
		read_flag[i],
		if_seed[i],
		ppath[i],
		pdirect[i],
		cand_flag[i]);
		gg.write((char*)&s,sizeof(s));
	}
	gg.close();
}
void in_info(int read_flag[],bool cand_flag[])
{
	ifstream gg("info_read",ios::in | ios::binary);
	read_info s;
	int i=0;
	while(gg.read((char*)&s,sizeof(s)))
	{
		read_flag[i]=s.read_flag0;
		if_seed[i]=s.if_seed0;
		ppath[i]=s.ppath0;
		pdirect[i][0]=s.pdirect0[0];
		pdirect[i][1]=s.pdirect0[1];
		cand_flag[i]=s.cand_flag0;
		
	}
}
void out_contig(int contig_id,vector<int> &left_read,vector<int> &right_read,string &Consensus)
{
	
	ofstream ff;
	ff.open("step1",ofstream::app);
	ff<<contig_id<<"\n";
	for (auto &v:left_read)
		ff<<v<<"\t";
	ff<<"\n";
	for (auto &v:right_read)
		ff<<v<<"\t";
	ff<<"\n";
	ff<<Consensus<<"\n-\n";
	ff.close();
	ff.open("step1.fasta",ofstream::app);
	ff<<">contig_"<<contig_id<<"\n";
	ff<<Consensus<<"\n";
	ff.close();
}
void overlap(vector<vector<_pair>> &minimizer, vector<string> &read, vector<vector<pair<int, int>> >& mini_hash, unordered_map<int, vector<rdid > >& read_id, char* output)
{
	
	ofstream f(output);
	inff seed,seed0;
	vector<rdid> A;
	vector<int> Path, path,_path;
	unordered_map<int, err>error;
	end_flag_1=false;
	end_flag_2=false;
	int max_rare=8;
	num_init(max_rare);
	num_update(num_,40000,max_rare,50);
	int i = 0,read_count=real_read_size,mx,itt_=241;
	int* read_flag;
	bool *cand_flag;
	int times=0,seqs=0;
	read_flag = new int[read_count];
	if_seed= new bool[read_count];
	ppath = new int[read_count];
	pdirect=new int* [read_count];
	cand_flag=new bool[read_count+1000];
	palign=new int[real_read_size];
	for (int i=0;i<read_count;i++)
	{
		pdirect[i]=new int[2];
		pdirect[i][0]=-1;
		pdirect[i][1]=-1;
		if_seed[i]=true;
		cand_flag[i]=false;
		palign[i]=1;
	}
	for (int i=0;i<read_count+1000;i++)
		cand_flag[i]=false;
	vector<pair<int,int> > read_size_list;
	pair<int,int> _T;
	clock_t  st,ed;
	for (int u = 0; u < read_count; u++)
	{
		read_flag[u] = 1;
		ppath[u] = -1;
	}
/*	for (int i=0;i<real_read_size;i++)
	{
		if (read[i].size()>2000)
			if (if_palindrome(read[i]))
				palign[i]=0;
	}
	ofstream hh("palign");
	for (int i=0;i<real_read_size;i++)
	{
		hh<<palign[i]<<"\n";
	}
*/
	ifstream hj("../../../date/palign");
	int ii=0;
	while(!hj.eof())	
	{
		hj>>palign[ii++];
	}
/*
	for (int i=0;i<real_read_size;i++)
	{
		if (palign[i]==0)
			read_flag[i]=0;
	}*/
	for (int u=0;u<real_read_size;u++)
	{
		_T.first=read[u].size();
		_T.second=u;	
		read_size_list.push_back(_T);
	}
	sort(read_size_list.begin(),read_size_list.end(),cmp5);
/*	for (int i=0;i<read_size_list.size();i++)
	{
		if (read_size_list[i].second!=3619)
			if_seed[read_size_list[i].second]=false;
		else
			break;
	}*/
	for (int i=0;i<16;i++)
		threads.push_back(thread(run,i));
//	while(read_size_list[itt_].second!=676303&&itt_<real_read_size) 
//			itt_++;
	for (int i=0;i<16;i++)
		threads.push_back(thread(run0,i));
	int contig_id=0;
	while(1)
	{
	correct_flag=false;
	int kkll=0;
	unrepeat=true;
	seed_err=false;
	rbra.clear();
	lbra.clear();
	while(!read_flag[read_size_list[itt_].second]&&itt_<real_read_size) itt_++;
	if (itt_==real_read_size)
		break;
	i=read_size_list[itt_].second;
	if (!if_seed[i])
	{
		itt_++;
		continue;
	}
	if (read[i].size()<1000)
		break;
	times++;
	st=clock();
	_buf.clear();
	cout<<"seed:"<<i<<" size:"<<read[i].size()<<endl;
	vector<int> left_read;
	vector<int> right_read;
	end_read.clear();
	seed.in(minimizer[i], read[i], i, A,seed.str.size());
	mxseed=i;
	bool act=seed_check_update(read,seed,minimizer,mini_hash,read_id,read_flag);
	if (!act)
	{
		itt_++;
		continue;
	}
	seed0=seed;
	rf_change2(read_flag, i);
	if (!Path.empty())
		Path.clear();
	if (!path.empty())
		path.clear();
	if (!readchange.empty())
		readchange.clear();
	seed0.re_str();
	string right_Consensus;
	right_Consensus.insert(0,seed0.str);
	extend(ppath,error,seed0, read_flag, minimizer, read, mini_hash,read_id, Path, 0,pdirect,cand_flag,right_Consensus);
	right_read.swap(end_read);
	right_Consensus=rev(right_Consensus);
	if (right_Consensus.size()<seed.str.size())
	{
		itt_++;
		continue;
	}
	seed.str=right_Consensus.substr(0,seed.str.size());
	if (seed_err)
	{
		times--;
		continue;
	}
	if (!unrepeat)
	{
		times--;
		if_seed[i]=false;
		if_seed[twin(i)]=false;
		rf_back(read_flag);	
		continue;
	}
	if (!readchange2.empty())
		{
		seed.in(minimizer[mxseed], read[mxseed], mxseed, A,0);
		for (const auto& vv:readchange2)
			rf_back(read_flag,vv);
		cout<<"!!----------------------------------!!\n\tstart_change:"<<mxseed<<"\n";
		readchange2.clear();
		}
	for (i = Path.size() - 1; i >= 0; i--)
		path.push_back(Path[i]);
	if (!Path.empty())
		Path.clear();
	mxseed_locat=path.size();
	path.push_back(mxseed);
	//i=mxseed;
	//seed.in(minimizer[i], read[i], i, A,seed.str.size());
	cout<<"_____________------------____________"<<endl;
	_buf+="----------Right=0-------------\n";
	_lastcand.clear();
	//seed.re_str();
	string left_Consensus;
	extend(ppath,error,seed, read_flag, minimizer, read, mini_hash, read_id, Path, 0,pdirect,cand_flag,left_Consensus);
	left_read.swap(end_read);
	string Consensus=left_Consensus+right_Consensus;
	out_contig(contig_id,left_read,right_read,Consensus);
	contig_id++;
	_buf+="----------Right=1-------------\n";
	_buf+="end________________________\n";
	for (i = 0; i < Path.size(); i++)
		path.push_back(Path[i]);
	cout<<"	len:"<<path.size();
	ed=clock();
	cout<<" time:"<<(double)(ed-st)/CLOCKS_PER_SEC<<endl;
	if (!Path.empty())
		Path.clear();
	_P=path;
//	correct_flag=true;
//	path_correct(ppath,error, read_flag, minimizer, read, mini_hash, read_id,path,0);
	out_info(read_flag,cand_flag);
	if (path.size()<20) 
	{
		seqs++;
		readtem2.clear();
		continue;
	}
	else 
		seqs=0;
	cf_change(readtem2,cand_flag);
/*	if (!readtem2.empty())
	{
		for (auto vv:readtem2)
			rf_change2(read_flag,vv);
		readtem2.clear();
	}*/
	f << "begin:\n\nPath:";
	for (auto v : path)
	{
		if (v > read.size() / 2) v = read.size() / 2 - v;
		f << v << "->";
	}
	f << "\n" << endl;
	fP<<_buf;
	consensus_total.clear();
	rf_change2(read_flag,mxseed);
	}
	thread_exit=true;
	thread_snp_exit=true;
	start_the();
	start_snp_the();
	for(int i=0;i<32;i++)
		threads[i].join();
	for (int i=0;i<real_read_size;i++)
	{
		if(read_flag[i]==1)
			f<<i<<":"<<read[i].size()<<"\n";
	}
	num_delete(max_rare);
}
	
void path_correct(int* ppath,unordered_map<int, err>& error, int* read_flag, vector<vector<_pair> >& minimizer, vector<string>& read, vector<vector<pair<int, int > > > & mini_hash, unordered_map<int, vector<rdid > >& read_id, vector<int>& path,int Right,bool cand_flag[])
{
	int err_=0,I,Check=-1,i=0,count_=0,count=0,p_locat=0;
	vector<int>::iterator _pa,_pb;
	vector<int> Path;
	vector<Path_>pPath;
	vector<int> path_insrt;
	vector<pair<int,int>> err_count;
	pair<int,int> t;
	inff seed;
	vector<inff> cand0;
	vector<rdid> A;
	correct_start=1;
	while(1)
	{
	if (!err_count.empty())
		err_count.clear();
	for (int i=0;i<path.size();i++)
	{	
		if (error[path[i]].repeat>=1)
		{
			count=0;
			for (auto v:error[path[i]].clique)
			{
				count_=1;
				for (auto vv:v)
					if (read_flag[vv]==0)
					{
						count_=0;
						break;
					}
				count+=count_;
			}
			t.first=count;
			t.second=i;
			if (count>=1)
				err_count.push_back(t);
		}
	}
	sort(err_count.begin(),err_count.end(),cmp5);
	if (Check==err_count.size()) break;
	Check=err_count.size();
	for (int ii=err_count.size()-1;ii>=0;ii--)
	{
		if (err_count[ii].first==0) 
			break;
		i=err_count[ii].second;
		if (i>=mxseed_locat)
			Right=1;
		else Right=0;
		I=path[i];
		end_id=-1;
		int kk=0,kt=0,tag=0,j=0,_flag=0,Len=0,Len_=-1;;
		if (error[I].repeat>0)
		{
			for (j=0;j<error[I].repeat;j++)
			{	
				_flag=0;
				kt=0;
				kk=0;
				
				for (auto vv:error[I].clique[j])
				{	
				if (read_flag[vv]==1) 
					kk++;
				kt++;
				}
			
				if (kk!=kt)
				{
//					for (auto vv:error[I].clique[j])
//						rf_change(read_flag,vv);
//					error[I].repeat--;
					continue;
				}
    				for (auto v:error[I].clique[j])
				{
					if (read_flag[v]==1)
						if(Len<read[v].size())
						{	
							Len=read[v].size();
							Len_=v;
						}
				}	
				if (!creadtem.empty())
					creadtem.clear();
				_flag=0;
				seed.in(minimizer[Len_], read[Len_], Len_, A,Len);
				if (!Path.empty())
					 Path.clear();
				if (!readtem.empty())
					readtem.clear();
				readtem.push_back(Len_);
				Path.push_back(seed.id);
				mxseed=I;
				correct_suc=false;
				string Consensus;
				extend(ppath,error,seed, read_flag, minimizer, read, mini_hash,read_id, Path, Right,pdirect,cand_flag,Consensus);
				if (Path.size()>=10&&Path.back()!=-1)
					for (auto vv:path)
					{
						if (ppath[Path.back()]==vv)
							_flag=1;
					}
				else
					Path[Path.size()-1]=Path[Path.size()-2];
				if (_flag!=1)
				{
					pp_back(ppath,pdirect);
					for (auto v_:readtem)
						rf_back(read_flag,v_);
					 continue;
				}
				if (correct_suc)
				{
					_pa=path.begin();
					advance(_pa,i);
					if (Right)
					{
					path.insert(_pa, Path.begin(), Path.end());
					}
					else
					{
					reverse(Path.begin(),Path.end());
					path.insert(_pa+1, Path.begin(), Path.end());
					}	
					
				}		
				
		/*		p_locat=0;
				for (_pa=path.begin();_pa<path.end();_pa++)
				{
					
					if (*_pa == ppath[Path.back()])
						break;
					p_locat++;
				}
				if (_pa == path.end()||p_locat<5||path.size()-p_locat<5)
				{
					cout << "error";
					Path.clear();
				 	for (auto v_:readtem)
						rf_back(read_flag,v_);	
					if (!readtem2.empty())
						readtem2.clear();
					pp_back(ppath,pdirect);
					continue;
				}
				_pb=path.begin();
				advance(_pb,i);
				if (Right==0)
				{
					if (distance(_pb,_pa)<0)
					{
						for (auto v_:readtem) 
                                                       	 rf_back(read_flag,v_);
					if (!readtem2.empty())
						readtem2.clear();
						Path.clear();
						pp_back(ppath,pdirect);
						continue;
					}
					for (int k=i;path[k]!=*(_pa+1);k++)
						path_insrt.push_back(path[k]);
					reverse(Path.begin(),Path.end());
					_pa=path.begin();
					advance(_pa,i);
				}
				else
				{
					if (distance(_pa,_pb)<0)
                                        {
						for (auto v_:readtem)
                                               		 rf_back(read_flag,v_);
					if (!readtem2.empty())
						readtem2.clear();
                                                Path.clear();
						pp_back(ppath,pdirect);
                                                continue;
                                        }
					for (vector<int>::iterator p1=_pa;*p1!=path[i+1];p1++)
						path_insrt.push_back(*p1);
				}
				path_insrt.insert(path_insrt.end(),Path.begin(),Path.end());
				path.insert(_pa, path_insrt.begin(), path_insrt.end());*/
				Path.clear();
				path_insrt.clear();
				for (auto vv_:error[I].clique[j])
					rf_change2(read_flag,vv_);
				if (!readtem2.empty())
				{
					for (auto vv:readtem2)
						rf_change2(read_flag,vv);
					readtem2.clear();
				}
			}
		
		}
	}
	}

}
void extend(int* ppath,unordered_map<int, err>& error,inff seed, int* read_flag, vector<vector<_pair> >& minimizer, vector<string>& read, vector< vector<pair<int, int > > > & mini_hash, unordered_map<int, vector<rdid > >& read_id, vector<int> &Path, int Right,int ** pdirect,bool cand_flag[],string &Consensus)
{
	int i,count=0;
	vector<table> t;
	vector <int> snp_sum;
	vector<string> candidate;
	vector<vector<_pair>> candidate_mini;
	vector<vector<pair<char, char> > > alnment;
	vector<Path_> pPath;
	unordered_map<int, vector<rdid > > read_temporary;
	vector<rdid> A;
	rdid p;
	inff x;
	target=0;
	inff nextseed;
	vector<inff> cand;
	vector<inff> candb;
	vector<inff> candtem;
	vector<inff> cand0;
	vector<INFO> align_info;
	stack<pair<int,int>> branch;
	pair<int,int> _x;
	vector<char> fla_a;
	vector<diver> read_set;
	vector<int> path;
	void lastcand_update(int* read_flag);
	int RE=0;
	int C,next=-1,next_=-1,c_len=0,candc_len=0;
	int FLG,num_;
	bool if_end=false;
	Path_ _p;
//	alnment.clear();
	int _maxlen=0,_end=0,K_,len,tag;
	char c_str[1000000];
	vector<Path_>::iterator ppo;
	ofstream ff("close");
	ofstream ff5("close_num");
	bool repeat = 0;
	if (_lastcand.empty())
	_lastcand.clear();
	end_flag_1=false;
	end_flag_2=false;
	bool _start=true;
	int palign_target=0;
	while (1)
	{
		if (_start)
		{
		//		seed_rchose(seed,523);
				//seed_rchose(seed,2174);
	//		seed_rchose(seed,1075);
			_start=false;
		}
		bool jp=false;
		_buf0.clear();	
		seed_is_correct=true;
		next_ = -1;
		count++;
		RE=0;
		c_len=seed.c_len;
		end_flag_1=false;
		bool mark=false;
		int ec_len=0;
		if (!candtem.empty())
			candtem.clear();
		if (!read_id.empty())
			read_id.clear();
		if (1)
		{
			inff seedt=seed;
			cand_select(read,seed,candtem,mini_hash,read_id,read_flag);
			_mini_map(seed,candtem,Right,read,minimizer,true,Consensus);
			if (seed.is_consensus==true)
			{
//				seed=seedt;
//
//
//
//
//				centromare//
//
				seed.is_consensus=true;
				candtem.clear();
				read_id.clear();
				cand_select(read,seed,candtem,mini_hash,read_id,read_flag);
				_mini_map(seed,candtem,Right,read,minimizer,false,Consensus);
				vector<inff> cand0;
				if (seed.is_consensus)
				{
					for (auto &v:candtem)
						end_read.push_back(v.id);
				}
				seed.is_consensus=false;
				BH<<">read_"<<BHNUM++<<"\n"<<seed.str<<"\n";
				return ;
				}
		}
		else 
		{
			palign_target++;
			seed_is_correct=false;
		}
		int score=0;
		palign_target=0;
		int _i;
		_i=candtem.size();
		int f_;
/*		for (auto v:lastcand)
		{
			f_=1;
			for (int ii=0;ii<_i;ii++)
			{
				if (v.id==candtem[ii].id)
				{
					f_=0;
					break;
				}
			}	
			if (f_)
			{
				v.if_cand=0;
				v.if_extend=0;
				candtem.push_back(v);
			}
		}*/
	//	for (auto v:candtem)
	//	{
	//		if_seed[v.id]=false;
	//		if_seed[twin(v.id)]=false;
	//	}
//		for (int i=0;i<_lastcand.size();i++)
//		{
//			if (_lastcand[i].first<real_read_size)
//				rf_change(read_flag,_lastcand[i].first);
//		}
/*		if (seed.id!=mxseed&&_lastcand.size()<50)
		{
			for (const auto v:candtem)
				if (v.is_rextend&&!v.is_lextend)
				{
					rf_change2(read_flag,v.id);
				}
		}
*/		_lastcand.clear();
		if (!correct_flag)
		for (const auto& v:candtem)
		{
			if (v.id==seed.id)
				continue;
			if_seed[v.id]=false;
			if_seed[twin(v.id)]=false;
			//if (Right==0&&v.id==mxseed)
			if (Right==0&&v.id==mxseed&&v.is_lextend)
			{
	//		int opl=0;
	//		opl++;	
				unrepeat=false;
			cout<<"--------new----------"<<endl;
				return;
			}
		}
		if (correct_flag)
		{
			if (Right)
			{
				for (const auto& v:candtem)
				{
					if (v.is_rextend&&v.id==mxseed)
					{
						correct_suc=true;
						
						return;
					}
				}
			}
			else
			{
				for (const auto& v:candtem)
				{
					if (v.is_lextend&&v.id==mxseed)
					{
						correct_suc=true;
						return;
					}
				}
			}
		}
/*		if (!correct_flag&&seed.id==mxseed&&Right==0&&(coverage[0]>_deep||coverage[1]>_deep))
		{
			unrepeat=false;
			cout<<"--------new0----------"<<endl;
			return;
		}*/
		if (Right)
		{
			if (coverage[0]-extendsize[1]>_deep)
				mark=true; 
			else 
				mark=false;
		}
		else 
		{
			if (coverage[1]<-extendsize[0]>_deep)
				mark=true;
			else 
				mark=false;
		}
		static int now_num=0;
		vector<rdid> AA;
		cand0.clear();
		cout<<"\rseed:len:"<<Consensus.size()<<"        "<<now_num++<<flush;
		vector<inff> cand_;
		if (correct_flag==true)
			cand.clear();
		for (const auto& v:candtem)
		{
			if (v.is_lextend||v.is_rextend||v.is_insert)
				pp_change(ppath,v.id,seed.id);
		}
			
		if_seed_change(candtem);
		for (const auto& v:candtem)
            	{
			if (v.id>=real_read_size)
				continue;
			if (ppath[v.id]==seed.id)
			{
				if (v.extend[0]<200&&v.extend[1]>200)
				{
					if (pdirect[v.id][0]==-1)
					{
					pdirect[v.id][0]=seed.id;
					pdirect[twin(v.id)][0]=seed.id;
					}
				}
				if (v.extend[1]<200&&v.extend[0]>200)
				{
					if (pdirect[v.id][1]==-1)
					{
					pdirect[v.id][1]=seed.id;
					pdirect[twin(v.id)][1]=seed.id;
		     			}
				}
			}
		}
		int branch_count=0;
		if (!candb.empty())
			candb.clear();
                len=candtem.size();
		for (int i=0;i<candtem.size();i++)
			if (candtem[i].if_cand==1)
				cand.push_back(candtem[i]);
		for (const auto& v:cand)
			if (v.id!=seed.id&&v.if_extend==1)
				cand_.push_back(v);
		for (const auto& v:cand)
		{
			if_seed[v.id]=false;
			if_seed[twin(v.id)]=false;
		}
		int done=0;
		done=snp_find(seed, candtem, alnment, read_flag, ff, read_set, ff5,Right,ec_len,Consensus);
		if (done==-1)
			break;
		seed.is_consensus=false;
		consensus_len=0;
		seed.confirm=0;
		seed.cons_len.clear();
		Consensus.insert(0,seed.str);
		int len=10000;
		seed.str=Consensus.substr(0,len);
		seed.row_str=seed.str;
		BH<<">read_"<<BHNUM++<<"\n"<<seed.str<<"\n";
		seed.msa.clear();
		seed.cov.clear();
		_p.push(candtem,nextseed,Right);
		pPath.push_back(_p);
		if (!read_id.empty())
			read_id.clear();
		if (!candtem.empty())
			candtem.clear();
		if (!alnment.empty())
		{
			for (auto v:alnment)
				vector<pair<char, char> >().swap(v);
			alnment.clear();
		}	
		if (!cand.empty())
			vector<inff>().swap(cand);
		if (!t.empty())
			for (auto v:t)
				v.cle();
		
	}
}
bool cmp1(rdid x, rdid y)
{
	return x.d < y.d;
}
bool cmp2(rdid x, rdid y)
{
	return x.seed_locat < y.seed_locat;
}
bool cmp3(inff x, inff y)
{
	return x.match_start < y.match_start; 
}
bool cmp4(inff x, inff y)
{
	return x.match_end < y.match_end;
} 
bool cmp5(pair<int,int> x, pair<int,int> y)
{
	return x.first > y.first; 
}
bool cmp6(pair<int,int> x,pair<int,int>y)
{
	return x.second<y.second;
}
bool cmp7(pair<int,int> x,pair<int,int>y)
{
	return x.second>y.second;
}
bool cmp8(int x,int y)
{
	return x>y;
}
bool cmp9(int x,int y)
{
	return x<y;
}
bool cmp10(pair<int,int> x, pair<int,int> y)
{
	return x.first < y.first; 
}
bool cmp14(_pair x, _pair y)
{
	return x.second < y.second; 
}
bool cmp15(pair<int,double> x, pair<int,double> y)
{
	return x.second < y.second; 
}
bool cmp16(pair<int,double> x, pair<int,double> y)
{
	return x.second > y.second; 
}

int find_count(vector< rdid> &T, vector<rdid>& A)
{
	int m_=656,m_kmer_num=0,max_count=0,max_len=0,len=0;
	int d=T[0].d;
	len=T[0].self_locat;
	for (int i=1,size=T.size();i<size;i++)
	{
		if (T[i].d-d<m_)
			++m_kmer_num;
		else
		{
			len=T[i].self_locat-len;
			max_len=max(max_len,len);
			max_count=max(m_kmer_num,max_count);
			m_kmer_num=0;
			len=T[i].self_locat;
		}
	}
	len=T.back().self_locat-len;
	max_len=max(max_len,len);
	max_count=max(m_kmer_num,max_count);
	if (max_len>500&&max_count>2)
		return 1;
	return 0;
}
void cand_check2(inff seed, vector<inff>& cand,vector<inff>& cand_2left,vector<inff>& cand_2right,vector<inff>& candleft,vector<inff>& candright)
{
	int len = cand.size();
	int* K = new int[len];
	for (int i = 0; i < len; i++)
		K[i] = i;
	sort(cand.begin(), cand.end(), cmp3);
	for (int i = 0; i < len; i++)
	{
		if (cand[i].str.size()<4000) continue;
		if (cand[i].match[0].seed_locat < 1000&&cand[i].match.back().seed_locat-cand[i].match[0].seed_locat>2000&&cand[i].match[0].self_locat>2000)
		{candleft.push_back(cand[i]);continue;}
		if (seed.str.size() - cand[i].match.back().seed_locat < 1000 && cand[i].match.back().seed_locat - cand[i].match[0].seed_locat > 2000 && cand[i].str.size()-cand[i].match.back().self_locat>2000)
		{candright.push_back(cand[i]);continue;}
		if ((cand[i].match.back().self_locat - cand[i].match[0].self_locat) < 0.5 * (int)cand[i].str.size())
		{
			if ((cand[i].match[0].self_locat > 1000)&&(cand[i].str.size() - cand[i].match.back().self_locat > 2000)) 
			cand_2right.push_back(cand[i]);
			if ((cand[i].match[0].self_locat > 1000)&&(cand[i].str.size() - cand[i].match.back().self_locat < 2000)) 
			cand_2left.push_back(cand[i]);
		}
	}
}
void cand_check(inff seed,vector<inff>& cand,int& num_)
{
	int len=cand.size();
	int* K=new int[len];
	for (int i=0;i<len;i++)
		K[i]=i;
	sort(cand.begin(), cand.end(), cmp3);
	for (int i=0;i<len;i++)
	{
		if (cand[i].match_start < 10||(seed.str.size()-cand[i].match_end<10)) continue;
		if ((cand[i].self_end-cand[i].self_start)<0.8*(int)cand[i].str.size())
		num_++;
	}
}
int seed_select(vector<string> read,int* read_flag,vector<pair<int,int> > read_size_list)
{
	int m = 0,tar=-1;
	for (int i = 0; i < read.size(); i++)
	{
		if (m < read[i].size()&&read_flag[i]==1&&read[i].size()>4000)
		{
			m = read[i].size();
			tar = i;
		}
	}
	if (m<10000) sig_end=0;
	return tar;
}
void rf_change(int* read_flag, int i)
{
	if (read_flag[i] == 1)
	{
		read_flag[i] = 0;
		readtem.push_back(i);
	}
}
void rf_change2(int* read_flag, int i)
{
	if (read_flag[i] == 1)
	{
		read_flag[i] = 0;
		readchange.push_back(i);
		readtem.push_back(i);
	}
}
void rf_back(int* read_flag )
{
	for (const auto& i:readchange)
	if (read_flag[i] == 0)
	{
		read_flag[i] = 1;
	}
}
void rf_back(int* read_flag, int i)
{
	if (read_flag[i] == 0)
	{
		read_flag[i] = 1;
	}
}
void rf_back1(int* read_flag, int i)
{
	if (read_flag[i] == 1)
	{
		read_flag[i] = 2;
	}
}
void rf_back2(int* read_flag, int i)
{
	if (read_flag[i] == 0)
	{
		read_flag[i] = 2;
	}
}
int _back(vector<string> read, int* read_flag, vector<int>& Path, unordered_map<int, err>& error, int Right, int start)
{
	int i = Path.size() - 1, mx = 0, mn = 9999999, an = -1, j = 0;
	int tag=-1;
	int I=i;
	for (; i > 0; i--)
	{
		if (error[Path[i]].repeat != 0)
			break;
	}
	if (i < 0)
		return -1;
	i = Path[i];
		for (int vv=0;vv<error[i].clique.size();vv++)
		{
			if (read_flag[error[i].clique[vv][0]]==0) continue;
			for (const auto& v:error[i].clique[vv])
			{
				if (read_flag[v] == 0) continue;
				if (mx < read[v].size())
				{
					mx = read[v].size();
					an = v;	
					tag=vv;
				}
			}
		}
	if (an == -1) return an;
	if (tag!=-1)
		{
		error_change(read_flag,error[i],tag);
		error[i].repeat--;
		}	
		while (Path[j] != start)
			j++;
		for (int k = j; k <= I; k++)
			Path.push_back(Path[k]);
		return an;
		
}
void out_seed(inff seed)
{
	ofstream ggg("seed.fasta");
	ggg<<">"<<seed.id<<"\n"<<seed.str<<endl;
}
void out_cand(vector<inff> &cand,vector<int> &b)
{
	int x;
	ofstream ffg("lis");
	for (const auto& v:b)
	{
		x=cand[v].id;
		if (x>real_read_size)
			x-=real_read_size;
		ffg<<x<<"\n";
	}
	ofstream ggg("cand.fasta");
	for (const auto& v:b)
	{
		ggg<<">"<<cand[v].id<<"\n"<<cand[v].str<<endl;
	}
		
} 
void out_cand_l(vector<inff> &cand,unsigned int tag)// le re in lb rb
{
	int x;
	ofstream ffg("lis");
	for (const auto& v:cand)
	{
		if (!(v.tag&tag))
			continue;
		x=v.id;
		if (x>real_read_size)
			x-=real_read_size;
		ffg<<x<<"\n";
	}
	ofstream ggg("cand.fasta");
	for (const auto& v:cand)
	{
		if (!(v.tag&tag))
			continue;
		ggg<<">"<<v.id<<"\n"<<v.str<<endl;
	}
		
} 
void out_cand_palin(vector<inff> &cand)
{
	int x;
	ofstream ffg("lis");
	for (const auto& v:cand)
	{
		x=v.id;
		if (x>real_read_size)
			x-=real_read_size;
		ffg<<x<<"\n";
	}
	ofstream ggg("cand.fasta");
	for (const auto& v:cand)
	{
		if (v.is_palindrome)
			ggg<<">"<<v.id<<"_"<<v.sim<<"\n"<<v.str<<endl;
	}
		
} 
void out_cand(vector<inff> &cand)
{
	int x;
	ofstream ffg("lis");
	for (const auto& v:cand)
	{
		x=v.id;
		if (x>real_read_size)
			x-=real_read_size;
		ffg<<x<<"\n";
	}
	ofstream ggg("cand.fasta");
	for (const auto& v:cand)
	{
		ggg<<">"<<v.id<<"_"<<v.sim<<"\n"<<v.str<<endl;
	}
		
} 
void out_table(vector<table> &t)
{
	for (const auto& v:t)
	{
	cout<<v.tab<<"\n";
	for (const auto& vv:v.rd)
		cout<<vv<<" ";
	cout<<"\n";
	int x;
	ofstream ffg("lis");
	for (auto vv:t)
	{	
		for (auto v:vv.rd)
		{
			x=v;
			if (x>real_read_size)
				x-=real_read_size;
			ffg<<x<<"\n";
		}
	}
	}
}
void out_err(unordered_map<int,err> error ,int i)
         {
                 cout<< error[i].repeat<<"\n";
                 for (auto v:error[i].clique)
                         {
                         for (auto vv:v)
                                 cout<<vv<<" ";
                         cout<<"\n";
                         }
         }
void out_err(unordered_map<int,err> error ,int i,int* read_flag)
         {
                 cout<< error[i].repeat<<"\n";
                 for (auto v:error[i].clique)
                         {
                         for (auto vv:v)
                                 cout<<vv<<":"<<read_flag[vv]<<" ";
                         cout<<"\n";
                         }
         }

void out_err(unordered_map<int,err> error ,vector<int> Path)
	{
		ofstream fgh("error");
		for (auto v:Path)
		{
		fgh<<"read:"<<v<<"\n";
		fgh<<"repeat:"<<error[v].repeat<<"\n";
		}
	}
void error_change(int* read_flag,err error,int tag)
{
	for (auto v:error.clique[tag])
		rf_change2(read_flag,v);
}			
void out_p(vector<int> &p)
{
	ofstream f("lis");
	for (auto v:p)
	{
		if (v>real_read_size)
			v=v-real_read_size;		
		f<<v<<"\n";
	}
}
void pp_change(int* ppath,int i,int k)
{
	if (ppath[i]<0)
	{
//		if (i>=real_read_size)
//			ppath[i-real_read_size]=k;
//		else
//			ppath[i+real_read_size]=k;
		ppath[i]=k;
		if (correct_start==1)
		creadtem.push_back(i);
	}
}
int find(vector<pair<int,int>> A,int a)
{
	for (int i=0,size=A.size();i<size;i++)
		if (A[i].first==a)
			return 1;
	return 0;
}
void find0(vector<pair<int,int>> A,int a)
{
	for (int i=0,size=A.size();i<size;i++)
		if (A[i].first==a)
			cout<<i<<endl;
	cout<<"no find"<<endl;
}
int find(vector<int> A,int a)
{
	for (auto v:A)
		if (v==a)
			return 1;
	return 0;
}
void pp_back(int* ppath,int** pdirect)
{
	for (auto i:creadtem)
	{
/*		if (i>=real_read_size)
		{
			ppath[i-real_read_size]=-1;
			pdirect[i-real_read_size][0]=-1;
			pdirect[i-real_read_size][1]=-1;
		}
		else
		{
			ppath[i+real_read_size]=-1;
			pdirect[i+real_read_size][0]=-1;
			pdirect[i+real_read_size][1]=-1;
		}
*/		ppath[i]=-1;
		pdirect[i][1]=-1;
		pdirect[i][0]=-1;
	}
}
void check(vector<inff>& cand,int* ppath)
{
	for(auto v:cand)
	{
		cout<<"id:"<<v.id<<"-seed:"<<ppath[v.id]<<endl;
		cout<<"left_overhang:"<<v.self_start<<"\tright_overhang:"<<v.str.size()-v.self_end<<endl;
		cout<<"left_locat:"<<v.match_start<<"\tright_locat:"<<cand[0].str.size()-v.match_end<<endl;
		cout<<"lextend:"<<v.is_lextend<<"\tinsert:"<<v.is_insert<<"\trextend:"<<v.is_insert<<endl;
		cout<<"\nsize:"<<v.str.size()<<endl;
	}
}
void relist(vector<int>& Path,pair<int,int>& _x)
{
	vector<int> t0;
	vector<int> t1;
	vector<int>::iterator t,tt;
	t0=Path;
	Path.clear();
	for (t=t0.begin();t!=t0.end();t++)
		if (*t==_x.first)
			break;
	for (tt=t+1;tt!=t0.end();tt++)
		Path.push_back(*tt);
	Path.push_back(mxseed);
	for (tt=t0.begin();tt<=t;tt++)
		Path.push_back(*tt);
	_x.second--;
}
void check(vector<inff> &cand)
{
	class cc{
	public:
	vector<int> id;
	int count=0;
	vector<int> locat;
	};
	unordered_map<long long int,cc> tt; 
	for (auto v:cand)
		for(auto vv:v.match)
			{
			tt[vv.kmer].id.push_back(v.id);
			tt[vv.kmer].count++;
			tt[vv.kmer].locat.push_back(vv.self_locat);
			}
	ofstream f("kmer_match");
	for (auto v:tt)
		{
		if (v.second.count>1)
		{	f<<n2c(v.first)<<"\t"<<v.second.count;
			f<<"\n";
			for (auto vv:v.second.id)
			f<<"\t"<<vv;
			f<<"\n";
			for (auto vv:v.second.locat)
			f<<"\t"<<vv;
			f<<"\n";
			
		}
		}
}
void soil_km_count(inff seed,vector<inff> cand)
{
	soil x;
	int ii;
	pair<int,int> y;
	vector<soil> list;
	vector<pair<int,int> > yy;
	unordered_map<unsigned long long ,int>countk;
	int *nu=new int[cand.size()];
	for (ii=0;ii<cand.size();ii++)
		nu[ii]=0;
	int minpo=9999999,maxpo=0;
	for (const auto& v:cand)
	{
		if (minpo>v.match[0].seed_locat)
			minpo=v.match[0].seed_locat;
		if (maxpo<v.match.back().seed_locat)
			maxpo=v.match.back().seed_locat;
	}
	for (const auto& v:seed.mini)
	{
		if (v.second<minpo||v.second>maxpo)
			continue;
		x.kmer=v.first;
		x.count=0;
		x.lo=v.second;
		countk[v.first]++;
		if (!yy.empty())
			yy.clear();
		for (int i=0;i<cand.size();i++)
		{
			for (const auto& vv:cand[i].match)
			{
				if (v.first==vv.kmer)
				{
					x.count++;
					y.first=i;
					y.second=vv.self_locat;
					yy.push_back(y);
				}
			}
		}
		x.locat=yy;
		list.push_back(x);
	}
	for (const auto& vv:countk)
	{
		if (vv.second>1)
			cout<<"\nkmer:"<<n2c(vv.first)<<"\tcount:"<<vv.second;
	}	
	sort(list.begin(),list.end(),cmp_list);	
	for(ii=0;ii<list.size();ii++)
	{
		for (const auto& v:list[ii].locat)
		{
			nu[v.first]++;
		}
		if (list[ii].count<4)
			break;
	}
	delete[] nu;	
}
bool cmp_list(soil x, soil y)
 {
         return x.count > y.count;
 }
void outsoil(soil x)
{
	cout<<"\nkmer:"<<n2c(x.kmer)<<"locat:"<<x.lo;
	for (auto v:x.locat)
	{
		cout<<"\ncandid:"<<v.first<<"\tlocat:"<<v.second;
	}
	cout<<endl;
}	
void outaln(vector<pair<char,char> > aln)
{
	int j=aln.size(),i=0;
	for (int ii=0;ii<=aln.size()/80;ii++)
	{
		for (i=0+ii*80;i<aln.size()&&i<(ii+1)*80;i++)
			cout<<aln[i].first; 		
		cout<<"\n";
		for (i=0+ii*80;i<aln.size()&&i<(ii+1)*80;i++)
			cout<<aln[i].second; 	
		cout<<"\n-----------------\n"<<endl;
	}
}	
void outMM(vector<vector<char> > M,vector<Snp> &_snp)
{
	ofstream fgg("111111");
	int j=M.size(),k=M[0].size(),i=0,ii=0,jj=j;
	for (auto v:_snp)
	{
		if (v.locat-10>0)
			ii=v.locat-10;
		else ii=0;
		for (int yy=0;yy<k;yy++)	
		{	
			for (i=ii;i<v.locat+10&&i<j;i++)
				fgg<<M[i][yy]; 		
			fgg<<"\n";
		}
		fgg<<"\n-----------------\n"<<endl;
	}
}	
void outMM(vector<vector<char> > M)
{
	ofstream fgg("111111");
	int j=M.size(),k=M[0].size(),i=0;
	for (int ii=0;ii<=j/80;ii++)
	{
		for (int yy=0;yy<k;yy++)	
		{	for (i=0+ii*80;i<j&&i<(ii+1)*80;i++)
				fgg<<M[i][yy]; 		
			fgg<<"\n";
		}
		fgg<<"\n-----------------\n"<<endl;
	}
}	
void lastcand_update(int* read_flag)
{
	vector<inff> t_last;
	for (const auto& v:lastcand)
	{
		if (read_flag[v.id]==1)
			t_last.push_back(v);
	}
	lastcand.clear();
	lastcand=t_last;
	t_last.clear();
}
bool find(vector<inff> &lastcand,int id)
{
	if (lastcand.empty())
		return false;
	for(const auto& v:lastcand)
	{
		if (id==v.id)
			return true;
	}
	return false;
}
int twin(int x)
{
	if (x>=real_read_size)
		return x-real_read_size;
	else
		//return x+real_read_size;
		return x;
}
bool prefix_check(vector<inff> &cand,int seed_len,int Right,inff &seed,vector<vector<_pair>> &minimizer,vector<string> &read)
{       //only left,prepare to change others
	int cand_len,_start,_end;
	vector<inff> cand0;
	cand0.reserve(cand.size());
	for (int i=0,size=cand.size();i<size;i++)
	{
		if (cand[i].match_start==-1||cand[i].match_end-cand[i].match_start<500) continue;
		cand_len=cand[i].str.size();
		if (cand[i].match_start>exact&&cand[i].self_start>exact&&(seed_len-cand[i].match_end>exact)&&(cand_len-cand[i].self_end>exact))
			continue;
		cand0.push_back(move(cand[i]));
	}
	cand.swap(cand0);
	vector<inff>().swap(cand0);
	int ss0,ss1,ss2,ss3;
	bool Fd=false;
	for (int i=0,size=cand.size();i<size;i++)
	{
		ss0=cand[i].match_start;
		ss1=cand[i].self_start;
		ss2=seed_len-cand[i].match_end;
		ss3=cand[i].str.size()-cand[i].self_end;
		if (ss0<exact&&ss2<exact&&ss1-ss0>200&&ss3-ss2>200)
		{
			cand[i].is_cover=true;
			if (ppath[cand[i].id]!=-1)
			{	cand[i].is_rextend=true;
				cand[i].is_lextend=true;
			}
		}
		if (ss1<exact&&ss3<exact&&ss0>exact&&ss2>exact)
		{
			cand[i].is_insert=true;
		}
		if (ss2<exact)
		{
			coverage[1]++;
			if (cand[i].extend[1]<=0) continue;
			if ((ss1<exact))
			{
				cand[i].is_rextend=true;
				extendsize[1]++;
			}
		}
		if ((ss1<exact||ss0<exact)&&ss2>exact&&ss3>exact)
			cand[i].is_rbranch=true;
		if (ss0<exact)
		{
			coverage[0]++;
			if ((ss3<exact)&&cand[i].extend[0]>0)
			{
				cand[i].is_lextend=true;
				extendsize[0]++;
			}
		}
		if ((ss2<exact||ss3<exact)&&ss1>exact&&ss0>exact)
			cand[i].is_lbranch=true;	
		if (Right)
		{
			if (cand[i].is_lbranch)
			{
		
				Fd=false;
				for (const auto& v:_lastcand)
					if (v.first==cand[i].id)
					{
						Fd=true;
						break;
					}
				if (Fd)
				{
					cand[i].is_lbranch=false;
					cand[i].is_lextend=true;
				}
			}		
		}
		else
		{
			if (cand[i].is_rbranch)
			{
				Fd=false;
				for (const auto& v:_lastcand)
					if (v.first==cand[i].id)
					{
						Fd=true;
						break;
					}
				if (Fd)
				{
					cand[i].is_rbranch=false;
					cand[i].is_rextend=true;
				}
			}
		}
		cand[i].tag_b();
	}
	vector <pair<int,int>> CAND;
	pair<int,int> _x;
	if(Right)
	{
		for(int i=0;i<cand.size();i++)
		{
			if (cand[i].is_rextend)
			{
				_x.first=i;
				_x.second=cand[i].match_end-cand[i].match_start;
				CAND.push_back(_x);
			}
		}
	}
	else
	{
		for(int i=0;i<cand.size();i++)
		{
			if (cand[i].is_lextend&&(cand[i].match_end-cand[i].match_start>min((int)seed.str.size()/2,5000)))
			{
				_x.first=i;
				_x.second=cand[i].match_end-cand[i].match_start;
				CAND.push_back(_x);
			}
		}
	}	
	sort(CAND.begin(),CAND.end(),cmp7);
	int align_num=0;
	int _locat=-1,I;
	vector<int> Locat_;
	if (CAND.size()>0)
	{
		int cand_count=0;
		align_num=0;
		cand_count=0;
		for (const auto& v:CAND)
		{
			cand_count++;
			if (cand[v.first].n_reg>1)
			align_num++;
		}	
		if (align_num<max(5,(int)(cand_count*0.1)))
			end_flag_1=false;
		else
			end_flag_1=true;
	}
	else
		end_flag_1=true;
	if (end_flag_1&&!seed.is_consensus)
	{	
		seed.is_consensus=true;
		consensus_len=seed.str.size();
		seed.cons_len.push_back(consensus_len);
	}
	return end_flag_1;
}
int candclean(vector<inff> &cand,int seed_len,int Right,inff &seed)
{
	int cand_len,_start,_end;
	vector<inff> cand0;
	cand0.reserve(cand.size());
	for (int i=0,size=cand.size();i<size;i++)
	{
		if (cand[i].match_start==-1||cand[i].match_end-cand[i].match_start<500) continue;
		cand_len=cand[i].str.size();
		if (cand[i].match_start>exact&&cand[i].self_start>exact&&(seed_len-cand[i].match_end>exact)&&(cand_len-cand[i].self_end>exact))
			continue;
		cand0.push_back(cand[i]);
	}
	cand.swap(cand0);
	vector<inff>().swap(cand0);
	coverage[0]=0;
	coverage[1]=0;
	extendsize[0]=0;
	extendsize[1]=0;
	int ss0,ss1,ss2,ss3;
	int* cover=new int[seed_len];
	is_repeat(cand,seed_len,cover);
	bool Fd=false;
	for (int i=0,size=cand.size();i<size;i++)
	{
		ss0=cand[i].match_start;
		ss1=cand[i].self_start;
		ss2=seed_len-cand[i].match_end;
		ss3=cand[i].str.size()-cand[i].self_end;
		if (ss0<exact&&ss2<exact&&ss1-ss0>200&&ss3-ss2>200)
		{	//if (ppath[cand[i].id]!=-1&&ppath[cand[i].id]!=ppath[cand[0].id])
			cand[i].is_cover=true;
			if (ppath[cand[i].id]!=-1)
			{	cand[i].is_rextend=true;
				cand[i].is_lextend=true;
			}
		}
		if (ss1<exact&&ss3<exact&&ss0>exact&&ss2>exact)
		{
			cand[i].is_insert=true;
		}
		if (ss2<exact)
		{
			coverage[1]++;
			if (cand[i].extend[1]<=0) continue;
			if ((ss1<exact))
			{
				cand[i].is_rextend=true;
				extendsize[1]++;
			}
		}
		if ((ss1<exact||ss0<exact)&&ss2>exact&&ss3>exact)
			cand[i].is_rbranch=true;
		if (ss0<exact)
		{
			coverage[0]++;
			if ((ss3<exact)&&cand[i].extend[0]>0)
			{
				cand[i].is_lextend=true;
				extendsize[0]++;
			}
		}
		if ((ss2<exact||ss3<exact)&&ss1>exact&&ss0>exact)
			cand[i].is_lbranch=true;	
		if (Right)
		{
			if (cand[i].is_lbranch)
			{
		
				Fd=false;
				for (const auto& v:_lastcand)
					if (v.first==cand[i].id)
					{
						Fd=true;
						break;
					}
				if (Fd)
				{
					cand[i].is_lbranch=false;
					cand[i].is_lextend=true;
					coverage[0]++;
					extendsize[0]++;
				}
			}		
		}
		else
		{
			if (cand[i].is_rbranch)
			{
				Fd=false;
				for (const auto& v:_lastcand)
					if (v.first==cand[i].id)
					{
						Fd=true;
						break;
					}
				if (Fd)
				{
					cand[i].is_rbranch=false;
					cand[i].is_rextend=true;
					coverage[1]++;
					extendsize[1]++;
				}
			}
		}
		cand[i].tag_b();
	}
	cout<<"left:"<<coverage[0]<<"\tright:"<<coverage[1]<<endl;
	cout<<"extendleft:"<<extendsize[0]<<"\textendright:"<<extendsize[1]<<endl;
	Nc=0;
	for (const auto& v:cand)
	{
		if (Right)
		{
			if (!find(_lastcand,v.id)&&v.is_rextend)
				Nc++;
		}
		else
		{
			if (v.is_lextend&&!find(_lastcand,v.id))
				Nc++;
		}
	}
	cout<<"real_extend:"<<Nc<<endl;
	vector<int> lbranch_locat;
	vector<int> rbranch_locat;
	vector<int> lbranch_num;
	vector<int> rbranch_num;
	vector<int> l_locat;
	vector<int> r_locat;
	for (int i=0,size=cand.size();i<size;i++)
	{
		if (cand[i].is_lbranch)
		{
			lbranch_locat.push_back(cand[i].match_start);
			lbranch_num.push_back(i);
		}
		if (cand[i].is_rbranch)
		{
			rbranch_locat.push_back(cand[i].match_end); 
			rbranch_num.push_back(i);
		}
	}
	sort(rbranch_locat.begin(),rbranch_locat.end(),cmp8);
	if (rbranch_locat.size()>10)
	while(1)
	{
		int t=0,i=0;
		bool br=true;
		for (i=0;i<rbranch_locat.size();i++)
		{
			if (rbranch_locat[i]!=-1)
			{
				t=rbranch_locat[i];
				br=false;
				break;
			}
		}
		if (br)
			break;
		rbranch_locat[i]=-1;
		int _count=1;
		for (i=0;i<rbranch_locat.size();i++)
		{
			if (abs(rbranch_locat[i]-t)<10)
			{
				_count++;
				rbranch_locat[i]=-1;
			}
		}
		if (_count>10)
			r_locat.push_back(t);
	}
	sort(lbranch_locat.begin(),lbranch_locat.end(),cmp9);
	if (lbranch_locat.size()>10)
	while(1)
	{
		int t=0,i=0;
		bool br=true;
		for (i=0;i<lbranch_locat.size();i++)
		{
			if (lbranch_locat[i]!=-1)
			{
				t=lbranch_locat[i];
				br=false;
				break;
			}
		}
		if (br) break;
		lbranch_locat[i]=-1;
		int _count=1;
		for (i=0;i<lbranch_locat.size();i++)
		{
			if (abs(lbranch_locat[i]-t)<10)
			{
				_count++;
				lbranch_locat[i]=-1;
			}
		}
		if (_count>10)
			l_locat.push_back(t);
	}
	if (!lbranch_locat.empty()||!rbranch_locat.empty())
	{
		int yyy=0;
		yyy++;
	}
	if (coverage[0]<_deep*rth&&coverage[1]<_deep*rth)
	{
		if (Right)
		{
			for (int i=0,size=cand.size();i<size;i++)
			{
				if (cand[i].is_insert||cand[i].is_lextend)
					//cand[i].is_shutdown=true;
					;
			}
		}
		else 
		{
			for (int i=0,size=cand.size();i<size;i++)
			{
				if (cand[i].is_insert||cand[i].is_rextend)
					//cand[i].is_shutdown=true;
				;
			}
		}
	}
	else if (Right)
	{
		int Deep=coverage[0];
		for (int i=0,size=cand.size();i<size;i++)
		{
			if (!cand[i].is_insert&&!cand[i].is_lextend) continue;
			if (Deep<_deep)
			{
				//cand[i].is_shutdown=true;
				continue;
			}
		//	if (rbranch_locat.empty()) continue;
			if (r_locat.empty()) continue;
		//	for (auto vv:rbranch_locat)
			//for (auto vv:r_locat)
				//if (Deep>cover[vv+1]&&cand[i].match_start<vv&&cand[i].match_end>vv&&cover[vv+1]<_deep)
			//	if (cand[i].match_start<vv&&cand[i].match_end>vv)
	//				Deep=cover[vv+1];
	//		if (Deep<_deep)
				//cand[i].is_shutdown=true;
			//	;
		}
	}
	else
	{
		int Deep=coverage[1];
		for (int i=0,size=cand.size();i<size;i++)
		{
			if (!cand[i].is_insert&&!cand[i].is_rextend) continue;
			if (Deep<_deep)
			{
				//cand[i].is_shutdown=true;
				continue;
			}
			if (rbranch_locat.empty()) continue;
			if (l_locat.empty()) continue;
			for (const auto& vv:l_locat)
			for (const auto& vv:lbranch_locat)
				if (Deep>cover[vv-1]&&cand[i].match_start<vv&&cand[i].match_end>vv&&cover[vv-1]<_deep)
				if (cand[i].match_start<vv&&cand[i].match_end>vv)
					Deep=cover[vv-1];
		//	if (Deep<_deep)
		//		cand[i].is_shutdown=true;
				//;
		}
	}	
	delete[] cover;
	vector <pair<int,int>> CAND;
	pair<int,int> _x;
	if(Right)
	{
		for(int i=0;i<cand.size();i++)
		{
			if (cand[i].is_rextend)
			{
				_x.first=i;
				_x.second=cand[i].match_end-cand[i].match_start;
				CAND.push_back(_x);
			}
		}
	}
	else
	{
		for(int i=0;i<cand.size();i++)
		{
			if (cand[i].is_lextend)
			{
				_x.first=i;
				_x.second=cand[i].match_end-cand[i].match_start;
				CAND.push_back(_x);
			}
		}
	}	
	sort(CAND.begin(),CAND.end(),cmp7);
	int align_num=0;
	int _locat=-1,I;
	vector<int> Locat_;
	if (Right)
		for (const auto& v:l_locat)
			Locat_.push_back(seed_len-v);
	else 
		for (const auto& v:r_locat)
		Locat_.push_back(v);
	if (Locat_.size()>0)
	{
		sort(Locat_.begin(),Locat_.end(),cmp9);
		int Count_=0,j=0,ii;
		for (ii=CAND.size()-1;ii>=10&&j<Locat_.size();ii--)
		{
			int overthr=Locat_[j];	
			if (CAND[ii].second<overthr)
				continue;
			else
			{
				if (ii>10)
				{
					_locat=overthr;
					j++;
				}
				else
					break;
			}		
		}
		seed.prefix=Locat_.back();
	}

	if (_locat==-1)
	{
		if (CAND.size()>100)
			_locat=CAND[100].second-1;
			
		else if (CAND.size()>10)
		{
			if (CAND[10].second>1000)
			_locat=1000;
			else
			_locat=CAND[10].second;
		}
		else
		{
			_locat=100;
		}
	}
	int cut_locat=999999;
	for (int i=0;i<CAND.size();i++)
	{
		if (CAND[i].second>_locat)
		{
			cand[CAND[i].first].if_cand=1;
			cand[CAND[i].first].if_extend=1;
			if (CAND[i].second<cut_locat)
				cut_locat=CAND[i].second;
		}
	}
	for (int i=0,size=cand.size();i<size;i++)
	{				
		int len=cand[i].match_end-cand[i].match_start;
		if (cand[i].match_start<exact&&seed_len-cand[i].match_end<exact&&cand[i].self_start>exact&&cand[i].str.size()-cand[i].self_end>exact);
		//	cand[i].if_extend=1;
			
		if (Right)
		{
			if (cand[i].self_start<exact&&cand[i].str.size()-cand[i].self_end>exact&&seed_len-cand[i].match_end>exact)
			{
				cand[i].if_rbranch=1;
				//cand[i].rbranch_locat=cand[i].c_len-(int)cand[i].str.size()+cand[i].self_end;
				cand[i].rbranch_locat=cand[i].match_end;
			}
			if (len>_locat&&seed_len-cand[i].match_end<exact)
			{
		//		cand[i].if_cand=1;
				if (cand[i].extend[Right]>exact&&cand[i].self_start<exact);
	//				cand[i].if_extend=1;
			}
		}
		else
		{
			if (cand[i].str.size()-cand[i].self_end<exact&&cand[i].self_start>exact&&cand[i].match_start>exact)
			{
				cand[i].if_lbranch=1;
				cand[i].lbranch_locat=cand[i].match_start;
			}
			if (len>_locat&&cand[i].match_start<exact)
			{
		//		cand[i].if_cand=1;
				if (cand[i].str.size()-cand[i].self_end<exact&&cand[i].extend[Right]>exact);
	//				cand[i].if_extend=1;
			
			}
		}
	}
	int ec_count=0,ex_count=0;
	if (Right)	
	{
		int locatec=seed_len-_locat;
		for (const auto& v:cand)
		{
			if (v.is_rextend||v.is_insert||v.is_lextend||v.is_lbranch)
				if (v.match_start<locatec&&v.match_end>locatec)
					ec_count++;
			if (v.if_extend==1)
				ex_count++;
		}
	}
	else
	{
		int locatec=_locat;
		for (const auto& v:cand)
		{
			if (v.is_rextend||v.is_insert||v.is_lextend||v.is_rbranch)
				if (v.match_start<locatec&&v.match_end>locatec)
					ec_count++;
			if (v.if_extend==1)
				ex_count++;
		}

	////////////////////////////////////
//	if (extendsize[Right]<10)
//		return -1;
	if (ec_count!=0)
		return 40*ex_count/ec_count;	
	else 
	{
		if (end_flag_1)
			return 0;
		seed_is_correct=false;
		return -1;
	}///////////////////////////////////
	}
	vector<int>().swap( lbranch_locat);
	vector<int>().swap( rbranch_locat);
	vector<int>().swap( lbranch_num);
	vector<int>().swap( rbranch_num);
	vector<int>().swap( l_locat);
	vector<int>().swap( r_locat);
}
		
	
bool is_repeat(vector<inff> &candtem,int seed_len,int cover[])
{
	int leftmin=9999999,rightmax=0,size=seed_len;
	if (candtem.size()==0)
		return false;	
	for (int i=0;i<size;i++)
		cover[i]=0;
	int maxd=0;
	int start,end,c_start,c_end;
	for (const auto& v:candtem)
	{
		if (v.is_align){
		c_start=v.match_start;
		c_end=v.match_end;
		for (int i=c_start;i<c_end;i++)
			cover[i]++;}
	}
//	for (int i=0;i<2000&&i<size;i++)
//		if (maxd<cover[i])
//			maxd=cover[i];
//	coverage[0]=maxd;
//	maxd=0;
//	for (int i=((size-2000)>0?seed.str.size()-2000:0);i<seed.str.size();i++)
//		if (maxd<cover[i])
//			maxd=cover[i];
//	coverage[1]=maxd;
//	cout<<"left:"<<coverage[0]<<"\tright:"<<coverage[1]<<endl;
//	cout<<"deep:"<<maxd<<endl;
//	delete cover;
	return false;
}
bool seed_check(int seedid,inff &seed,vector<inff> &cand)
{
	_align _align_(-1,-1);	
	int gd=0;
	int a,b;
	for (const auto& v:cand)
	{
		if (v.id==seedid)
			continue;
	//	if (v.is_lextend||v.is_insert||v.is_rextend)
		{
			a=v.match_start;
			b=v.match_end;
	//		_align_.add(v.match_start,v.match_end,100);
			_align_.add(a,b,1000);
			if (gd==1)
			_align_.out();
		}
	}
	if (_align_.size==1&&_align_.next->left<500&&(int)seed.str.size()-_align_.next->right<500)
	{
		_align_.del();
		return true;
	}
	else
	{
		_align_.del();
		return false;
	}		 
}
void cand_select(vector<string> &read,inff &seed,vector<inff> &candtem,vector< vector<pair<int, int > > > & mini_hash,unordered_map<int, vector<rdid > >& read_id,int read_flag[])
{
	rdid p;
	vector<rdid> A;
	unordered_map<int, vector<rdid > > read_temporary;
	inff x;
	vector<_pair> seed_mini;
	seed.mini.clear();
	string ts;
	ts=seed.str;
	string s=rev(ts);
	mizer_scan(s,ts,seed.mini,15,2);
	int i;
	for (const auto& v : seed.mini)
	{
		if (seed.is_centro_consensus||seed.is_consensus||mini_hash[v.first].size()<400)
		for (const auto& vv : mini_hash[v.first])
		{
			p.kmer = v.first;
			p.d = v.second - vv.second;
			p.self_locat = vv.second;
			p.seed_locat = v.second;
			read_temporary[vv.first].push_back(p);
		}
	}
	for (const auto& v : read_temporary)
	{
		if (v.second.size() > threshold)
			read_id[v.first] = v.second;
	}
	read_temporary.clear();
	for (const auto& vvv : read_id)
	{
//		match_len=find_count(vvv.second,A);
//		if (match_len<200||A.size()<match_len/656)
//			continue;
//		if (find_count(vvv.second,A)==0)
//			continue;
		i=vvv.first;
//		if (palign[i]==0) continue;
//		if (read_flag[i]==0) continue;
		if (read[i].size()<1000)
			continue;
		x.in(vvv.first,A,seed.c_len);
		x.str=read[i];
		candtem.push_back(x);
		x.cl();
	}
}
bool seed_check0(inff &seed,vector<inff> &cand,int Right)
{
	int min_overhang=1<<30,r_overhang=1<<30,l_overhang=1<<30;
	for (const auto& v:cand)
	{
			if ((v.self_start<200||v.match_start<200))
			{
				if(r_overhang>seed.str.size()-v.match_end)
				{
					r_overhang=seed.str.size()-v.match_end;
				}
			}
			if ((v.str.size()-v.self_end<200||seed.str.size()-v.match_end<200))
			{
				if ( l_overhang>v.match_start)
				{
					l_overhang=v.match_start;
				}
			}
	}
	min_overhang=max(r_overhang,l_overhang);
	if (min_overhang>200)
		return false;
	return true;	
}
bool seed_check(int seedid,int min_dgr,int min_ove,int seed_len,vector<inff> &cand)
{
	if (cand.size()<=min_dgr)
		return false;
	int size=cand.size();
	int a=0,b=0,j=0;
	vector<int> start;
	vector<int> end;
	vector<pair<int,int>> match_block;
	vector<pair<int,int>> tem_block;
	vector<pair<int,int>> block;
	pair<int,int> t;
	int ss0,ss1,ss2,ss3;
	int count=0;
	for (const auto& v:cand)
	{
                ss0=v.match_start;
                ss1=v.self_start;
                ss2=seed_len-v.match_end;
                ss3=v.str.size()-v.self_end;
		if (ss0>exact&&ss1>exact) continue;
		if (ss2>exact&&ss3>exact) continue;
		t.first=v.match_start;
		t.second=v.match_end;
		match_block.push_back(t);
	}

/*	sort(match_block.begin(),match_block.end(),cmp10);
	start.push_back(match_block[0].first);
	end.push_back(match_block[0].second);
	j=0;
	for (int i=1,size=match_block.size();i<size;i++)
	{
		if (match_block[i].second<match_block[j].second)
			continue;
		start.push_back(match_block[i].first);
		end.push_back(match_block[i].second);
		j=i;
	}*/
	for (int i=0,size=match_block.size();i<size;i++)
	{
		start.push_back(match_block[i].first);
		end.push_back(match_block[i].second);
	}
	sort(start.begin(),start.end());
	sort(end.begin(),end.end());
	j=0;
	for(int i=min_dgr-1;i<start.size();i++)
	{
		a=start[i];
		b=end[j++];
		if (a>b)
			continue;
		t.first=a;
		t.second=b;
		tem_block.push_back(t);
	}
	if (tem_block.size()>0)
	{
	t=tem_block[0];	
	for (int k=1,size=tem_block.size();k<size;k++)
	{
		
		if (tem_block[k].first==t.first||tem_block[k].first<200)
		{
			t.second=t.second>tem_block[k].second?t.second:tem_block[k].second;
			continue;
		}
		else if (tem_block[k].first>t.first&&tem_block[k].second<=t.second)
			continue;

		if (tem_block[k].first<t.second-min_ove||seed_len-t.second<200)
			t.second=tem_block[k].second;
		else
		{
			block.push_back(t);
			t=tem_block[k];
		}
	}
	block.push_back(t);
	}
	if (0)
	{
		if (block.size()==1)
			return true;
		else 
			return false;
	}
	else
	{
		if (block.size()==1&&block[0].first<200&&block[0].second>seed_len-200)
			return true;
		else 
			return false;
	}
}
bool seed_check(vector<Path_> &pPath,vector<inff> &candtem,inff &seed,int Right)
{
	if (pPath.size()==0)
		return false;
	int len=seed.str.size(),Len=0;
	unordered_map<int,bool> Check_;
	int size=pPath.size();
	int cand_count=0,find_count=0;
	for (int i=size-1;i>=0;--i)
	{
		for (auto v:pPath[i].cand)
			Check_[v]=1;
		Len+=pPath[i].len;
		if (Len>len)
			break;
	}
	if (Right)
	{
		for (const auto& v:candtem)
		{
			if (v.is_rextend)
			{
				++cand_count;
				if (v.is_rev&&Check_.find(v.id)!=Check_.end())
					++find_count;
			}
		}
	}
	else
	{
		for (const auto& v:candtem)
		{
			if (v.is_lextend)
			{
				++cand_count;
				if (v.is_rev&&Check_.find(v.id)!=Check_.end())
					++find_count;
			}
		}
	}
	if ((double)find_count/cand_count>3)
		return true;
	else 
		return false;
}
void LCS2(vector<_pair> &s1,vector<_pair> &s2,int **c,int **b,int x,int r)
{
	int len1=s1.size(),len2=s2.size();
	for (int i=0;i<=len1;i++)
		for (int j=0;j<=len2;j++)
		{
			c[i][j]=0;
			b[i][j]=0;
		}
	for (int i=1;i<=len1;++i)
	{
		for (int j=1;j<=len2;++j)
		{
			if (s1[i-1].first==s2[j-1].first&&(s1[i-1].second-s2[j-1].second>x-r&&s1[i-1].second-s2[j-1].second<x+r))
			{
				c[i][j]=c[i-1][j-1]+1;
				b[i][j]=1;
			}
			else 
			{
				if (c[i][j-1]>=c[i-1][j])
				{
					c[i][j]=c[i][j-1];
					b[i][j]=2;
				}
				else
				{
					c[i][j]=c[i-1][j];
					b[i][j]=3;
				}
			}
		}
	}
}	
void LCS(vector<_pair> &s1,vector<_pair> &s2,int **c,int **b)
{
	int len1=s1.size(),len2=s2.size();
	for (int i=1;i<=len1;++i)
	{
		for (int j=1;j<=len2;++j)
		{
			if (s1[i-1].first==s2[j-1].first)
			{
				c[i][j]=c[i-1][j-1]+1;
				b[i][j]=1;
			}
			else 
			{
				if (c[i][j-1]>=c[i-1][j])
				{
					c[i][j]=c[i][j-1];
					b[i][j]=2;
				}
				else
				{
					c[i][j]=c[i-1][j];
					b[i][j]=3;
				}
			}
		}
	}
}	
void LCS_print(vector<pair<int,int>> &LCS_mini,int i,int j,int **b,int *num)
{
	if (i==0||j==0)
		return;
	if (b[i][j]==1)
	{
		pair<int,int> x(i,j);
		LCS_mini[*num].first=i-1;
		LCS_mini[*num].second=j-1;
		*num=*num+1;
		LCS_print(LCS_mini,i-1,j-1,b,num);
	}
	else if (b[i][j]==2)
		LCS_print(LCS_mini,i,j-1,b,num);
	else
		LCS_print(LCS_mini,i-1,j,b,num);
}
void rdid_clean(vector<rdid> &A,unordered_map<ull,int> &n_hash)
{
	vector<rdid>::iterator v=A.begin();
	int k=0,l=0;	
	for (;v!=A.end();v++)
	{
		if (v->kmer==-1) continue;
		auto w=v+1;
		int seed_locat=v->seed_locat;
		int self_locat=v->self_locat;
		int sl_len=v->sl_len;
		int Lenv=Len_(v->kmer);
		int n_v=n_hash[v->kmer];
		for (;w!=A.end();w++)
		{
			
			if (w->kmer==-1||w->seed_locat==seed_locat) 
			{
		//		w++;
				continue;
			}
			else if (sl_len>=w->seed_locat&&sl_len<=w->sl_len)
			{
				if ((w->seed_locat-seed_locat)==(w->self_locat-self_locat))
				{
					int over=sl_len-w->seed_locat+1;
					int Lenw=Len_(w->kmer);
					int n_w=n_hash[w->kmer];
					int b=w->sl_len-seed_locat+1;
					v->weight+=w->weight;
					v->com_num++;
					seed_locat=w->seed_locat;
					self_locat=w->self_locat;
					sl_len=w->sl_len;
					v->sl_len=sl_len;
					Lenv=Lenw;
					n_v=n_w;
					//A.erase(w);
					w->kmer=-1;
				}
				else;
				//	w++;
			}
			else if(sl_len<w->seed_locat)
				break;
			else;
			//	w++;
		}
	}
	vector<rdid> B;
	for (int i=0,sizeA=A.size();i<sizeA;i++)
		if (A[i].kmer!=-1)
			B.push_back(A[i]);
	A.swap(B);
}
void rdid_clean(vector<rdid> &A,vector<int> &n_hash)
{
	vector<rdid>::iterator v=A.begin();
	int k=0,l=0;	
	for (;v!=A.end();v++)
	{
		if (v->kmer==-1) continue;
		auto w=v+1;
		int seed_locat=v->seed_locat;
		int self_locat=v->self_locat;
		int sl_len=v->sl_len;
		int Lenv=Len_(v->kmer);
		int n_v=n_hash[v->s_l];
		for (;w!=A.end();w++)
		{
			
			if (w->kmer==-1||w->seed_locat==seed_locat) 
			{
		//		w++;
				continue;
			}
			else if (sl_len>=w->seed_locat&&sl_len<=w->sl_len)
			{
				if ((w->seed_locat-seed_locat)==(w->self_locat-self_locat))
				{
					int over=sl_len-w->seed_locat+1;
					int Lenw=Len_(w->kmer);
					int n_w=n_hash[w->s_l];
					int b=w->sl_len-seed_locat+1;
					v->weight+=1.0/n_hash[w->s_l]*(1+Lenw*0.16)-(1.0/n_v*(1+Lenv*0.16)+1.0/n_w*(1+Lenw*0.16))*0.5*(double)over/b;
					v->com_num++;
					seed_locat=w->seed_locat;
					self_locat=w->self_locat;
					sl_len=w->sl_len;
					Lenv=Lenw;
					n_v=n_w;
					//A.erase(w);
					v->sl_len=sl_len;
					w->kmer=-1;
				}
				else;
				//	w++;
			}
			else if(sl_len<w->seed_locat)
				break;
			else;
			//	w++;
		}
	}
	vector<rdid> B;
	for (int i=0,sizeA=A.size();i<sizeA;i++)
		if (A[i].kmer!=-1)
			B.push_back(A[i]);
	A.swap(B);
}
void _creat(vector<rdid> &A,vector<_pair>& seed_mini,vector<_pair>& cand_mini,int d)
{
	if (!A.empty())
		A.clear();
	for (int i=0,sizes=seed_mini.size();i<sizes;i++)
	{
		for (int j=0,sizec=cand_mini.size();j<sizec;j++)
		{
			if(seed_mini[i].first==cand_mini[j].first&&seed_mini[i].second-cand_mini[j].second>=d-10&&seed_mini[i].second-cand_mini[j].second<d+10)
			{
			rdid x(	seed_mini[i].first,
				seed_mini[i].second,
				cand_mini[j].second,i,j);
			x.sl_len=x.seed_locat+Len_(x.kmer);
			A.push_back(x);
			}
		}
	}

}
void _creat(vector<rdid> &A,vector<_pair>& seed_mini,vector<_pair>& cand_mini)
{
	if (!A.empty())
		A.clear();
	for (int i=0,sizes=seed_mini.size();i<sizes;i++)
	{
		for (int j=0,sizec=cand_mini.size();j<sizec;j++)
		{
			if(seed_mini[i].first==cand_mini[j].first)
			{
			rdid x(	seed_mini[i].first,
				seed_mini[i].second,
				cand_mini[j].second,i,j);
			x.sl_len=x.seed_locat+Len_(x.kmer);
			A.push_back(x);
			}
		}
	}

}
void rdid_creat(vector<rdid> &A,vector<_pair>& seed_mini,vector<_pair>& cand_mini,vector<pair<int,int>> &LCS_mini,int num)
{
	if (!A.empty())
		A.clear();
	for (int i=num-1;i>=0;--i)
	{
		assert(seed_mini[LCS_mini[i].first].first==cand_mini[LCS_mini[i].second].first);
		rdid x(	seed_mini[LCS_mini[i].first].first,
			seed_mini[LCS_mini[i].first].second,
			cand_mini[LCS_mini[i].second].second,0,0);
		A.push_back(x);
	}

}
void rdid_creat_t(vector<rdid> &A,vector<_pair>& seed_mini,vector<_pair>& cand_mini,int d=-1)
{
	if (d==-1)
		_creat(A,seed_mini,cand_mini);
	else
		_creat(A,seed_mini,cand_mini,d);
}
void rdid_creat(vector<rdid> &A,vector<_pair>& seed_mini,vector<_pair>& cand_mini,int d,int len_s,int len_c)
{
	vector<_pair> seed_mini0;
	vector<_pair> cand_mini0;
	if (len_s>=len_c)
	{
		int locat=len_c;
		for (const auto& v:seed_mini)
		{
			if (v.second<locat)
				seed_mini0.push_back(v);
		}
		cand_mini0=cand_mini;
	}
	else
	{
		int locat=len_c-len_s;
		for (const auto& v:cand_mini)
		{
			if (v.second>=locat)
				cand_mini0.push_back(v);
		}
		seed_mini0=seed_mini;
	}
	if (d==-1)
		_creat(A,seed_mini0,cand_mini0);
	else
		_creat(A,seed_mini0,cand_mini0,d);
}

double chain_clean(vector<_pair> &seed_,vector<_pair> &cand_,vector<pair<int,int>> &LCS_mini,vector<rdid> &A,int **c,int **b,int seed_len,int cand_len,int Right)
{
	vector<vector<pair<int,int>>> chain;
	vector<vector<rdid>> chain_;
	vector<pair<int,int>> pos_diffs;
	int pos_diff=0,count,num=0;
	const int size_x=seed_.size()+1;
	const int size_y=cand_.size()+1;
	pair<int,int> tem;
	for (int i=0,size=A.size();i<size;++i)
	{
		count=0;	
		if(A[i].d==-1) continue;
		pos_diff=A[i].d;
		A[i].d=-1;
		count=1;
		int maxr=100;
		for (int j=i+1;j<size;j++)
		{
			if (A[j].d<pos_diff-maxr||A[j].d>pos_diff+maxr)
				continue;
			++count;
			tem.first=pos_diff*(count-1)/count+A[j].d/(count);
			tem.second=maxr;
			A[j].d=-1;
		}
		if (count>2)
		pos_diffs.push_back(tem);
	}
	for (const auto& v:pos_diffs)
	{
		num=0;
		LCS2(seed_,cand_,c,b,v.first,v.second);
		LCS_print(LCS_mini,size_x-1,size_y-1,b,&num);
		rdid_creat(A,seed_,cand_,LCS_mini,num);
		chain_.push_back(A);
	}
	int length=0;
	A.clear();
	double score=0,score0;
	for (int i=0,size=chain_.size();i<size;++i)
	{
		if ((score0=chain_sc(chain_[i],17.0))>score)
		{
			A.swap(chain_[i]);
			score=score0;
		}
	}
	return score;	
}

double chain(vector<rdid>& A,int Right)
{
	vector<_pair> seed_,cand_;
	for (const auto& v:A)
	{
		seed_.push_back(make_pair(v.kmer,v.seed_locat));
		cand_.push_back(make_pair(v.kmer,v.self_locat));
	}
	sort(seed_.begin(),seed_.end(),cmp14);
	seed_.erase(unique(seed_.begin(),seed_.end()),seed_.end());
	sort(cand_.begin(),cand_.end(),cmp14);
	cand_.erase(unique(cand_.begin(),cand_.end()),cand_.end());
	return chain(seed_,cand_,A,(int)seed_.size(),(int)cand_.size(),Right);
	
}
double chain(vector<vector<_pair>> &minimizer,int seed_id,int cand_id, vector<rdid>& A,int seed_len,int cand_len,int Right)
{
	vector<_pair> seed_=minimizer[seed_id];
	vector<_pair> cand_=minimizer[cand_id];
	return chain(seed_,cand_,A,seed_len,cand_len,Right);	
}
double chain(vector<_pair> &seed_,vector<_pair> &cand_, vector<rdid>& A,int seed_len,int cand_len,int Righ)
{
	void LCS(vector<_pair> &s1,vector<_pair> &s2,int** c,int **b);
	void LCS_print(vector<pair<int,int>> &LCS_mini,int i,int j,int **b,int *num);
	const int size_x=seed_.size()+1;
	const int size_y=cand_.size()+1;
	int **c=new int*[size_x],**b=new int*[size_x];
	int num=0;
	for (int i=0;i<size_x;++i)
	{
		c[i]=new int[size_y];
		b[i]=new int[size_y];
	}
	for (int i=0;i<size_x;++i)
		c[i][0]=0;	
	for (int j=0;j<size_y;++j)
		c[0][j]=0;	
	LCS(seed_,cand_,c,b);
	vector<pair<int,int>> LCS_mini;
	LCS_mini.resize(size_x);
	LCS_print(LCS_mini,size_x-1,size_y-1,b,&num);
	rdid_creat(A,seed_,cand_,LCS_mini,num);
	for (int i=0;i<size_x;i++)
	{
		delete[] c[i];
		delete[] b[i];
	}
	delete[] c;
	delete[] b;
	//return score1_(A);
	return A.size();
}
int log_(int x)
{
	int i=1,j=0;
	while (i<x)
	{
		i=i>1;
		++j;
	}
	return j;
}
double chain_sc(vector<rdid> &A,double w)
{
	vector<double> f;
	double f1=w,max_f=0;
	f.push_back(f1);
	double score=f1;
	for (int i=1,size=A.size();i<size;i++)
	{
		max_f=0;
		for (int j=0;j<i;j++)
		{
			double d1=A[i].seed_locat-A[j].seed_locat;
			double d2=A[i].self_locat-A[j].self_locat;
			double d3=abs(d1-d2);
			if (d3<=1.0)
				d3=0;
			//f1=max(f[j]+min(min(d1,d2),w)-(100*d3+0.5*log2(d3+1)),w);
			f1=max(f[j]+min(min(d1,d2),w)-(0.5*log2(d3+1)),w);
			if (f1>max_f)
				max_f=f1;	
		}
		f.push_back(max_f);
		score+=max_f;
	}
	return score;
}
bool cand_count(vector<inff> &cand )
{
	int rextend=0,lextend=0,inextend=0;
	for (const auto& v:cand)
	{
		if (v.is_rextend)
			rextend++;
		if (v.is_insert)
			inextend++;
		if (v.is_lextend )
			lextend++;
	}
	cout<<"left:"<<lextend<<"\tinsert:"<<inextend<<"\tright:"<<rextend<<endl;
	return true;
}
void if_seed_change(vector<inff> &cand)
{
	for (const auto& v:cand)
	{
		if (v.is_rextend||v.is_insert||v.is_lextend)
		{
			if_seed[v.id]=false;
			readtem2.push_back(v.id);
		}
	}
}
void cf_change(vector<int> &readtem2,bool cand_flag[])
{
	for (const auto& v:readtem2)
		cand_flag[v]=true;
}
int cluster_len(vector<rdid> &cluster)
{
	int max=0,min=1<<30-1;
	for (const auto& v:cluster)
	{
		if (max<v.seed_locat)
			max=v.seed_locat;
		if (min>v.seed_locat)
			min=v.seed_locat;
	}
	return max-min;
}
double single_link(vector<rdid> &A,vector<_pair> &seed_)
{
	const int len=A.size();
	vector<vector<rdid>> cluster;
	vector<rdid> ty;
	int **dis_mtx=new int*[len];
	bool flag=false;
	int i=0,j=0,inf=1<<30-1;
	for (i=0;i<len;i++)
		dis_mtx[i]=new int[len];
	
	for (i=0;i<len;i++)
		for (j=0;j<len;j++)
			dis_mtx[i][j]=inf;
	for (i=0;i<len;i++)
		for (j=0;j<len;j++)
		{
			dis_mtx[i][j]=abs(A[j].d-A[i].d);
			dis_mtx[i][j]=abs(A[j].d-A[i].d);
		}
	for (i=0;i<len;i++)
	{
		ty.push_back(A[i]);
		cluster.push_back(ty);
		ty.clear();
	}
	int x,y,min_dis_mtx;
	while(1)
	{
		min_dis_mtx=1<<30-1;
		for (i=0;i<len;i++)
			for (j=i+1;j<len;j++)
			{
				if (min_dis_mtx>dis_mtx[i][j])
				{
					x=i;
					y=j;
					min_dis_mtx=dis_mtx[i][j];
				}
			}
		if (min_dis_mtx>100||min_dis_mtx==1<<30-1)
			break;
		for (const auto& v:cluster[y])
			cluster[x].push_back(v);
		cluster[y].clear();
			for (j=0;j<len;j++)
			{
				dis_mtx[x][j]=min(dis_mtx[x][j],dis_mtx[y][j]);
				dis_mtx[j][x]=min(dis_mtx[j][x],dis_mtx[j][y]);
				dis_mtx[y][j]=inf;
				dis_mtx[j][y]=inf;
			}
	}
	for (auto i=cluster.begin();i!=cluster.end();)
	{
		if (i->size()<2)
			i=cluster.erase(i);
		else if (cluster_len(*i)<100)
			i=cluster.erase(i);
		else
		{
			chain(*i,0);
			i++;
		}
	}
	double score=0,score0=0;
	pair<int,int> chain_check(vector<rdid> &A,vector<vector<rdid>> &cluster,vector<_pair> &seed_);
	pair<int,int> id_;
	for (int i=0;i<cluster.size();i++)
	{
		if (cluster[i].size()!=0)
			if ((score0=cluster[i].size())>score)
			{
	//			AA=cluster[i];
				score=score0;
				id_.first=i;
				id_.second=score;
			}
	}
	//id_=chain_check(A,cluster,seed_);
	for (i=0;i<len;i++)
		delete[] dis_mtx[i];
	delete[] dis_mtx;
	A.swap(cluster[id_.first]);
	return id_.second;	
}
pair<int,int> chain_check(vector<rdid> &A,vector<vector<rdid>> &cluster,vector<_pair> &seed_)
{
	const int len=cluster.size(),len1=seed_.size();
	int *chain_c=new int[len1];
	int *chain_s=new int[len];
	for (int i=0;i<len1;i++)
	{
		chain_c[i]=0;
	}
	for (int i=0;i<len;i++)
	{
		chain_s[i]=0;
		if (cluster[i].size()!=0)
		{
			for (int j=0;j<cluster[i].size();j++)
			{
				for (int k=0;k<len1;k++)
				{
					if (cluster[i][j].seed_locat==seed_[k].second&&cluster[i][j].kmer==seed_[k].first)
					{
						chain_c[k]++;
						break;
					}
				}
			}
		}
	}
	for (int i=0;i<len;i++)
	{
		if (cluster[i].size()!=0)
		{
			for (int j=0;j<cluster[i].size();j++)
			{
				for (int k=0;k<len1;k++)
				{
					if (cluster[i][j].seed_locat==seed_[k].second&&cluster[i][j].kmer==seed_[k].first)
						if (chain_c[k]==1)
							chain_s[i]++;
				}
			}
		}
	}
	pair<int,int> maxs;
	maxs.first=0;
	maxs.second=0;
	for (int i=0;i<len;i++)
	{
		if (maxs.second<chain_s[i])
		{
			maxs.second=chain_s[i];
			maxs.first=i;
		}
	}
	A.clear();
	for (int i=0;i<cluster[maxs.first].size();i++)
		A.push_back(cluster[maxs.first][i]);
	delete[] chain_c;
	delete[] chain_s;
	return maxs;			
}
int Len_(unsigned long long kmer)
{
	return kmer>>(2*k_lenmax);           //1209
}
double count_c(vector<rdid> &A)
{
	vector<double> p;
	double q;
	int ds,dc,d,size=A.size();
	for (int i=0;i<size;i++)
		for (int j=i+1;j<size;j++)
		{
			ds=A[j].seed_locat-A[i].seed_locat;
			dc=A[j].self_locat-A[i].self_locat;
			if (ds>0&&dc>0)
			{
				q=(double)abs(A[j].d-A[i].d)/min(ds,dc);
				p.push_back(q);
			}
		}
	sort(p.begin(),p.end());
	size=p.size();	
	if (size>=4)
	return p[size*0.5]+p[size*0.75]-p[size*0.25];
	else
	return 0.2;
}
void out_match(vector<inff> &cand_,unordered_map<ull,vector<_pair>> &hash,int i)
{
	inff cand=cand_[i];	
	ofstream hj("match_info");
	for (const auto& v:cand.match)
	{
		ull t=v.kmer;
		hj<<t<<"\t"<<n2c(t,Len_(t))<<"\n";
		hj<<"size_0:"<<hash[t].size()<<"\tsize_1:"<<size_count(hash[t])<<endl;
		hj<<"seed_locat:"<<v.seed_locat<<"\tself_locat:"<<v.self_locat<<"\ts_len:"<<v.sl_len<<"\td:"<<v.d<<endl;
	}
}
void out(vector<inff> &cand,unordered_map<ull,vector<_pair>> &hash,ull t)
{
	ofstream hj("hash_c");
	hj<<"size_0:"<<hash[t].size()<<"\tsize_1:"<<size_count(hash[t])<<endl;
	hj<<n2c(t,Len_(t))<<endl;
	if (hash[t].size()!=0)
		for (const auto& v:hash[t])
		{
			hj<<cand[v.second].id<<"\tlocat:\t"<<v.first<<endl;
		}
}
void erase(vector<inff> &cand,int i)
{
	auto v=cand.begin();
	v+=i;
	cand.erase(v);
}
void find0_cand(vector<inff> &cand,int id)
{
	int i=0;
	for (;i<cand.size();i++)
		if (cand[i].id==id)
		{
			cout<<i<<"\t"<<cand[i].score_<<endl;
			break;
		}
	if (i==cand.size())
		cout<<"no find"<<endl;
}
void inff::flag_build(int seed_len)
{
	int ss0=match_start;
	int ss1=self_start;
	int ss2=seed_len-match_end;
	int ss3=str.size()-self_end;
	if (ss0<exact&&ss2<exact&&ss1-ss0>200&&ss3-ss2>200)
	{	
		is_cover=true;
		is_rextend=true;
		is_lextend=true;
	}
	if (ss1<exact&&ss3<exact&&ss0>exact&&ss2>exact)
		is_insert=true;
	if (ss2<exact&&ss1<exact&&extend[1]>0)
		is_rextend=true;
	if (ss0<exact&&ss3<exact&&extend[0]>0)
		is_lextend=true;
	if ((ss1<exact||ss0<exact)&&ss2>exact&&ss3>exact)
		is_rbranch=true;
	if ((ss2<exact||ss3<exact)&&ss1>exact&&ss0>exact)
		is_lbranch=true;	
	tag_b();
	
}
void cand_erase(vector<string> &cand,int id)
{
	vector<string> t;
	for (int i=0;i<cand.size();i++)
		if (i!=id)
			t.push_back(move(cand[i]));
	cand.swap(t);
}
void cand_erase(vector<inff> &cand,int id)
{
	vector<inff> t;
	for (int i=0;i<cand.size();i++)
		if (i!=id)
			t.push_back(move(cand[i]));
	cand.swap(t);
}
_pair rare_kmer_serch(_pair &t,unordered_map<ull,vector<_pair>> &hash,int &filter,unordered_map<ull,int> &real_count,unordered_map<ull,int> &rare_mini_count,vector<int> &mini_n_rare,int up_bd,int down_bd)
{
	_pair ti=t;
	int n;
	rare_serch_r(ti,hash,rare_mini_count,n);
	ull V=(1LL<<(2*k_lenmax))-1;
	ull kmer=t.first&V,t_kmer0;
	ull Kmer;
	rare_serch_l(ti,hash,rare_mini_count,n);
	if (n<1)
		n=1;
	real_count[ti.first]=n;
	mini_n_rare.push_back(max(1,n/rare_mini_count[ti.first]));
	return ti;
}
void inff::mini_update(unordered_map<ull,vector<_pair>> &hash,int down_bd,int up_bd,int len,bool comp)
{
	vector<pair<ull,int>> tem;
	int locat_t=-1;
	int filter=0;
	vector<int> t_n_rare;
	for (auto v:mini2)
	{
		if (v.second<=locat_t) 
			continue;
		if ((v.first>>(2*k_lenmax))!=k_lenmax)
			break;
		{
			_pair ty=rare_kmer_serch(v,hash,filter,real_count,rare_mini_count,t_n_rare,up_bd,down_bd);
			if (ty.second!=-1)
			{
				tem.push_back(ty);
				locat_t=ty.second;
			}
		}
	}
	if (!comp)
	{
		mini2=tem;
		mini.swap(tem);
		mini_n_rare.swap(t_n_rare);
		for (int i=0,size_tem=tem.size();i<size_tem;i++)
		{
			tem[i].second+=len;
			mini.push_back(tem[i]);
			mini_n_rare.push_back(t_n_rare[i]);
		}
	}
	else
	{
		int k=str.size()-cons_len.back();
		for (int i=0,size_tem=tem.size();i<size_tem;i++)
		{
			tem[i].second+=k;
			mini.push_back(tem[i]);
			mini_n_rare.push_back(t_n_rare[i]);
		}
		mini2=tem;
	}
}
/*
int size_count(vector<_pair> &v)
{
	if (v.size()==0)
		return 0;
	int size=1,id=v[0].second;
	for (auto vv:v)
		if (vv.second!=id)
		{
			size++;
			id=vv.second;
		}
	return size;
}*/

inline int size_count(vector<_pair> &v)
{
	return v.size();
}
void out_rdid(vector<vector<rdid>> &T)
{
	for (auto v:T)
	{
		cout<<"------\n";
		out_rdid(v);
	}
}
void out_rdid(vector<rdid> &T)
{
	for (auto v:T)
		cout<<"seed:"<<v.seed_locat<<"self:"<<v.self_locat<<"d:"<<v.d<<"sl_len:"<<v.sl_len<<endl;
}
void out_id0(vector<pair<int,double>> &id0,vector<inff> &cand)
{
	for (auto v:id0)
	{
		cout<<cand[v.first].id<<endl;
	}
}
void update_cons_len(vector<int> &T,int k)
{
	while (k>0)
	{
		if (T.size()==1)
			break;
		for (auto v=T.end()-1;v!=T.begin();v--)
		{
			if (k-*v>=0)
			{
				k=k-*v; 
				T.erase(v);
				break;
			}
			else
			{
				*v=*v-k;
				k=0;
				break;
			}
		}
	}
}
void consensus_totallink(string &q,int Right)
{
	if (consensus_total.size()==0)
	{
		consensus_total=q;
		return;
	}
	if (Right==0)
	{
		int len=q.size();
		string ts=consensus_total.substr(0,len);
		string tc=q.substr(q.size()-500,500);
		EdlibAlignResult result = edlibAlign(tc.c_str(),tc.size(),ts.c_str(),ts.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0 ));
		char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
		string ciga=cigar;
		cigar2string_link(ciga,tc,ts,*result.startLocations);
		consensus_total=q+consensus_total.substr(*result.startLocations,consensus_total.size());
		free(cigar);
	}
}
void cigar2string_link(string &cigar,string &t,string &c,int start)
{
	int num=0,t_l=0,c_l=0,total_len=0,insert,delet,mismatch,match;
	insert=delet=mismatch=match=0;
	string c_align,t_align;
	for (int i=0,size=cigar.size();i<size;i++)
	{
		if (cigar[i]>='0'&&cigar[i]<='9')	
		{
			num=num*10+cigar[i]-'0';
			continue;
		}
		else if (cigar[i]=='M')
		{
			c=t.substr(0,start)+c;
			break;
		}
		else if (cigar[i]=='I')
		{
			c=c.substr(num,c.size());
			c=t.substr(0,start)+c;
			break;
		}
		else if (cigar[i]=='D')
		{
			c=t.substr(0,start+num)+c;
			break;
		}
	}
}
void size_(unordered_map<ull,vector<_pair>>& hash,ull i)
{
	cout<<hash.size()<<"\t"<<hash[i].size()<<endl;
}
bool check(vector<pair<int,int>> &_lastcand,vector<inff> &cand)
{
	if (_lastcand.size()==0)
		return true;
	unordered_map<int,int> T;
	T.reserve(cand.size());
	for (const auto& v:cand)
	{
		if (v.is_insert||v.is_lextend||v.is_rextend)
			T[v.id]=1;
	}
	double num=0,total=0;
	for (const auto& v:_lastcand)
	{
		if (T.find(v.first)!=T.end())
			num++;
		if (v.first!=-1)
			total++;
	}
	if ((double)num/total<0.3)
		return false;
	else
		return true;		
}
void out(vector<string> &T)
{
	ofstream kl("seed.fasta");
	int i=0;
	for (;i<T.size();i++)
		kl<<">read_"<<i<<"\n"<<T[i]<<"\n";
}
void out(vector<string> &T,int i)
{
	ofstream kl("checking");
	for (int ii=0;ii<T.size();ii++)
	kl<<T[ii][i+1];
	kl<<endl;
}
void out_comp(vector<int> &T)
{
	char h[256];
	h['A']='A';
	h['T']='T';
	h['C']='C';
	h['G']='G';
	ofstream kl("checking");
	for (int ii=0;ii<T.size();ii++)
		kl<<h[T[ii]&255];
	kl<<endl;
}	
  double n_hash_score(unordered_map<ull,int> &n_hash,int i,ull t)
{
	double q=1.0/(n_hash[t])*(1+0.16*Len_(t));
}
double n_hash_score(vector<int> &n_hash,int i,ull t)
{
	double q=1.0/n_hash[i]*(1+0.16*Len_(t));
}
void out(unordered_map<ull,int> &T,ull i)
{
	cout<<T[i]<<endl;
}
void seed_rchose(inff &seed,int i)
{
	ifstream kj("seed0.fasta");
	vector<string> t;
	string s;
	while(!kj.eof())
	{
		getline(kj,s);
		getline(kj,s);	
		t.push_back(s);	
	}
	seed.str=t[i];
}
void out(vector<_pair> &T,unordered_map<ull,int> &c,unordered_map<ull,vector<_pair>> &hash)
{
	ofstream k("rare_kmer_num");
	for (int i=0;i<T.size();i++)
	{
		k<<n2c(T[i].first,Len_(T[i].first))<<"\t"<<T[i].first<<"\t"<<c[T[i].first]<<"\t"<<hash[T[i].first].size()<<"\t"<<T[i].second<<"\n";

	}
}
double varance(double *T,int left,int right)
{
	if (right==left) return 0;
	double x=0,y=0;
	for (int i=left;i<=right;i++)
		y+=T[i];
	y/=right-left+1;
	for (int i=left;i<=right;i++)
		x+=pow(y-T[i],2);
	x/=right-left+1;	
	return x;
}
int peek_serch(double *T,int w_size)
{
	int start,end;
	for (int i=k_len;i<k_lenmax+3;i++)
		if (T[i]!=-1)
		{
			start=i;break;
		}
	for (int i=k_lenmax+2;i>k_len-1;i--)
		if (T[i]!=-1)
		{
			end=i;break;
		}
	if (end<start+3) 
		return start; 
	vector<double> t;
	int k=-1;
	for (int i=start;i<=end-2;i++)
		t.push_back(varance(T,i,i+2));
	for (int i=t.size()-1;i>=0;i--)
		if (t[i]<0.33)
			k=i;
		else
			break;
	return k+start;
}
int rare_serch_r(_pair &t,unordered_map<ull,vector<_pair>> &hash,unordered_map<ull,int> &rare_mini_count,int &n)
{
	ull i=0;
	ull V=(1LL<<(2*k_lenmax))-1;
	ull kmer=t.first&V;
	ull Kmer;
	ull vb;
	double *nr=new double[k_lenmax+3];
	double x_size;
	double x[25]={-1};
	int T[25]={-1};
	for (int k=0;k<k_lenmax+3;k++)
		nr[k]=-1;
	int len=Len_(t.first);
	bool flag=false;
	for (i=len;i>=k_len;i--)
	{
		ull t_kmer=(kmer>>((len-i)*2))|(i<<(2*k_lenmax));        //1207
		vb=t_kmer;
		if (hash[t_kmer].size()>3)
		{
			flag=true;
			double y=hash[t_kmer].size();
			T[i]=y;
			double y1=y-1/6*(CAND_SIZE/HOR_LEN-y/P_[i][0])*P_[i][1];
//			double y1=y-check(t_kmer,hash);
			x_size=y1/(P_[i][0]*40);
			x[i]=x_size;
//			if (x_size/rare_mini_count[t_kmer]>=1)
				nr[i]=x_size;
				//nr[i]=x_size/rare_mini_count[t_kmer];
//			else
//				nr[i]=1;
		}
	}
	for (int i=k_lenmax+2;i>=10;--i)
		if (nr[i]!=-1)
		{
			len=i;
			break;
		}
	for (int i=k_len;i<=len;i++)
		nr[i]/=nr[len];
	nr[len+2]=nr[len+1]=nr[len];
	ull k=len;
/*	for (int i=len+2;i>0;i--)
	{
		if (nr[i-1]<nr[i])
			nr[i-1]=nr[i];
	}
*/	if (flag)
		k=peek_serch(nr,3);
	t.first=(kmer>>((k_lenmax-k)*2))|(k<<(k_lenmax*2));        //1207
	int y=hash[t.first].size();
	n=x[k]+1;
	delete[] nr;
}
int rare_serch_l(_pair &t,unordered_map<ull,vector<_pair>> &hash,unordered_map<ull,int> &rare_mini_count,int &n)
{
	ull i=0;
	ull V=(1LL<<(2*k_lenmax))-1;
	ull kmer=t.first&V;
	ull Kmer;
	ull vb;
	double *nr=new double[k_lenmax+3];
	double x_size;
	double x[25]={-1};
	int T[25]={-1};
	for (int k=0;k<k_lenmax+3;k++)
		nr[k]=-1;
	int len=Len_(t.first);
	int len0=len;
	bool flag=false;
	for (i=len;i>=k_len;i--)
	{
		ull t_kmer=(kmer&((1LL<<(2*i))-1))|(i<<(k_lenmax*2));
		vb=t_kmer;
		if (hash[t_kmer].size()>3)
		{
			flag=true;
			double y=hash[t_kmer].size();
			T[i]=y;
			double y1=y-1/6*(CAND_SIZE/HOR_LEN-y/P_[i][0])*P_[i][1];
			//double y1=y-check(t_kmer,hash);
			x_size=y1/(P_[i][0]*40);
			x[i]=x_size;
//			if (x_size/rare_mini_count[t_kmer]>=1)
				nr[i]=x_size;
				//nr[i]=x_size/rare_mini_count[t_kmer];
//			else
//				nr[i]=1;
		}
	}
	for (int i=k_lenmax+2;i>=10;--i)
		if (nr[i]!=-1)
		{
			len=i;
			break;
		}
	for (int i=k_len;i<=len;i++)
		nr[i]/=nr[len];
	nr[len+2]=nr[len+1]=nr[len];
/*	for (int i=len+2;i>0;i--)
	{
		if (nr[i-1]<nr[i])
			nr[i-1]=nr[i];
	}
*/	ull k=len;
	if(flag)
		k=peek_serch(nr,3);
	t.first=(kmer&((1LL<<(2*k))-1))|(k<<(2*k_lenmax));
	int y=hash[t.first].size();
	n=x[k]+1;
	t.second+=len0-k;
	delete[] nr;
}
