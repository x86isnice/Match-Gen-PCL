#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <ctime>
#include <windows.h> 
#include <sstream>
#include <float.h>
#include <pcl/point_types.h>
//#include <pcl/io/pcd_io.h>

#include <Eigen/Dense>  
#include <vector>
#include <math.h>

using namespace Eigen;
using namespace Eigen::internal;  
using namespace Eigen::Architecture;  
using namespace pcl;
using namespace std;

const int ERROR_RANGE = 5;
const int COUNT_RANGE = 10;
int TTTTT = 0;
int KKKKK = 0;


inline  int string2int(std::string s)
{
    int n = atoi(s.c_str());
	return n;
}

inline float string2float(std::string str)
{
	float result = atof(str.c_str());
	return result;
}

inline double string2double(std::string str)
{
	double result;
	stringstream sstr(str);
	sstr >> result;
	return result;
}

inline long string2long(std::string str)
{
   long result = atol(str.c_str());
   return result;
}

std::string double2String(double d) 
{  
	std::ostringstream os;  
	 if (os << d)
	    return os.str();
	 else
		 return "double2string error!";
}  

inline void cutstring(std::string s ,int a[])
{ 
	std::string flag = ":";  
	std::string::size_type position=0; 
	int i=0;  
	while((position=s.find_first_of(flag,position)) != std::string::npos)  
	{  
		//	  std::cout<< "position  " << i << " : " << position << std::endl;  
		a[i] = position;
		position++;  
		i++;  
	}  
}

//Akima方法
class  akima
{
private: 
	int n, k;
	double  *x, *y, t, z, s[4];
public:
	akima (int nn)
	{
		n = nn;
		x = new double[n];    //动态分配内存
		y = new double[n];
	}
	void input(const std::vector<double> &X, const std::vector<double> &Y);          //由文件读入n个数据点(x, y)
	double interp(double);  //计算插值点t所在子区间上的三次多项式
	//并计算t点的近似值z
	void output(); //输出三次多项式系数以及t点的近似值z到文件并显示
	~akima()
	{
		delete []x;  
	    delete []y;
	}
};

void akima::input(const std::vector<double> &X, const std::vector<double> &Y)
{
    for (int i = 0; i < X.size(); i++)
	{
	 x[i] = X[i];
	 y[i] = Y[i];
	}
}

double akima::interp(double tt)    ////计算插值点t所在子区间上的三次多项式
{ 
	//并计算t点的近似值z
	int m,kk;
	double u[5],p,q;
	t = tt;
	z=0.0; s[0]=0.0; s[1]=0.0; s[2]=0.0; s[3]=0.0;
	if (n<1) { k = 0; return z; }
	if (n==1) { k = 0; s[0]=y[0]; z=y[0]; return z;}
	if (n==2)
	{ 
		k = 0;
		s[0]=y[0]; s[1]=(y[1]-y[0])/(x[1]-x[0]);
		z=(y[0]*(t-x[1])-y[1]*(t-x[0]))/(x[0]-x[1]);
		return z;
	}
	if (t<=x[1]) k=0;
	else if (t>=x[n-1]) k=n-2;
	else
	{ 
		k=1; m=n;
		while (((k-m)!=1)&&((k-m)!=-1))
		{ 
			kk=(k+m)/2;
			if (t<x[kk-1]) m=kk;
			else k=kk;
		}
		k=k-1;
	}
	u[2]=(y[k+1]-y[k])/(x[k+1]-x[k]);
	if (n==3)
	{ 
		if (k==0)
		{ 
			u[3]=(y[2]-y[1])/(x[2]-x[1]);
			u[4]=2.0*u[3]-u[2];
			u[1]=2.0*u[2]-u[3];
			u[0]=2.0*u[1]-u[2];
		}
		else
		{ 
			u[1]=(y[1]-y[0])/(x[1]-x[0]);
			u[0]=2.0*u[1]-u[2];
			u[3]=2.0*u[2]-u[1];
			u[4]=2.0*u[3]-u[2];
		}
	}
	else
	{ 
		if (k<=1)
		{ 
			u[3]=(y[k+2]-y[k+1])/(x[k+2]-x[k+1]);
			if (k==1)
			{ 
				u[1]=(y[1]-y[0])/(x[1]-x[0]);
				u[0]=2.0*u[1]-u[2];
				if (n==4) u[4]=2.0*u[3]-u[2];
				else u[4]=(y[4]-y[3])/(x[4]-x[3]);
			}
			else
			{ 
				u[1]=2.0*u[2]-u[3];
				u[0]=2.0*u[1]-u[2];
				u[4]=(y[3]-y[2])/(x[3]-x[2]);
			}
		}
		else if (k>=(n-3))
		{ 
			u[1]=(y[k]-y[k-1])/(x[k]-x[k-1]);
			if (k==(n-3))
			{ 
				u[3]=(y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
				u[4]=2.0*u[3]-u[2];
				if (n==4) u[0]=2.0*u[1]-u[2];
				else u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
			}
			else
			{ 
				u[3]=2.0*u[2]-u[1];
				u[4]=2.0*u[3]-u[2];
				u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
			}
		}
		else
		{ 
			u[1]=(y[k]-y[k-1])/(x[k]-x[k-1]);
			u[0]=(y[k-1]-y[k-2])/(x[k-1]-x[k-2]);
			u[3]=(y[k+2]-y[k+1])/(x[k+2]-x[k+1]);
			u[4]=(y[k+3]-y[k+2])/(x[k+3]-x[k+2]);
		}
	}
	s[0]=fabs(u[3]-u[2]);
	s[1]=fabs(u[0]-u[1]);
	if ((s[0]+1.0==1.0)&&(s[1]+1.0==1.0))
		p=(u[1]+u[2])/2.0;
	else p=(s[0]*u[1]+s[1]*u[2])/(s[0]+s[1]);
	s[0]=fabs(u[3]-u[4]);
	s[1]=fabs(u[2]-u[1]);
	if ((s[0]+1.0==1.0)&&(s[1]+1.0==1.0))
		q=(u[2]+u[3])/2.0;
	else q=(s[0]*u[2]+s[1]*u[3])/(s[0]+s[1]);
	s[0]=y[k];
	s[1]=p;
	s[3]=x[k+1]-x[k];
	s[2]=(3.0*u[2]-2.0*p-q)/s[3];
	s[3]=(q+p-2.0*u[2])/(s[3]*s[3]);
	p=t-x[k];
	z=s[0]+s[1]*p+s[2]*p*p+s[3]*p*p*p;
	return z;
}

void akima::output ()//输出三次多项式系数以及t点的近似值z到文件并显示
{
	char str2[20];
	cout <<"\n输出文件名:  ";
	cin >>str2;
	ofstream fout (str2, ios::app);
	if (!fout)
	{ cout <<"\n不能打开这个文件 " <<str2 <<endl; exit(1); }
	fout <<endl;  cout <<endl;
	fout <<k <<":" <<endl;
	fout <<s[0] <<"   " <<s[1] <<"   " <<s[2] <<"   " <<s[3] <<endl;
	cout <<k <<":" <<endl;
	cout <<s[0] <<"   " <<s[1] <<"   " <<s[2] <<"   " <<s[3] <<endl;
	fout <<endl <<t <<"   " <<z <<endl;
	cout <<endl <<t <<"   " <<z <<endl;
	fout.close ();
}
/*
void main ()      //主函数
{
	akima  solution(11); 
	solution.input ();          //由文件读入n个数据点(x, y)
	solution.interp (-0.85);         //执行Akima方法
	solution.output ();//输出三次多项式系数以及t点的近似值z到文件并显示
	solution.interp (0.15);         //执行Akima方法
	solution.output ();//输出三次多项式系数以及t点的近似值z到文件并显示	  
}
*/

double Assum(vector<double> &vec)
{
	std::size_t len = vec.size();
	double sum = 0;
	for (int i = 0; i < len; i++)
	{
	 sum += vec[i];
	}

	return (static_cast<double>(sum) / len);
}

void Seprate2GetTrans(const string recorder_str, 
	                     Matrix4d &rotation, 
					     Vector3d &Origanl_XYZ,
				    const string previous_time,
				    const string current_time)
{
	int a[10] = {0};
	cutstring(recorder_str,a);
	Vector4d q;
	Vector3d vel;
/*
    //////////////////////////////////////////////////////////		 
		 2𝑞02+2𝑞12−1  2𝑞1𝑞2−2𝑞0𝑞3   2𝑞1𝑞3+2𝑞0𝑞2 
		 2𝑞1𝑞2+2𝑞0𝑞3  2𝑞02+2𝑞22−1   2𝑞2𝑞3−2𝑞0𝑞1
		 2𝑞1𝑞3−2𝑞0𝑞2  2𝑞2𝑞3+2𝑞0𝑞1   2𝑞02+2𝑞32−1
*/
	q(0) = string2double(recorder_str.substr(a[0]+1, a[1] - a[0] - 1).c_str());
	q(1) = string2double(recorder_str.substr(a[1]+1, a[2] - a[1] - 1).c_str());
	q(2) = string2double(recorder_str.substr(a[2]+1, a[3] - a[2] - 1).c_str());
	q(3) = string2double(recorder_str.substr(a[3]+1, a[4] - a[3] - 1).c_str());

    int delta_time = string2int(current_time.substr(4).c_str()) - string2int(previous_time.substr(4).c_str());

	vel(0) =  string2double(recorder_str.substr(a[4]+1, a[5] - a[4] - 1).c_str());
	vel(1) =  string2double(recorder_str.substr(a[5]+1, a[6] - a[5] - 1).c_str());
	vel(2) =  string2double(recorder_str.substr(a[6]+1).c_str());

	Origanl_XYZ += static_cast<double>(delta_time)*vel / 3 ;
	double  w = q(0), x = q(1) , y = q(2) , z = q(3);

	rotation(0,0) = 2 * ( pow(q(0),2) + pow(q(1),2) ) - 1;
	rotation(0,1) = 2 * q(1) * q(2) - 2 * q(0) * q(3);
	rotation(0,2) = 2 * q(1) * q(3) + 2 * q(0) * q(2);
	rotation(0,3) = Origanl_XYZ(0);

	rotation(1,0) = 2 * q(1) * q(2) + 2 * q(0) * q(3);
	rotation(1,1) = 2 * ( pow(q(0),2) + pow(q(2),2) ) - 1;
	rotation(1,2) = 2 * q(2) * q(3) - 2 * q(0) * q(1);
	rotation(1,3) = Origanl_XYZ(1);

	rotation(2,0) = 2 * q(1) * q(3) - 2 * q(0) * q(2);
	rotation(2,1) = 2 * q(2) * q(3) + 2 * q(0) * q(1);
	rotation(2,2) = 2 * ( pow(q(0),2) + pow(q(3),2) ) - 1;
	//rotation(2,3) = Origanl_XYZ(2);
	rotation(2,3) = 0;
	

  /*
    | 1−2(y2+z2)  2(xy+zw)      2(xz−yw)
	| 2(xy−zw)    1−2(x2+z2)   2(yz+xw)
	| 2(xz+yw)     2(yz−xw)    1−2(x2+y2)
   */
	/*
	rotation(0,0) = - 2 * ( pow(q(2),2) + pow(q(3),2) ) + 1;
	rotation(0,1) = 2 * q(1) * q(2) - 2 * q(0) * q(3);
	rotation(0,2) = 2 * q(1) * q(3) - 2 * q(0) * q(2);
	rotation(0,3) = Origanl_XYZ(0);

	rotation(1,0) = 2 * q(1) * q(2) - 2 * q(0) * q(3);
	rotation(1,1) = - 2 * ( pow(q(1),2) - pow(q(3),2) ) + 1;
	rotation(1,2) = 2 * q(2) * q(3) + 2 * q(0) * q(1);
	rotation(1,3) = Origanl_XYZ(1);

	rotation(2,0) = 2 * q(1) * q(3) + 2 * q(0) * q(2);
	rotation(2,1) = 2 * q(2) * q(3) + 2 * q(0) * q(1);
	rotation(2,2) = - 2 * ( pow(q(1),2) + pow(q(2),2) ) + 1;
	rotation(2,3) = Origanl_XYZ(2);
	//rotation(2,3) = 0;
	*/
	////////////////////////////////////
	/*
	r(1,1)=1-2*y*y-2*z*z;  
	r(1,2)=2*x*y+2*w*z;  
	r(1,3)=2*x*z-2*w*y;  

	r(2,1)=2*x*y-2*w*z;  
	r(2,2)=1-2*x*x-2*z*z;  
	r(2,3)=2*z*y+2*w*x;  

	r(3,1)=2*x*z+2*w*y;  
	r(3,2)=2*y*z-2*w*x;  
	r(3,3)=1-2*x*x-2*y*y;  
	*/
	/*
	rotation(0,0) = 1 - ( 2*y*y + 2*z*z);  
	rotation(0,1) = 2*x*y - 2*w*z;  
	rotation(0,2) = 2*x*z - 2*w*y;  
	rotation(0,3) = Origanl_XYZ(0);

	rotation(1,0) = 2*x*y - 2*w*z;  
	rotation(1,1) = 1 + ( 2*x*x - 2*z*z );  
	rotation(1,2) = 2*z*y + 2*w*x;  
	rotation(1,3) = Origanl_XYZ(1);

	rotation(2,0) = 2*x*z + 2*w*y;  
	rotation(2,1) = 2*y*z + 2*w*x;  
	rotation(2,2) = 1 - ( 2*x*x + 2*y*y );  
	rotation(2,3) = Origanl_XYZ(2);
	//rotation(2,3) = 0;
	
	rotation(3,0) = 0;
	rotation(3,1) = 0;
	rotation(3,2) = 0;
	rotation(3,3) = 1;
	*/
}

void Seprate2GetPCvec(const string str, Vector4d &point,const int signal[5],const Vector3d &coor)
{
	point(0) = string2double(str.substr(0, signal[0]).c_str());
	point(1) = string2double(str.substr(signal[0] + 1, signal[1] - signal[0] -1).c_str());
    point(2) = 0;
 // point(2) = coor(2);
	point(3) = 1;
}

void Combine2PCL(std::map<std::string, std::string> &mymap,const char *str)
{
    if (str == NULL)
	{
	    std::cout << __FILE__<<":"<<__LINE__ <<"  File ERROR!"<<std::endl;
	    return ;
	}

	std::ifstream input(str);
	int number = -1;
	if (input.is_open())
	{
		std::string str = "";
		std::string prefix_timestamp = "";

		pcl::PointCloud<pcl::PointXYZI> TotalCloud;
		pcl::PointCloud<pcl::PointXYZI> midcloud;
		
		Matrix4d rotation;//旋转变换的矩阵，为4*4的矩阵，多一维是为了添加进入平移变换。
		Vector4d f; f(3) = 1; //每一个点在局部坐标系中的位置
		Vector3d current_position ;//第k帧时刻的坐标系距离原点的位置
		current_position(0) = current_position(1) = current_position(2) = 0.0;
		long int i = 0;

		string current_timestamp;
		string previous_timestamp;

		while ( getline(input,str) )
		{
			if (str.find("x86isnice") != std::string::npos)
			{
				number++;
				TotalCloud += midcloud;
				std::size_t count  = string2int(str.substr(0,str.find_first_of(" ")));
			
				midcloud.clear();
				pcl::PointCloud<pcl::PointXYZI>().swap(midcloud);

				i = 0;
				midcloud.height = 1;
				midcloud.width = count;
				midcloud.is_dense = true;
				midcloud.resize(midcloud.width  *  midcloud.height);

				getline(input ,str);
			
				int signal[5] = {0};
				cutstring(str,signal);
				prefix_timestamp = str.substr(signal[1]+1,signal[2] - signal[1] - 1);
				string  recorder_str = mymap[prefix_timestamp];
				/*分割出参数，获得转换矩阵的更新*///Seprate2GetTrans(string recorder_str, Vector4d transform, Vector3d Origanl_XYZ,string previous_time,string current_time)
				                              //Seprate2GetPCvec(string str, Vector4d point,int signal[5])
				if ( 0 == number )
				{
				   current_timestamp = previous_timestamp = prefix_timestamp;
				}
				else
				{
				   previous_timestamp = current_timestamp;
				   current_timestamp = prefix_timestamp;
				}

				Seprate2GetTrans(recorder_str,rotation, current_position, previous_timestamp,current_timestamp);
				/*分割出参数，获得转换矩阵的更新*/
				//Seprate2GetPCvec( str, f, signal);//QQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
				Seprate2GetPCvec( str, f, signal,current_position);//QQ

				f = rotation * f;
				midcloud.points[i].x = f(0);
				midcloud.points[i].y = f(1);
				midcloud.points[i].z = f(2);
				midcloud.points[i].intensity = string2double(str.substr(signal[2] + 1 , signal[3] - signal[2] - 1));
				i++;
			}
			else
			{   
				int signal[5] = {0}; 
				cutstring(str,signal);
				//Seprate2GetPCvec( str, f, signal);QQXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
				Seprate2GetPCvec( str, f, signal,current_position);//QQ
				f = rotation * f;
				midcloud.points[i].x = f(0);
				midcloud.points[i].y = f(1);
				midcloud.points[i].z = f(2);
				midcloud.points[i].intensity = string2double(str.substr(signal[2] + 1 , signal[3] - signal[2] - 1));
				i++;
			}
		}

		TotalCloud += midcloud;

		ofstream outpcd("../mypcd.pcd");

		outpcd << "# .PCD v.7 - Point Cloud Data file format"<<std::endl;
		outpcd << "VERSION 0.7" << std::endl;
		outpcd << "FIELDS x y z intensity" << std::endl;
		outpcd << "SIZE 4 4 4 4" <<std::endl;
		outpcd << "TYPE F F F F" <<std::endl;
		outpcd << "COUNT 1 1 1 1" <<std::endl;
		outpcd << "WIDTH " << TotalCloud.size() << std::endl;
		outpcd << "HEIGHT 1" <<std::endl;
		outpcd << "VIEWPOINT 0 0 0 1 0 0 0" <<std::endl;
		outpcd << "POINTS " << TotalCloud.size() << std::endl;
		outpcd << "DATA ascii" << std::endl;

		for (std::size_t pp = 0; pp < static_cast<std::size_t>(TotalCloud.size()); pp++)
		{
		     outpcd << TotalCloud.points[pp].x<<"  "<< TotalCloud.points[pp].y << "  " << TotalCloud.points[pp].z<< "  "<< TotalCloud.points[pp].intensity << std::endl;
		}
		outpcd.close();
	//	pcl::io::savePCDFile("../pp.pcd",TotalCloud);

		input.close();
		std::cerr<<"Research completed!  "<<TotalCloud.size() << std::endl;
	}
	else
	{
		std::cerr << "One of the file is not Open ,please check!" << endl; 
		exit(EXIT_FAILURE);
	}
}

#if 1

int main()//动态规划匹配
{
	 ifstream xsens;
	 ifstream laser;
	 std::map<string , string> recorder;

	 clock_t start , finish;

	 string laserfile = "E:/binding/1450857349-laser.txt";
	 string xsensfile = "E:/binding/1450857349-Xsnes.txt";
	// xsens.open("D:\\1448246245-Xsnes.txt");
	// laser.open("D:\\test.txt");
	 laser.open(laserfile);
	 xsens.open(xsensfile);

	 if (laser.is_open() && xsens.is_open())
	 {
		 start = clock();
	     string str;
		 vector<string> log;

		 vector<string> xsens_log; //record th log of xsens to interplot value @ "ERROR"

		 while ( getline(laser,str) )
		 {
		    if (str.find("x86isnice") != std::string::npos)
			{
				string xx;
				getline(laser,str);
				int laser_signal[5] = {0};
				cutstring(str,laser_signal);
				int laser_result = string2int(str.substr(laser_signal[1]+5,laser_signal[2] - laser_signal[1] - 5));//@@@@@@@@@@@@@@@@@@@@
			//	int laser_result = string2int(str.substr(laser_signal[1]+4,laser_signal[2] - laser_signal[1] - 5));
				int xsens_signal[10] = {0};
				int xsens_result = 0;
				int count = 0;

			    bool flag = 0;
		
				if ( !log.empty() )
				{
				  for (vector<string>::iterator it = log.begin(); it != log.end(); it++)
				  {  
					   cutstring(*it ,xsens_signal);
					   xsens_result = string2int(it->substr( 4, 10));
					   if ( abs(laser_result - xsens_result) <= ERROR_RANGE)
					   {
					    // recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1), it->substr(0,13)));
						   recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1), *it));
						   flag = 1;
						  break;
					   }
				  }
	              if ( flag == 1)
				  {
				//	  log.clear();
				//	  vector<string>().swap(log);
					  
					  vector<string>::iterator it = log.begin();
					  vector<string>::iterator xl = xsens_log.begin();
					  //while (it != log.end() && it->compare(xx.substr(0,13)) < 0)
					 // while (it != log.end() && it->compare(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1)) < 0)
					  while (it != log.end() && (it->substr(0,13)).compare(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1)) < 0)
					  {//(it->substr(0,13)).compare(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1);
						  it = log.erase(it);
						  xl = xsens_log.erase(xl);
					  }
					  //xsens_log to clear-------------------------------------------//
				  }
				}

				if ( flag == 0 )
				{	
					bool tag = 0;
					do 
					{
						count++;
					  
						getline( xsens, xx);
						
						if ( xx == "")
						{
							tag = 1;
							break;
						}

						cutstring( xx, xsens_signal);
						//log.push_back(xx.substr(0,13));
						log.push_back(xx);
						xsens_log.push_back(xx);//xsens_log的记录，并查集合//////////////////////////////////////////////

						xsens_result = string2int(xx.substr( 4, 10));
					} while (abs(laser_result - xsens_result) > ERROR_RANGE  && !xsens.eof() && count < COUNT_RANGE);

					//	if (count < COUNT_RANGE )
					if (count < COUNT_RANGE && tag != 1)			
					{  
					   //while (it != log.end() && it->compare(xx.substr(0,13)) < 0)
					    vector<string>::iterator it = log.begin();	
						vector<string>::iterator xl = xsens_log.begin();
						while (it != log.end() && it->compare(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1)) < 0)
						{
								it = log.erase(it);
						        xl = xsens_log.erase(xl);//////////////////////////////并查集合剔除
						}
					//   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
					//	recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1), xx.substr(0,13)));
					//   recorder.insert(pair<string,string>(str, xx));
						recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1), xx));
					}
					else
					{
						TTTTT++;
                     #if 0
						recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1),"ERROR"));
                     #else  
                        /////////////
						
					   	 for (int cheap = 0 ; cheap <= 3 && !xsens.eof(); cheap++)
						 {
						 getline( xsens, xx);

							 if ( xx == "")
							 {
								 break;
							 }

						 cutstring( xx, xsens_signal);
						 log.push_back(xx);

						 xsens_log.push_back(xx);//xsens_log的记录，并查集合//////////////////////////////////////////////
						 }
						
						////////////
						vector<vector<double>> YY(7);
						vector<double> XX;

						for (vector<string>::iterator i = xsens_log.begin(); i != xsens_log.end(); i++)
						{
						  int pp[10] = {0};
						  ::cutstring(*i , pp);

						  XX.push_back(static_cast<double>(::string2int(i->substr(4,10).c_str())));
						
						  YY[0].push_back(string2double(i->substr(pp[0]+1, pp[1] - pp[0] - 1).c_str()));
						  YY[1].push_back(string2double(i->substr(pp[1]+1, pp[2] - pp[1] - 1).c_str()));
						  YY[2].push_back(string2double(i->substr(pp[2]+1, pp[3] - pp[2] - 1).c_str()));
						  YY[3].push_back(string2double(i->substr(pp[3]+1, pp[4] - pp[3] - 1).c_str()));
						  YY[4].push_back(string2double(i->substr(pp[4]+1, pp[5] - pp[4] - 1).c_str()));
						  YY[5].push_back(string2double(i->substr(pp[5]+1, pp[6] - pp[5] - 1).c_str()));
						  YY[6].push_back(string2double(i->substr(pp[6]+1).c_str()));
						}//遍历获得数值

						double interplot[7] = {0.0};
						int num_of_element = XX.size();
						for (int i = 0; i < 7; i++)
						{
						  akima computing(num_of_element);

						  computing.input(XX,YY[i]);

						  double x = string2double(str.substr(laser_signal[1]+5,laser_signal[2] - laser_signal[1] - 5));
						  interplot[i] = computing.interp(x);
						  if (_isnan(interplot[i]))
						  {
						     interplot[i] = Assum(YY[i]);
						  }
						//  computing.~akima();
						}
					//	std::cout << interplot[0]<<" "<< interplot[1]<<" "<<interplot[2]<<" "<< interplot[3]<<" " << interplot[4]<<" " << interplot[5]<<" "<<interplot[6]<<endl;
						string result = str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1) +":"+ 
							                                               double2String(interplot[0]) +":"+
																		   double2String(interplot[1]) +":"+
																		   double2String(interplot[2]) +":"+
																		   double2String(interplot[3]) +":"+
																		   double2String(interplot[4]) +":"+
																		   double2String(interplot[5]) +":"+
																		   double2String(interplot[6]) ;
					  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
					  recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1),result));
				      //recorder.insert(pair<string,string>(str,result + " == x86isnice"));
						////////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
					 //   recorder.insert(pair<string,string>(str.substr(laser_signal[1]+1,laser_signal[2] - laser_signal[1] - 1),result + " == x86isnice"));

							
                     #endif
					}
				}

				str.clear();
			}
			else
			{  
				str.clear();
				continue;
			}
		 }

	 }
	 else
	 {
		 std::cerr << "One of the file is not Open ,please check!" << endl; 
		 exit(EXIT_FAILURE);
	 }

	 laser.close();
	 xsens.close();
	 
	 for (std::map<string,string>::iterator it = recorder.begin(); it != recorder.end(); it++)
	 {
	    std::cout << it->first << " <--> " << it->second << endl;
	 }
	 
	 finish = clock();
	 cout <<"recorder.size()::"<< recorder.size()<<"   "<<"NOT HIT:: " << TTTTT<<"  OUT_OF_RANGE:"<<KKKKK/7<<"  The Cost Of Time is:"<<finish - start<<" ms" <<endl;
	 
     Combine2PCL(recorder, laserfile.c_str());

	 std::cout << "Complete All" << std::endl;
	 return 0;
}

#else 
int main()
{
	ifstream laser;
	ofstream output("D:\\test.txt");
	laser.open("D:\\1448246245-laser.txt");

	if (laser.is_open())
	{
		string str;
		while ( getline(laser,str)  && output.good())
		{
			if (str.find("x86isnice") != std::string::npos)
			{
				output << str <<endl;

				getline(laser,str);

				output << str<< endl;
/*
				getline(laser,str);	
				output << str << endl;

				getline(laser,str);
				output << str << endl;
*/
			}
			else
				continue;
		}
	}
	else
	{
		std::cerr << "laser is not open!" <<endl;
		exit(EXIT_FAILURE);
	}

	laser.close();
	output.close();
   return 0;
}
#endif