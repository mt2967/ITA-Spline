//Author: Kazuya Mori
//Original development : 1990s
//Public release : 2026
// 
// Vct3.h: Vct3 クラスのインターフェイス
//
// 
//
// 使用方法例
//
//+		// 平面πと点(q)の距離
//+		//
//+		//   点pを通り、単位法線ベクトルnを持つ平面の方定式  
//+		//      <x-p,n>=0      <*,*>は内積です。
//+		//   
//+		//   点qとの距離(d)は
//+		//     d=|<p-q,n>|
//+		//   交点ベクトルeは
//+		//     e=q+<p-q,n>n
//+		//
//+		Vct3 p(10,20,30);
//+		Vct3 n(0,1,0);
//+		Vct3 q(500,600,700);
//+		
//+		double d=fabs((p-q)|n) ;    //d=|<p-q,n>|
//+		cout<<"距離＝"<<d<<endl;
//+
//+		Vct3 e;
//+		e=q+((p-q)|n)*n;			//e=q+<p-q,n>n
//+		e.print("交点＝");
//
//////////////////////////////////////////////////////////////////////
//
#include <vector>
#include <algorithm>

namespace Vct3Const {
	const double PI     = 3.1415926535897932;
	const double PI_2   = 1.5707963267948966;
	const double PI_4   = 0.7853981633974483;
	const double PI_180 = 0.0174532925199432957;   //  PAI/180.0 ラジアンにする
	const double ZERO_TOLERANCE = 0.000001;
	const double ABOUT_ZERO = 0.001;
} 

class Vct3  
{

public:	
	double x;
	double y; 
	double z; 

	Vct3();
	Vct3(double a,double b, double c);			//設定付きコンストラクタ
	virtual ~Vct3();
	
	//////int print(const std::string name);			//coutで出力
	bool set(const double buf[3]);				//設定
	bool set(const double a,const double b,const double c); //設定 
	bool get(double  buf[3]);					//取得
	bool get(double& a,double& b,double& c);		//取得

	Vct3 operator/(double s) const;				// スカラ割り算 例. v/5.0
	Vct3 operator+(Vct3 ob2) const;			// ベクトル和
	Vct3 operator-(Vct3 ob2) const;			// ベクトル引算
	Vct3 operator^(Vct3 ob2) const;			// 外積
	Vct3 operator+() const;					//+単項演算子  
	Vct3 operator-() const ;					//-単項演算子
	double absolute() const;						// 絶対値
	double length() const;						// 長さ(絶対値と同じ）
	double operator|(Vct3 ob2) const;			// 内積
	bool get_norm_vector(Vct3& w) const;		// 規格化ベクトル(自分自身は不変）
	bool is_parallel_to(const Vct3 w) const;	// 平行か？
	bool is_orthogonal_to(const Vct3 w) const;	// 直交か？

	friend Vct3 operator*(double s,Vct3 v);	// スカラ積  s*v

	static Vct3 outer_product(const Vct3 v,const Vct3 w);  //外積ベクトル
	static double inner_product(const Vct3 v,const Vct3 w);  //内積（スカラー）
 
	typedef  std::pair<Vct3,Vct3> Line;
	typedef  std::vector<Line> Lines; 
	typedef  std::vector<Lines> vecLines;	

};
