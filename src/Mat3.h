//Author: Kazuya Mori
//Original development : 1990s
//Public release : 2026
// 
// Mat3.h: Mat3 クラスのインターフェイス
//
//
// 使用方法例
//	 
//+		Mat3 A;  //インスタンス発生（=変数の初期化）
//+		A.set_axis_y_rot(30);	//Y軸回りの30度回転行列のセット
//+
//+		//行列の行列式
//+		double det=A.det();
//+		cout<<"行列式＝"<<det<<endl;
//+
//+		if (fabs(det)>0.01){
//+			Mat3 C=A.Inverse();  //逆行列
//+			A.print("A=");
//+			C.print("InvA=");
//+			(A*C).print("単位行列？＝");
//+		}
//
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
#include "pch.h"  


#include "Vct3.h"

class Mat3  
{
public:
	double det();			//行列式
	double Trace();			//対角和
	double m_11;				//(1,1)成分
	double m_12;				//(1,2)成分
	double m_13;				//(1,3)成分
	double m_21;				//(2,1)成分
	double m_22;				//(2,2)成分
	double m_23;				//(2,3)成分
	double m_31;				//(3,1)成分
	double m_32;				//(3,2)成分
	double m_33;				//(3,3)成分

	Mat3();
	Mat3(double f[9]);	//設定付きコンストラクタ
 
	Mat3(double a11,double a12,double a13,   //<---行列の1行目（横）
		   double a21,double a22,double a23,   //<---行列の2行目（横）
	       double a31,double a32,double a33);  //<---行列の3行目（横）
	virtual ~Mat3();

	Mat3 operator+(Mat3 ob2);       // 行列和
	Mat3 operator-(Mat3 ob2);		// 行列差
	Mat3 operator+();					// +単項演算子  
	Mat3 operator-();					// -単項演算子
	Mat3 operator*(Mat3 ob2);		// 行列積
	Vct3 operator*(Vct3 v);			// 行列と縦ベクトルの積 Mv
	friend Mat3 operator*(double s,Mat3 ob2);	// スカラ積  s*ob2
	 
    bool  set_with_3column_vector(
		double tate1[3],
		double tate2[3],
		double tate3[3]);// 3個の列ベクトルによる設定

    bool  set_with_3column_vector(
		const Vct3 f1,
		const Vct3 f2,
		const Vct3 f3);// 3個の列ベクトルによる設定	

    bool  set_SO3_with_3column_vector(
		double tate1[3],
		double tate2[3],
		double tate3[3]);// 3個の列ベクトルによる直交行列の設定
	
    bool  set_SO3_with_2column_vector(
		double tate1[3],double tate2[3]);// 2個の列ベクトルによる直交行列の設定

	bool is_SO3_matrix();	//回転行列か？

    Mat3 Transposed();				// 転置行列
    Mat3 Inverse();					// 逆行列

	// 基本軸に関する回転ベクトルの設定（角度の単位はdegreeで）
	bool set_axis_x_rot(const double deg); // X軸回転行列の設定
	bool set_axis_y_rot(const double deg); // Y軸回転行列の設定
	bool set_axis_z_rot(const double deg); // Z軸回転行列の設定

};
