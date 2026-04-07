//Author: Kazuya Mori
//Original development : 1990s
//Public release : 2026
// 
// Mat3.cpp: Mat3 クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
////#include "stdafx.h"
#include "pch.h"  


#include <iostream>

#include "Mat3.h"
#include "math.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

Mat3::Mat3()
{
	 m_11=0.0;
	 m_12=0.0;
	 m_13=0.0;
	 m_21=0.0;
	 m_22=0.0;
	 m_23=0.0;
	 m_31=0.0;
	 m_32=0.0;
	 m_33=0.0;
}

Mat3::~Mat3()
{

}

Mat3::Mat3(double f[9])
{
	 m_11=f[0];
	 m_12=f[1];
	 m_13=f[2];
	 m_21=f[3];
	 m_22=f[4];
	 m_23=f[5];
	 m_31=f[6];
	 m_32=f[7];
	 m_33=f[8];	
}

 
Mat3::Mat3(double a11,double a12,double a13,	//<---行列の1行目（横）
		       double a21,double a22,double a23,	//<---行列の2行目（横）
	           double a31,double a32,double a33)	//<---行列の3行目（横）
{
	 m_11=a11;
	 m_12=a12;
	 m_13=a13;
	 m_21=a21;
	 m_22=a22;
	 m_23=a23;
	 m_31=a31;
	 m_32=a32;
	 m_33=a33;	
}



Mat3 Mat3::operator +(Mat3 ob2)
{
  Mat3 tmp;

  tmp.m_11=m_11+ob2.m_11;
  tmp.m_12=m_12+ob2.m_12;
  tmp.m_13=m_13+ob2.m_13;
  tmp.m_21=m_21+ob2.m_21;
  tmp.m_22=m_22+ob2.m_22;
  tmp.m_23=m_23+ob2.m_23;
  tmp.m_31=m_31+ob2.m_31;
  tmp.m_32=m_32+ob2.m_32;
  tmp.m_33=m_33+ob2.m_33;

  return(tmp);

}

Mat3 Mat3::operator -(Mat3 ob2)
{
  Mat3 tmp;

  tmp.m_11=m_11-ob2.m_11;
  tmp.m_12=m_12-ob2.m_12;
  tmp.m_13=m_13-ob2.m_13;
  tmp.m_21=m_21-ob2.m_21;
  tmp.m_22=m_22-ob2.m_22;
  tmp.m_23=m_23-ob2.m_23;
  tmp.m_31=m_31-ob2.m_31;
  tmp.m_32=m_32-ob2.m_32;
  tmp.m_33=m_33-ob2.m_33;

  return(tmp);

}


Mat3 Mat3::operator +()
{
 
	return *this;

}

Mat3 Mat3::operator -()
{

	 Mat3 tmp;

	 tmp.m_11=-m_11;
	 tmp.m_12=-m_12;
	 tmp.m_13=-m_13;
	 tmp.m_21=-m_21;
	 tmp.m_22=-m_22;
	 tmp.m_23=-m_23;
	 tmp.m_31=-m_31;
	 tmp.m_32=-m_32;
	 tmp.m_33=-m_33;
	 
	 return(tmp);

}


//	Mat3 operator*(Vct3 v);			// 行列と縦ベクトルの積 Mv
Vct3 Mat3::operator*(Vct3 v)
{
  Vct3 tmp;

  tmp.x=m_11*v.x +m_12*v.y+m_13*v.z;
  tmp.y=m_21*v.x +m_22*v.y+m_23*v.z;  
  tmp.z=m_31*v.x +m_32*v.y+m_33*v.z;  

  return(tmp);

}

Mat3 Mat3::Transposed()	// 転置行列
{
  Mat3 tmp;

  tmp.m_11=m_11;
  tmp.m_12=m_21; 
  tmp.m_13=m_31;  
  tmp.m_21=m_12;  
  tmp.m_22=m_22;  
  tmp.m_23=m_32;  
  tmp.m_31=m_13;  
  tmp.m_32=m_23;  
  tmp.m_33=m_33;  

  return(tmp);
}


Mat3 Mat3::operator *(Mat3 ob2)
{
  Mat3 tmp;

  tmp.m_11=m_11*ob2.m_11+m_12*ob2.m_21+m_13*ob2.m_31;
  tmp.m_12=m_11*ob2.m_12+m_12*ob2.m_22+m_13*ob2.m_32; 
  tmp.m_13=m_11*ob2.m_13+m_12*ob2.m_23+m_13*ob2.m_33;  
  tmp.m_21=m_21*ob2.m_11+m_22*ob2.m_21+m_23*ob2.m_31;  
  tmp.m_22=m_21*ob2.m_12+m_22*ob2.m_22+m_23*ob2.m_32;  
  tmp.m_23=m_21*ob2.m_13+m_22*ob2.m_23+m_23*ob2.m_33;  
  tmp.m_31=m_31*ob2.m_11+m_32*ob2.m_21+m_33*ob2.m_31;  
  tmp.m_32=m_31*ob2.m_12+m_32*ob2.m_22+m_33*ob2.m_32;  
  tmp.m_33=m_31*ob2.m_13+m_32*ob2.m_23+m_33*ob2.m_33;  

  return(tmp);

}

//以下、フレンド関数
Mat3 operator *(double s,Mat3 ob2)
{
	Mat3 tmp;
 
	tmp.m_11=s*ob2.m_11;
	tmp.m_12=s*ob2.m_12;
	tmp.m_13=s*ob2.m_13;
	tmp.m_21=s*ob2.m_21;
	tmp.m_22=s*ob2.m_22;
	tmp.m_23=s*ob2.m_23;
	tmp.m_31=s*ob2.m_31;
	tmp.m_32=s*ob2.m_32;
	tmp.m_33=s*ob2.m_33;

	return(tmp);

}

double Mat3::Trace()  //トレース（対角成分の和）
{
	double f=m_11+m_22+m_33;
	return(f);
}


double Mat3::det()  //行列式
{
	double f=+m_11*m_22*m_33
			+m_12*m_23*m_31
			+m_13*m_21*m_32
			-m_11*m_23*m_32
			-m_12*m_21*m_33
			-m_13*m_22*m_31;
	return(f);
}
 
Mat3 Mat3::Inverse()	// 逆行列
{
  Mat3 tmp;
	
  double determinant= this->det() ;
  if ( fabs(determinant) <Vct3Const::ZERO_TOLERANCE ){		//行列式が０
	  return tmp;	//赤バッテンエラーでもよい 
  }  
 
  // 順番非常に注意！！（余因子の関係）
  tmp.m_11=+(m_22*m_33-m_23*m_32)/determinant;
  tmp.m_21=-(m_21*m_33-m_23*m_31)/determinant; 
  tmp.m_31=+(m_21*m_32-m_22*m_31)/determinant; 
  tmp.m_12=-(m_12*m_33-m_13*m_32)/determinant; 
  tmp.m_22=+(m_11*m_33-m_13*m_31)/determinant; 
  tmp.m_32=-(m_11*m_32-m_12*m_31)/determinant; 
  tmp.m_13=+(m_12*m_23-m_13*m_22)/determinant; 
  tmp.m_23=-(m_11*m_23-m_13*m_21)/determinant; 
  tmp.m_33=+(m_11*m_22-m_12*m_21)/determinant;   
  
  return(tmp);
}


bool Mat3::set_with_3column_vector
	(double tate1[3],double tate2[3],double tate3[3]) 	// 3個の列ベクトルによる設定
{
	this->m_11=tate1[0];
	this->m_21=tate1[1];
	this->m_31=tate1[2];

	this->m_12=tate2[0];
	this->m_22=tate2[1];
	this->m_32=tate2[2];
	
	this->m_13=tate3[0];
	this->m_23=tate3[1];
	this->m_33=tate3[2];

	return (0);
}

bool Mat3::set_axis_x_rot(const double deg)
{
	double radi=deg*Vct3Const::PI_180; //ラジアンにする

	//第1列
	this->m_11=1.0;
	this->m_21=0.0; 
	this->m_31=0.0;

	//第2列
	this->m_12=0.0; 
	this->m_22=cos(radi); 
	this->m_32=sin(radi); 
		
	//第3列
	this->m_13=0.0; 
	this->m_23=-sin(radi); 
	this->m_33=cos(radi); 

	return (0);

}


bool Mat3::set_axis_y_rot(const double deg)
{
	double radi=deg*Vct3Const::PI_180; //ラジアンにする

	//第1列
	this->m_11=cos(radi);
	this->m_21=0.0; 
	this->m_31=-sin(radi);

	//第2列
	this->m_12=0.0; 
	this->m_22=1.0; 
	this->m_32=0.0; 
		
	//第3列
	this->m_13=sin(radi); 
	this->m_23=0.0; 
	this->m_33=cos(radi); 

	return (0);

}


bool Mat3::set_axis_z_rot(const double deg)
{
 
	double radi=deg*Vct3Const::PI_180; //ラジアンにする

	//第1列
	this->m_11=cos(radi);
	this->m_21=sin(radi); 
	this->m_31=0.0;

	//第2列
	this->m_12=-sin(radi); 
	this->m_22=cos(radi); 
	this->m_32=0.0; 
		
	//第3列
	this->m_13=0.0; 
	this->m_23=0.0; 
	this->m_33=1.0; 

	return (0);

}

bool Mat3::set_SO3_with_3column_vector
	(double tate1[3],double tate2[3],double tate3[3]) 	// 3個の列ベクトルによる設定
{

	Mat3 tmp;
	tmp.set_with_3column_vector(tate1 ,tate2 ,tate3 );
	double determinant=tmp.det();

	if ( determinant>0.999 && determinant<1.001 )
	{
		Vct3 v1,v2,v3;
		v1.set(tate1);
		v2.set(tate2);
		v3.set(tate3);

		// この条件だけでは、直交行列とは限らない(? →将来課題）
		double s=((v1^v2)^v3).absolute();  // ＾は外積  
		if ( fabs(s)>Vct3Const::ABOUT_ZERO ) {
			return(false);//直交行列ではない
		}

		this->m_11=tate1[0];
		this->m_21=tate1[1];
		this->m_31=tate1[2];

		this->m_12=tate2[0];
		this->m_22=tate2[1];
		this->m_32=tate2[2];
		
		this->m_13=tate3[0];
		this->m_23=tate3[1];
		this->m_33=tate3[2];

		return (0);
	}
	else{
		return(false);  //直交行列ではない
	}

}

bool Mat3::set_SO3_with_2column_vector
	(double tate1[3],double tate2[3]) 	// 2個の列ベクトルによる設定
{

	Vct3 v1,v2,v3;
	v1.set(tate1);
	v2.set(tate2);
		
	double naiseki=v1|v2;
	if (fabs(naiseki)>Vct3Const::ABOUT_ZERO){
		return(false);//直交ベクトルではない
	}

	Vct3 n1,n2,n3;
	if ( v1.get_norm_vector(n1) !=0 ) {
		return(false);
	}
	if ( v2.get_norm_vector(n2) !=0 ) {
		return(false);
	}
 
	n3=n1^n2;  //外積

	this->m_11=n1.x;
	this->m_21=n1.y;
	this->m_31=n1.z;

	this->m_12=n2.x;
	this->m_22=n2.y;
	this->m_32=n2.z;
		
	this->m_13=n3.x;
	this->m_23=n3.y;
	this->m_33=n3.z;

	return (0);
 
}


bool Mat3::is_SO3_matrix()
{

	double determinant=this->det();
	if ( fabs(determinant-1.0)>Vct3Const::ABOUT_ZERO){
		return false;
	}

	//回転行列の判断  転置行列との積が、単位行列になるか？
	Mat3 A,B,C;
	A=*this;
	B=A.Transposed();
	C=A*B;
	if (fabs(C.m_11-1.0)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_22-1.0)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_33-1.0)>Vct3Const::ABOUT_ZERO)return(false);

	if (fabs(C.m_12)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_21)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_13)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_31)>Vct3Const::ABOUT_ZERO)return(false);	
	if (fabs(C.m_23)>Vct3Const::ABOUT_ZERO)return(false);
	if (fabs(C.m_32)>Vct3Const::ABOUT_ZERO)return(false);		

	return true;
}



bool Mat3::set_with_3column_vector(
		const Vct3 f1,
		const Vct3 f2,
		const Vct3 f3) 	// 3個の列ベクトルによる設定
{
	this->m_11=f1.x;
	this->m_21=f1.y;
	this->m_31=f1.z;

	this->m_12=f2.x;
	this->m_22=f2.y;
	this->m_32=f2.z;
	
	this->m_13=f3.x;
	this->m_23=f3.y;
	this->m_33=f3.z;

	return (0);
}