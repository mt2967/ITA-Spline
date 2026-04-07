//Author: Kazuya Mori
//Original development : 1990s
//Public release : 2026
// 
// 
// Vct3.cpp: Vct3 クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
//////#include "stdafx.h"
#include "pch.h"  



#include "Vct3.h"
//#include <iostream>
#include "math.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

Vct3::Vct3()
{
	x=0.0;
	y=0.0;
	z=0.0;
}

Vct3::~Vct3()
{

}

double Vct3::absolute() const
{
	return(sqrt(x*x+y*y+z*z));
}
double Vct3::length() const
{
	return(this->absolute());
}

bool Vct3::get_norm_vector(Vct3& w) const
{
 
 	double r=sqrt(x*x+y*y+z*z);
	  
	if ( fabs(r) <= Vct3Const::ZERO_TOLERANCE ) {
		return(false);
	}
 
	w.x=x/r;
	w.y=y/r;
	w.z=z/r;

	return(true);
}

double Vct3::operator|(Vct3 ob2) const      // ベクトルどおしの内積
{
  	double r=x*ob2.x+y*ob2.y+z*ob2.z;
  	return(r);
}

double Vct3::inner_product(const Vct3 v, const Vct3 w)
{
 
  	double r=v.x*w.x+v.y*w.y+v.z*w.z;
  	return(r);
 
}

Vct3 Vct3::operator^(Vct3 ob2) const     // ベクトルどおしの外積
{
 	Vct3 g;
  	g.x=y*ob2.z-z*ob2.y;
 	g.y=z*ob2.x-x*ob2.z;
  	g.z=x*ob2.y-y*ob2.x;
  
  	return(g);
}

Vct3 Vct3::outer_product(const Vct3 v,const Vct3 w)
{
 
 	Vct3 g;
  	g.x=v.y*w.z-v.z*w.y;
 	g.y=v.z*w.x-v.x*w.z;
  	g.z=v.x*w.y-v.y*w.x;
  
  	return(g);
 
}

Vct3::Vct3(double a, double b, double c)
{
	x=a;
	y=b;
	z=c;
}

Vct3 Vct3::operator/(double s) const
{

	if (fabs(s) <Vct3Const::ZERO_TOLERANCE ) {return *this;}

    Vct3 tmp;
    tmp.x=x/s;
    tmp.y=y/s;
    tmp.z=z/s;

    return(tmp);
}

Vct3 Vct3::operator +(Vct3 ob2) const
{
  Vct3 tmp;
  tmp.x=x+ob2.x;
  tmp.y=y+ob2.y;
  tmp.z=z+ob2.z;

  return(tmp);
}

Vct3 Vct3::operator -(Vct3 ob2) const
{
  Vct3 tmp;
  tmp.x=x-ob2.x;
  tmp.y=y-ob2.y;
  tmp.z=z-ob2.z;

  return(tmp);

}

Vct3 Vct3::operator +() const
{
	return *this;
}

Vct3 Vct3::operator -() const
{
	Vct3 tmp;
	tmp.x=-x;
	tmp.y=-y;
	tmp.z=-z;
	return(tmp);

}

bool Vct3::set(const double buf[3])
{

	x=buf[0];
	y=buf[1];
	z=buf[2];

	return(true);

} 

bool Vct3::set(const double a,const double b,const double c)
{

	x=a;
	y=b;
	z=c;

	return(true);

} 
 
bool Vct3::get(double buf[3])	//取得
{
	buf[0]=this->x;
	buf[1]=this->y;
	buf[2]=this->z;

	return(true);
}

bool Vct3::get(double& a,double& b,double& c)  //取得
{
	a=this->x;
	b=this->y;
	c=this->z;

	return(true);
}

bool Vct3::is_parallel_to(const Vct3 w) const  // 平行か？
{
	// 平行の定義
	// 外積ベクトルの大きさが、０ならば平行とする。
	double len=((*this)^w).length();
	if ( fabs(len) < Vct3Const::ABOUT_ZERO  ) {
		return(true);
	}
	else{
		return(false);
	}
}
bool Vct3::is_orthogonal_to(const Vct3 w) const // 垂直か？
{
	// 垂直 
	double naiseki=(*this)|w;
	if ( fabs(naiseki) < Vct3Const::ABOUT_ZERO ) {
		return(true);
	}
	else{
		return(false);
	}
}
	//以下、フレンド関数
 

Vct3 operator *(double s,Vct3 v)
{
  Vct3 tmp;
  tmp.x=s*v.x;
  tmp.y=s*v.y;
  tmp.z=s*v.z;

  return(tmp);

}