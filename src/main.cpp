//Author: Kazuya Mori
//Original development : 1990s
//Public release : 2026
//

#include "pch.h" 

#include "Vct3.h"


#include <iostream> 

#include "math.h" 
#include <algorithm>
#include "iostream"

#define PI 3.141592653

// 円弧
struct strct_arc {
	Vct3 Org;			//円弧中心（親円の中心）   
	Vct3 Vm;			//主軸ベクトル（親円の接線すべてに垂直な方向）
	Vct3 Vn;			//円弧補助軸ベクトル（円弧中心点から、弧の始点方向のベクトル　規格化）
	double Ang_rad;     //円弧の中心角度（弧度法）
	double Radius;		//円弧半径
};



Vct3 V_rot(Vct3 v1, Vct3 v2, double theta_RAD) {
	// v2を v1軸回りにθ回転した結果を返す

	Vct3 e1 = v1 / v1.length();
	Vct3 e2 = v2 / v2.length();

	Vct3 vOut = (e1 | e2) * (1 - cos(theta_RAD)) * e1 + cos(theta_RAD) * e2 + sin(theta_RAD) * (e1 ^ e2);
	vOut = v2.length() * vOut;

	return(vOut);

}


int ITA_ArcConnectPoint(Vct3 Pf, Vct3 Vf, Vct3 Pb, Vct3 Vb, Vct3& Pcnt) {
	//(in) Pf,Vf  頂点、方向ベクトル(forward)
	//(in) Pb,Vb  頂点、方向ベクトル(back)
	//(out)Pcnt   ２円弧の接続点


	Vct3 Vd = Pf - Pb;

	double rr = (Vf | Vb);  //内積

	double r2 = (Vd | Vd);

	double t;

	if (rr <= -0.999) {

		double r1 = (Vd | Vb);
		t = r2 / (4.0 * r1);  //r1 = 0は、この場合ありえない
	}
	else {

		Vct3 V = Vf - Vb;

		double r1 = (Vd | V);
		double r3 = 1.0 + rr;

		t = (r1 + sqrt(r1 * r1 + 2.0 * r2 * r3)) / (2.0 * r3);
	}

	Vct3 P1 = Pf + t * Vf;
	Vct3 P2 = Pb + t * Vb;
	Pcnt = (P1 + P2) / 2.0;

	return(0);
}

int ArcConv(Vct3 Ps, Vct3 Vdir, Vct3 Vtan, Vct3& Pt, Vct3& Vm, Vct3& Vn, double& ang_rad, double& radius) {
	//(in)
	//  Ps 円弧の端の点
	//  Vdir Psからもう一方の端の点への方向ベクトル
	//  Vtan Psの接線ベクトル
	//(out)
	//  Pt 円弧の中心 
	//  Vm Vn 主軸、補助軸
	//  ang_rad 円弧の中心角度（ラジアン）
	//  radius 円弧半径[km]
	//
	Vtan = Vtan / Vtan.length();

	double r_1 = (Vdir | Vtan);
	double r_2 = (Vdir | Vdir);
	double r_3 = (Vtan | Vtan);

	if (abs(r_2) < 0.001) {
		ang_rad = 0;
		radius = 0;

		return(0);
	}

	//  (Vdir|Vtan)=|Vdir|*|Vtan|*cos_p_v   
	double cos_p_v = r_1 / (sqrt(r_2) * sqrt(r_3));  //角度は０～１８０度

	cos_p_v = std::min(1.0, cos_p_v);

	ang_rad = 2.0 * acos(cos_p_v);   //   anr_radは、0～2π


	if (ang_rad < 0.001) {
		ang_rad = 0;
		radius = 0;

		return(0);
	}

	double  sin_2 = 1 - cos_p_v * cos_p_v;

	//トーラス（円弧）の中心位置(arc_cent)
	Vct3 Vr = (-r_1 / (2 * sin_2 * r_3)) * Vtan + Vdir / (2.0 * sin_2);

	Pt = Ps + Vr;

	//トーラス（円弧）の曲率半径
	radius = Vr.length();

	Vn = -Vr / Vr.length();
	Vm = (Vn ^ Vtan) / (Vn ^ Vtan).length();

	return(0);
}

int ITA_ArcConv(Vct3 Ps, Vct3 Vdir, Vct3 Vtan, strct_arc& arc) {
	//
	// 円弧端点１、円弧端点２から円弧端点１までの差ベクトル、端点１での接ベクトルから、円弧構造体に変換する
	//(in) Ps 円弧端点１
	//     Vdir 円弧端点２から円弧端点１までの差ベクトル
	//	   Vtan 円弧端点１での接ベクトル（規格化）
	//(out) 円弧構造体
	
	Vct3 _Pt;
	Vct3 _Vm;
	Vct3 _Vn;

	double _ang_rad;
	double _radius;

	int ret = ArcConv(Ps, Vdir, Vtan, _Pt, _Vm, _Vn, _ang_rad, _radius);

	arc.Org = _Pt;
	arc.Vm = _Vm;
	arc.Vn = _Vn;
	arc.Ang_rad = _ang_rad;
	arc.Radius = _radius;

	return(ret);

}

double ArcLength(std::vector<strct_arc> vecArc) {
	//(戻り値）円弧曲線長
	//　 円弧列の長さを計算する
	//
	double len = 0;

	for (int i = 0; i < vecArc.size(); i++) {

		strct_arc arc = vecArc[i];
		len = len + arc.Radius * arc.Ang_rad;
	}

	return(len);
}

int Point_in_Arcs(std::vector<strct_arc> vecArc, double curve_length, Vct3& Pt,double &r1) {
	// 円弧列に対して、始点からの距離を与え、その地点の座標と曲率半径を出力する
	// (in) vecArc 円弧列
	//      curve_length　始点からの距離（全円弧長以下であること）
	// (out)Pt 座標
	//      r1 この座標の曲率半径
	// 
	double len = 0;

	for (int i = 0; i < vecArc.size(); i++) {

		strct_arc arc = vecArc[i];

		len = len + arc.Radius * arc.Ang_rad;

		if (len >= curve_length) {

			double beta = (len - curve_length) / arc.Radius;
			double alfa = arc.Ang_rad;

			Pt = arc.Org + arc.Radius * (cos(alfa - beta))*arc.Vn + arc.Radius * (sin(alfa - beta))*((arc.Vm)^(arc.Vn));
			r1 = arc.Radius;

			goto L100;
			
		}
	}

	//　エラー
	return(-1);

L100:
	return(0);
}


int ITA_Points_to_Arcs(std::vector<Vct3> vecPt,int Type, std::vector <strct_arc> &vecArc) {
	// 空間内の点列からITAスプライン円弧列を取得する。
	//(in) vecPt  3次元点列
	//(in) Type   10=対辺平行法＆開曲線,11=対辺平行法＆閉曲線   30=外心法＆開曲線,31=外心法＆閉曲線
	//(out)vecArc 円弧列

	bool IsCloseCurve = Type % 2;

	int numP = vecPt.size();

	for (int j = 0; j < numP - 1; j++) {


		Vct3 Pf = vecPt[j];
		Vct3 Pb = vecPt[j + 1];

		Vct3 Vf, Vb;

		int FLG = Type / 10;

		if (FLG == 1) {
			//垂心点法（従来）
			Vct3 Pff, Pbb;

			if (j == 0) {

				Pff = vecPt[numP - 2];
				Pbb = vecPt[j + 2];


				if (IsCloseCurve == 0) {
					Vf = Pb - Pf;
					Vb = -(Pbb - Pf);
				}
				else {
					Vf = Pb - Pff;
					Vb = -(Pbb - Pf);
				}

			}
			else if (j == (numP - 2)) {

				Pff = vecPt[j - 1];
				Pbb = vecPt[1];

				if (IsCloseCurve == 0) {
					Vf = Pb - Pff;
					Vb = -(Pb - Pf);
				}
				else {
					Vf = Pb - Pff;
					Vb = -(Pbb - Pf);
				}

			}
			else {

				Pff = vecPt[j - 1];
				Pbb = vecPt[j + 2];

				Vf = Pb - Pff;
				Vb = -(Pbb - Pf);

			}
		}
		else if (FLG == 2) {
			//重心点法
			Vct3 Pff, Pbb;

			if (j == 0) {

				Pff = vecPt[numP - 2];
				Pbb = vecPt[j + 2];


				if (IsCloseCurve == 0) {
					Vf = Pb - Pf;
					Vb = -((Pbb - Pb) / (Pbb - Pb).length() - (Pf - Pb) / (Pf - Pb).length());
				}
				else {
					Vf = (Pb - Pf) / (Pb - Pf).length() - (Pff - Pf) / (Pff - Pf).length();
					Vb = -((Pbb - Pb) / (Pbb - Pb).length() - (Pf - Pb) / (Pf - Pb).length());
				}

			}
			else if (j == (numP - 2)) {

				Pff = vecPt[j - 1];
				Pbb = vecPt[1];

				if (IsCloseCurve == 0) {
					Vf = (Pb - Pf) / (Pb - Pf).length() - (Pff - Pf) / (Pff - Pf).length();
					Vb = -(Pb - Pf);
				}
				else {
					Vf = (Pb - Pf) / (Pb - Pf).length() - (Pff - Pf) / (Pff - Pf).length();
					Vb = -((Pbb - Pb) / (Pbb - Pb).length() - (Pf - Pb) / (Pf - Pb).length());
				}

			}
			else {

				Pff = vecPt[j - 1];
				Pbb = vecPt[j + 2];

				Vf = (Pb - Pf) / (Pb - Pf).length() - (Pff - Pf) / (Pff - Pf).length();
				Vb = -((Pbb - Pb) / (Pbb - Pb).length() - (Pf - Pb) / (Pf - Pb).length());
			}
		}
		else if (FLG == 3) {
			//外心点法
			Vct3 Pff, Pbb;

			if (j == 0) {

				Pff = vecPt[numP - 2];
				Pbb = vecPt[j + 2];


				if (IsCloseCurve == 0) {

					Vf = Pb - Pf;

					{
						Vct3 a = Pf - Pb;
						Vct3 b = Pbb - Pb;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vb = V_rot(a ^ b, h_, -PI * 0.5);
					}
				}
				else {
					{
						Vct3 a = Pff - Pf;
						Vct3 b = Pb - Pf;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vf = V_rot(a ^ b, h_, +PI * 0.5);
					}
					{
						Vct3 a = Pf - Pb;
						Vct3 b = Pbb - Pb;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vb = V_rot(a ^ b, h_, -PI * 0.5);
					}
				}

			}
			else if (j == (numP - 2)) {

				Pff = vecPt[j - 1];
				Pbb = vecPt[1];

				if (IsCloseCurve == 0) {
					{
						Vct3 a = Pff - Pf;
						Vct3 b = Pb - Pf;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vf = V_rot(a ^ b, h_, +PI * 0.5);
					}
					Vb = -(Pb - Pf);
				}
				else {
					{
						Vct3 a = Pff - Pf;
						Vct3 b = Pb - Pf;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vf = V_rot(a ^ b, h_, +PI * 0.5);
					}
					{
						Vct3 a = Pf - Pb;
						Vct3 b = Pbb - Pb;
						Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
						Vb = V_rot(a ^ b, h_, -PI * 0.5);
					}
				}

			}
			else {

				Pff = vecPt[j - 1];
				Pbb = vecPt[j + 2];

				{
					Vct3 a = Pff - Pf;
					Vct3 b = Pb - Pf;
					Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
					Vf = V_rot(a ^ b, h_, +PI * 0.5);
				}
				{
					Vct3 a = Pf - Pb;
					Vct3 b = Pbb - Pb;
					Vct3 h_ = -(b | b) * ((a | b) - (a | a)) * a - (a | a) * ((a | b) - (b | b)) * b;
					Vb = V_rot(a ^ b, h_, -PI * 0.5);
				}

			}
		}

		Vct3 Pcnt;

		Vf = Vf / Vf.length();
		Vb = Vb / Vb.length();

		ITA_ArcConnectPoint(Pf, Vf, Pb, Vb, Pcnt);

		Vct3 Pt, Vm, Vn;
		double ang_rad, radius;
		strct_arc arc_f, arc_b;

		ITA_ArcConv(Pf, Pcnt - Pf, Vf, arc_f);
		ITA_ArcConv(Pb, Pcnt - Pb, Vb, arc_b);

		vecArc.push_back(arc_f);

		// 連続円弧にするため、２番目の円弧の主軸を反転します。
		strct_arc arc2;
		arc2.Ang_rad = arc_b.Ang_rad;
		arc2.Org = arc_b.Org;
		arc2.Radius = arc_b.Radius;
		arc2.Vm = -arc_b.Vm;  //ここ
		arc2.Vn = V_rot(arc_b.Vm, arc_b.Vn, arc_b.Ang_rad);  //arc_b.Vmの回りに、arc_b.Vnを回転する


		vecArc.push_back(arc2);


	}

	return(0);

}



int main()
{

	std::cout << "------------------------------------------------------------ " << std::endl;
	std::cout << "Welcome to ITA Spline Test Program       2026-04-01 k.mori   " << std::endl;
	std::cout << "------------------------------------------------------------ " << std::endl;
	std::cout << std::endl;

	int Kind;
	int Type;  //11  30  31

	std::cout << "Input Kind(1=Helix/ 2=Tornado/ 3=Pentagon/ 4=Torus knot) :";
	std::cin >> Kind;

	std::cout << "Input Type(10=orthocenter&open/ 11=orthocenter&close/ 30=circumcenter&open/ 31=circumcenter&close) :";
	std::cin >> Type;

	std::cout << "　Kind=" << Kind << " Type=" << Type << std::endl;

    std::vector<Vct3> vecPt;

	std::vector <strct_arc> vecArc;
 
	std::string strFileName0 = "./ITA_points.csv";
	//std::string strFileName0 = "C:/MISSILE_PATH/LogFolder/ITA_points.csv";

	FILE* fp0 = fopen(strFileName0.c_str(), "wt");  //出力ファイル

	

	if (Kind == 1) {

		// へリックス（螺旋）

		int PARTITION_NO_max100 = 5;   //3,4,5,
	
		int cnt = 0;

		for (double t = 0.0; t < 4 * (2 * PI); t = t + 2 * PI / PARTITION_NO_max100) {

			//if (cnt % 5 < 3) {
				Vct3 v1;
				v1.x = cos(t);
				v1.y = sin(t);
				v1.z = t / 4;
				//v1.z = t / 12;

				vecPt.push_back(v1);

				fprintf(fp0, "%8.3f,%8.3f,%8.3f\n", v1.x, v1.y, v1.z);
			//}
			cnt++;

		}

	}
	else if (Kind == 2) {

		// 半径線形竜巻(Tornado)

		int PARTITION_NO_max100 =  5;  // 6;   //3,4,5,
		double a = 1;		//半径初期値
		double b = 0.20;	// 0.25;
		double h = 0.25;	//ｚ方向ピッチ

		int M = 5;

		for (double t = 0.0; t < M * (2 * PI); t = t + 2 * PI / PARTITION_NO_max100) { 

			double r = a + b * t;
			
			Vct3 v1;
			v1.x = r * cos(t);
			v1.y = r * sin(t);
			v1.z = h * t;

			vecPt.push_back(v1);

			fprintf(fp0, "%8.3f,%8.3f,%8.3f\n", v1.x, v1.y, v1.z);

		}

	}
	else if (Kind == 3) {

		// 五角形(pentagon)

		double a = 0;
	 
		for (int i = 0; i < 5; i++) {

			Vct3 v1;
			double b = (i + 2) * 0.1 * PI;
			a = a + b;

			v1.set(cos(a), sin(a), 0);
			vecPt.push_back(v1);

		}
		 

		for (int j = 0; j < vecPt.size(); j++) {

			fprintf(fp0, "%8.3f,%8.3f,%8.3f\n", vecPt[j].x, vecPt[j].y, vecPt[j].z);

		}

	}
	else if (Kind == 4) {

		// トーラスノット(Torus knot)

		int PARTITION_NO_max100 = 12;   //3,4,5,

		int p = 3;
		int q = 4;

		//int p = 2;
		//int q = 3;

		int N = 1.0 / PARTITION_NO_max100;

		for (int i = 0; i < PARTITION_NO_max100; i++) {

			double t = t + 1.0 / PARTITION_NO_max100;

			Vct3 v1;

			double r = 1 + 0.5 * cos(2 * PI * p * t);
			
			double b = 2 * PI * q * t;

			double z = 0.5 * sin(2 * PI * p * t);

			v1.x = r * cos(b);
			v1.y = r * sin(b);
			v1.z = z;

			vecPt.push_back(v1);

			fprintf(fp0, "%8.3f,%8.3f,%8.3f\n", v1.x, v1.y, v1.z);

		}
	}
	else {

		// Direct Points

		Vct3 v1;
		v1.set(0, 0, 0); vecPt.push_back(v1);
		v1.set(3, 0, 0); vecPt.push_back(v1);
		v1.set(3, 4, 0); vecPt.push_back(v1);
		v1.set(3, 4, 5); vecPt.push_back(v1);
		v1.set(0, 4, 5); vecPt.push_back(v1);
		v1.set(0, 4, 0); vecPt.push_back(v1);

		for (int j = 0; j < vecPt.size(); j++) {

			fprintf(fp0, "%8.3f,%8.3f,%8.3f\n", vecPt[j].x, vecPt[j].y, vecPt[j].z);

		}

	}


	fclose(fp0);

	bool IsCloseCurve = Type % 2;

	if (IsCloseCurve == 1) {
		//終点＝始点
		Vct3 P0 = vecPt[0];
		vecPt.push_back(P0);
	}

	ITA_Points_to_Arcs(vecPt, Type, vecArc);


	std::string strFileName = "./ITA_spline.csv";
	//std::string strFileName = "C:/MISSILE_PATH/LogFolder/ITA_spline.csv";
 
	FILE* fp1 = fopen(strFileName.c_str(), "wt");  //出力ファイル

	int n1 = vecArc.size();
	printf("vecArc.size()= %5d\n", n1);


	double aL = ArcLength(vecArc);

	printf("ArcLength= %8.2f\n",aL);

	for (double len = 0; len < aL; len = len + aL / 1000) {

		Vct3 Pt;
		double Radius;

		int err = Point_in_Arcs(vecArc, len, Pt, Radius);

		if (err < 0) {
			fprintf(fp1, "ERROR!!!");
			goto L100;
		}

		fprintf(fp1, "%8.3f,%8.3f,%8.3f,%8.3f\n", Pt.x, Pt.y, Pt.z, Radius);


	}


L100:


	fclose(fp1);



	std::string strFileName2 = "./ITA_arcs.csv";
	//std::string strFileName2 = "C:/MISSILE_PATH/LogFolder/ITA_arcs.csv";

	FILE* fp2 = fopen(strFileName2.c_str(), "wt");  //出力ファイル

	for (int i = 0; i < n1; i++) {

		
		fprintf(fp2, "%3d,%8.3f,%8.3f\n",i, vecArc[i].Radius, vecArc[i].Ang_rad);
	}

	fclose(fp2);
}

 
