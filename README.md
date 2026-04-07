# ITA Spline (Isosceles Triangle Arc Spline)

A 3D arc-based interpolation method developed in the 1990s and 
publicly released in 2026.  
It connects N points in 3D space using 2(N-1) circular arcs with C1 continuity.

## Features
- Pure arc-based interpolation
- Analytical solution using elementary geometry
- Two direction-vector models (Opposite-edge / Circumcenter)
- Reproduces a perfect circle when points lie on a circle
- C++ implementation included

## Documents
- ITA_Spline_JP_rev3.pdf

## Source Code
IDE Microsoft VC++2022

## Author
Kazuya Mori  
Original development: 1990s  
Public release: 2026
====================================================================
(Japanese)
コンピュータグラフィック分野で任意曲線を扱う際、曲線を定義する点列を補完する手法の
代表的な例に、ラグランジュ補間、Bスプライン補完などがあり、CADシステム等で古くか
ら利用されている。本稿では、３次元空間内のN個の点列に対して、２（N―１）個の曲率
半径の異なる円弧で滑らかに（C1級）つなぎ合わせる手法（ITA Spline）を紹介する。構成
がシンプルな割には、外見上の不自然さがなく、曲線の⾧さの計算も簡単であり、点列が正
多角形の頂点の場合、円を誤差なく表現できることに特徴がある。本手法は、目的、計算手
法、結果がいずれもシンプルであり、理解も容易である。既存のBiarc法との関係も説明する。 
