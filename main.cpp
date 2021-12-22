/*
 * PMC(Sakuratani)
 *
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <chrono>
#include <string>
#include <random>

#include <Vector3d.h>
#include <memory_allocate.h>
#include <vector_gauss_quadrature.h>
//#include <command.h>
#include "Date.h"
#include "Geoparameter.h"
#include "Geocoordinate.h"
#include "pmc.h"
#include "pmc_simulation.h"

/* 設定するパラメタ */
//int Day_of_Year { 190 };
//int Minute_of_Day { 21 * 60 }; /* UT */

/* 設定するパラメタ(ここまで) */

std::ofstream ofs_log;

int main(int argc, char **argv){


  double lat=0,lon=0;
  int side=0;
  std::string str_lat(argv[1]);//./main 70の70がargv[1]なのでstr_latに代入
  std::string str_side(argv[2]);/*1がアメリカ側、-1がヨーロッパ側*/
  std::string str_doy(argv[3]);
  std::string str_time(argv[4]);
  int Day_of_Year { std::stoi(str_doy) };
  int Minute_of_Day { std::stoi(str_time) * 60 }; /* UT */
  side = std::stoi(str_side);
  //str_latをdoubleに変換
  lat = std::stod(str_lat);


  /************************************************************
   * 初期化
   ************************************************************/

  /* 調査する日の設定 */
  std::string str_hour;
  /*
  if ( argc > 1 ){
    str_hour = argv[1];
    Minute_of_Day = int( std::stod(str_hour) * 60 + 0.5 );
  }
  */
  Date date(Day_of_Year, Minute_of_Day);

  AndoLab::Vector3d <double> r_geo { Rgeo, 0.0, 0.0 }; /* ひまわりの位置 */

  constexpr int NUM_ALPHA { 1 };
  double** fitted_latlon = AndoLab::allocate_memory2d(2, NUM_ALPHA, 0.0);

  lon = search_latlon(lat,side);
  std::cout << "lat=" << lat << "lon=" << lon << std::endl;


  fitted_latlon[IDX_LATITUDE][0] = lat;
//  fitted_latlon[IDX_LONGITUDE][0] = 77.1;
  fitted_latlon[IDX_LONGITUDE][0] = lon;

  /* phi は正になるので、α < 0 */
  double geo_phi { phi(lat2theta(fitted_latlon[IDX_LATITUDE][0])) };

  /* ひまわりの東経より低い、（西経もあるので、ひまわりから見て左半球（ヨーロッパ側） */
  if ( (fitted_latlon[IDX_LONGITUDE][0] < Longitude_of_himawari) &&
      (fitted_latlon[IDX_LONGITUDE][0] > 180.0 - Longitude_of_himawari) ){
    geo_phi *= -1.0;
  }

  AndoLab::Vector3d <double>
  r0(Radius_of_Earth,
      lat2theta(fitted_latlon[IDX_LATITUDE][0]), geo_phi,
      AndoLab::coordinate::Spherical);
  std::cout << "theta = " << lat2theta(fitted_latlon[0][0])
      << ", phi = " << phi(lat2theta(fitted_latlon[0][0])) << std::endl;

  Geocoordinate Geo_r0( r0 ); /* 地表面の正接点 */

  double *alpha = new double [NUM_ALPHA];
  alpha[0] = Geo_r0.alpha();
  std::cout << "Alpha = " << alpha[0] << std::endl;

  /************************************************************
   * 初期化（ここまで）
   ************************************************************/

  /* 観測データの取得 */
  double ***Observed_data
  = AndoLab::allocate_memory3d(Num_Lambda, Number_of_Obsrvd_data_Latitude, N_alt, 0.0);
  double **Obsrvd_LatLon = AndoLab::allocate_memory2d(2, Number_of_Obsrvd_data_Latitude, 0.0);
  get_observation_data(N_alt, Altitude_min, dAlt, Observed_data, Obsrvd_LatLon);

  /* 観測データのうち、フィッティングに使用する緯度経度インデックス配列 */
  int *idx_fitting_data = new int [NUM_ALPHA];
  set_fitting_latlon(idx_fitting_data, Obsrvd_LatLon, NUM_ALPHA, fitted_latlon);

  /* レイリー散乱のみの計算(こちらはフィッティングする緯度経度のみ計算) */
  double ***Rayleigh = AndoLab::allocate_memory3d(Num_Lambda, NUM_ALPHA, N_alt, 0.0);
  get_rayleigh(date, NUM_ALPHA, alpha, N_alt, Altitude_min, dAlt, Rayleigh);

  /* フィッティング(観測値を見てから) */
  double **DL = AndoLab::allocate_memory2d(Num_Lambda, NUM_ALPHA, 0.0);
  double **C = AndoLab::allocate_memory2d(Num_Lambda, NUM_ALPHA, 0.0);
  std::cout << "a" << std::endl;
  fitting_rayleigh(NUM_ALPHA, alpha,
      N_alt, Altitude_min, dAlt, Rayleigh, Observed_data, idx_fitting_data,
      DL, C, lat);

  std::ofstream ofs( ("data/C"+str_side+".dat").c_str() , std::ios::app);
  ofs << C[0][0] << " " << lat << " " << DL[0][0] << " " << Rayleigh[0][0][0]*C[0][0]+DL[0][0] << "\n";
  ofs.close();

  delete [] idx_fitting_data;
  AndoLab::deallocate_memory3d(Observed_data);
  AndoLab::deallocate_memory2d(Obsrvd_LatLon);
  AndoLab::deallocate_memory2d(fitted_latlon);
  AndoLab::deallocate_memory3d(Rayleigh);
  AndoLab::deallocate_memory2d(DL);
  AndoLab::deallocate_memory2d(C);

  return 0;
}
