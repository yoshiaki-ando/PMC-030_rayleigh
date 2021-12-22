/*
 * PMC
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
int Day_of_Year { 190 };
int Minute_of_Day { 21 * 60 }; /* UT */

/* 設定するパラメタ(ここまで) */

std::ofstream ofs_log;

int main(int argc, char **argv){

  /************************************************************
   * 初期化
   ************************************************************/

  /* 調査する日の設定 */
//  std::string str_hour;
//  if ( argc > 1 ){
//    str_hour = argv[1];
//    Minute_of_Day = int( std::stod(str_hour) * 60 + 0.5 );
//  }
  Date date(Day_of_Year, Minute_of_Day);

  AndoLab::Vector3d <double> r_geo { Rgeo, 0.0, 0.0 }; /* ひまわりの位置 */

  constexpr int NUM_ALPHA { 1 };
  double** fitted_latlon = AndoLab::allocate_memory2d(2, NUM_ALPHA, 0.0);
//  fitted_latlon[IDX_LATITUDE][0] = 70.0;
//  fitted_latlon[IDX_LONGITUDE][0] = 204.3;
  if ( argc < 2 ){
    std::cout << "give lat & long as arguments.\n";
    exit(0);
  }
  std::string str_lat = argv[1];
  std::string str_lon = argv[2];
  fitted_latlon[IDX_LATITUDE][0] = std::stod( str_lat );
  fitted_latlon[IDX_LONGITUDE][0] = std::stod( str_lon );

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
  std::cout << "Lat, Lon = " << Geo_r0.latitude() << " , " << Geo_r0.longitude() << std::endl;

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
  fitting_rayleigh(NUM_ALPHA, alpha,
      N_alt, Altitude_min, dAlt, Rayleigh, Observed_data, idx_fitting_data,
      DL, C);
  std::ofstream ofs_c("C.dat", std::ios::app);
  ofs_c << fitted_latlon[IDX_LATITUDE][0] << " "
      << fitted_latlon[IDX_LONGITUDE][0] << " ";
  for(int lam = 0; lam < 3; lam++){
    ofs_c << C[lam][0] << " ";
  }
  ofs_c << std::endl;
  ofs_c.close();

  delete [] idx_fitting_data;
  AndoLab::deallocate_memory3d(Observed_data);
  AndoLab::deallocate_memory2d(Obsrvd_LatLon);
  AndoLab::deallocate_memory2d(fitted_latlon);
  AndoLab::deallocate_memory3d(Rayleigh);
  AndoLab::deallocate_memory2d(DL);
  AndoLab::deallocate_memory2d(C);

  return 0;
}


