#include <iostream>
#include <cmath>

#include <Vector3d.h>
#include "Date.h"
#include "pmc.h"
#include "Geoparameter.h"

AndoLab::Vector3d <double> r_s(Date d){ /* Date d はUTで与えられる */
  constexpr double inclination_of_axis = 23.43; /* [deg] 公転の面からの地軸の傾き */
  constexpr double deg2rad = M_PI / 180.0; /* [rad/deg] */
  constexpr double gamma = inclination_of_axis * deg2rad;

  constexpr double minutes_of_a_day = 24. * 60; /* 1日の分数 */
  /*
   * ひまわりの経度で南中するUT
   * 0°E では UT12時で南中、15°E では UT11時で南中、
   * 同様に東経15°ごとに UT12時より 1時間ずつ早くなる
   */
  constexpr double minute_offset = (12. - Longitude_of_himawari/15.0) * 60.0; /* 2h37m12s */

  /*
   * 指定時刻における太陽の経度の角度φs
   *  1分あたり、1/(1日の分24*60) * 2π[rad]だけ
   * 反対(-φ方向)に周る。
   *
   * Formula A (Azimuthal angle of the sun)
   */
  const double phi_s
  = -2.0*M_PI * (d.minute_day() - minute_offset) / minutes_of_a_day;

  //  std::cout << phi_s / M_PI * 180.0 << std::endl;
  /*
   * 夏至からの日数、1年で規格化する
   */
  const double s = (d.doy_from_solstice() + d.minute_day()/minutes_of_a_day) / 365.;

  /*
   * Formula B (Polar angle of the sun)
   */
  const double cos_s_sin_g = std::cos(2.0*M_PI*s) * std::sin( gamma );
  const double th = std::atan2(
      std::sqrt( 1 - cos_s_sin_g * cos_s_sin_g ), cos_s_sin_g
      );
  return AndoLab::Vector3d <double> (1.0, th, phi_s, AndoLab::coordinate::Spherical );
}

/*
 * '21 9/3: E140.8°に修正。コメントを追加
 */
