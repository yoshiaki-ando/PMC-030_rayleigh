#include <iostream>
#include <cmath>
#include <Vector3d.h>
#include "Date.h"
#include "Msis.h"
#include "pmc_simulation.h"



bool is_illuminated(AndoLab::Vector3d <double> r, Date date){

  double rQ = AndoLab::abs( r_s(date) * (r * r_s(date)) );

  return ( r%r_s(date) >= 0.0 ) ||
      ( (r%r_s(date) < 0.0) && (rQ >= Radius_of_Earth) );
}

bool is_in_shadow_zone(AndoLab::Vector3d <double> r, Date date){
  return !is_illuminated(r,date);
}

/* r の点から太陽方向へ移動して大気圏界面に達する点 */
AndoLab::Vector3d <double> upper_atmosphere_to_solar(
    AndoLab::Vector3d <double> r, Date date){

  return cross_point_at_altitude(r, r_s(date).n(),
      Altitude_of_Atmosphere );
}

bool is_inside_PMC_region(AndoLab::Vector3d <double> r){
  return ( (r.abs() > R_lower_pmc) && (r.abs() < R_upper_pmc) );
}


double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const double lambda,
    const double th, /* 散乱角 */
    Msis &msis,
    double &delta,
    const int num_pmc,
    PMC *pmc,
    const double Interval){

  double I { 0.0 };

  /* 寄与を計算するベクトル */
  AndoLab::Vector3d <double> length = r_to - r_from;

  /* 分割点数。Interval毎 ===>> 要調査 */
  const int Num = int( length.abs() / Interval + 1 );
  AndoLab::Vector3d <double> dr = length / double(Num); /* 分割距離 */

  AndoLab::Vector3d <double> pre_rp = r_from;

  for(int i_dr = 0; i_dr < Num; i_dr++){
    AndoLab::Vector3d <double> r_p = r_from + (i_dr+0.5)*dr;

    /* 視線方向の光学的深さは、視線から遠い方向へと計算してゆけば
     * 先の計算が使える */
    delta += Optical_depth(lambda, date, pre_rp, r_p, num_pmc, pmc);

    if ( is_illuminated(r_p, date) ){
      msis.set_position( Geocoordinate(r_p) );

      /* r の点がシャドウ領域でなければ
       * 1. 太陽方向の大気圏界面の点 r0 を見つける
       * 2. exp( - delta(r0,r) - delta(r,r2geo) ) * σ(θ) */
      AndoLab::Vector3d <double> r_i =
          upper_atmosphere_to_solar(r_p, date);

      /* 散乱係数 */
      double total_beta = sigma(th, lambda, msis) * msis.N();
      if ( is_inside_PMC_region( r_p ) ){
        for(int npmc = 0; npmc < num_pmc; npmc++){
          total_beta += pmc[npmc].bn_th() * pmc[npmc].dist( r_p );
        }
      }

      I += std::exp(
          -delta - Optical_depth(lambda, date, r_i, r_p, num_pmc, pmc) ) *
          total_beta * dr.abs();
    }
    pre_rp = r_p;
  }
  delta += Optical_depth(lambda, date, pre_rp, r_to, num_pmc, pmc );

  return I;
}
