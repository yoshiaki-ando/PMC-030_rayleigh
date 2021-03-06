#include <iostream>
#include <fstream>
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

constexpr double Reflection_Coefficient { 0.3 };
constexpr int M_alpha { 30 };
constexpr int N_beta { 30 };
constexpr double Dbeta { 2.0*M_PI / N_beta };
constexpr double R0 { Radius_of_Earth };

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

  /* 反射のための各種設定 */
  AndoLab::Vector3d <double> yv( 0.0, 1.0, 0.0 );
  AndoLab::Vector3d <double> zv( 0.0, 0.0, 1.0 );
  AndoLab::Vector3d <double> r_gs = (r_from - r_to).n();

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

      /* 反射の検討 */
      double reflection = 0.0; /* 単一分子当たりの反射光 */
      const double alpha_max = std::acos( R0 / r_p.r() );
      const double Dalpha = alpha_max / M_alpha;
      const double theta_s = r_p.theta();
      const double phi_s = r_p.phi();

//      for(int m = 0; m < M_alpha; m++){
//        const double alpha = (m + 0.5) * Dalpha;
//        const double Smn = R0*R0 * std::sin(alpha) * Dalpha * Dbeta;
//
//        for(int n = 0; n < N_beta; n++){
//          AndoLab::Vector3d <double> rp_R_mn(R0, alpha, (n+0.5)*Dbeta, AndoLab::coordinate::Spherical);
//          AndoLab::Vector3d <double> r_R_mn = ( rp_R_mn.rotate(theta_s, yv) ).rotate(phi_s, zv); /* 反射点 */
//
//          /* 反射点への太陽光入射判定 */
//          if ( r_R_mn%r_s(date) <= 0.0 ){
//            continue;
//          }
//
//          AndoLab::Vector3d <double> r_Rs_mn = r_p - r_R_mn; /* 反射点→散乱点ベクトル */
//          AndoLab::Vector3d <double> r_Ri_mn =
//              upper_atmosphere_to_solar(r_R_mn, date); /* 反射点に入射する光の大気圏入射点 */
//
//          /* 反射点到達光 */
//          double I_r_mn = std::exp( - Optical_depth(lambda, date, r_Ri_mn, r_R_mn, num_pmc, pmc) );
//          double Ip_S_mn = Smn / M_PI / r_Rs_mn.abs() / r_Rs_mn.abs()
//              * I_r_mn * r_R_mn.n()%r_p.n() * r_Rs_mn.n()%r_R_mn.n(); /* 無損失散乱光 */
//          /* 損失散乱光 */
//          double I_S_mn = Ip_S_mn * std::exp( - Optical_depth(lambda, date, r_p, r_R_mn, num_pmc, pmc) );
//          double theta_s_mn = AndoLab::angle_between( r_Rs_mn, r_gs  ); /* 反射光散乱角 */
//          reflection += I_S_mn * sigma(theta_s_mn, lambda, msis);
//          reflection = 0.0;
//
//        }
//      }

      /* 反射の検討(ここまで) */

      I += std::exp(-delta) *
          (std::exp( - Optical_depth(lambda, date, r_i, r_p, num_pmc, pmc) ) * total_beta
              + reflection * msis.N() )
           * dr.abs();
    }
    pre_rp = r_p;
  }
  delta += Optical_depth(lambda, date, pre_rp, r_to, num_pmc, pmc );

  return I;
}
