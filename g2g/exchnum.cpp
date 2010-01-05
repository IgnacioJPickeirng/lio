#include <iostream>
#include "common.h"
#include "init.h"
#include "timer.h"
#include "cuda/cuda_extra.h"
using namespace G2G;
using namespace std;


template<bool compute_exc, bool compute_y2a> void cpu_pot(float dens, float& ex, float& ec, float& y2a);

template<bool compute_energy, bool do_forces>
void cpu_compute_density_forces(float* energy, float* point_weights, uint points, float* rdm, float* rmm_output,
  float* function_values, float4* gradient_values, float4* forces, uint* nuc, uint nucleii_count, uint m, Timer& t, Timer& trmm)
{
  if (!compute_energy) {
    trmm.start();
     for (uint i = 0; i < m; i++) {
       for (uint j = i; j < m; j++) {
         rmm_output[COALESCED_DIMENSION(m) * j + i] = 0.0f;
       }
     }
    trmm.pause();
  }

  for (uint point = 0; point < points; point++) {
    t.start();
    //cout << "punto" << endl;
    float partial_density = 0.0f;

    HostMatrixFloat w(fortran_vars.nco, 1);
    for (uint i = 0; i < fortran_vars.nco; i++) w.get(i) = 0.0f;

    for (uint j = 0; j < m; j++) {
      float f = function_values[m * point + j];
      for (uint i = 0; i < fortran_vars.nco; i++) {
        float r = rdm[COALESCED_DIMENSION(fortran_vars.nco) * j + i];
        w.get(i) += f * r;
      }	// TODO: usar rmmt
    }

    for (uint i = 0; i < fortran_vars.nco; i++) { partial_density += w.get(i) * w.get(i); }
    partial_density *= 2;

    /* compute energy / functional */
    float exc, corr, y2a;
    float point_weight = point_weights[point];
    float factor = 0;

    if (compute_energy) {
      cpu_pot<true, false>(partial_density, exc, corr, y2a);
      energy[point] = (partial_density * point_weight) * (exc + corr);
      t.pause();
    }
    else {
      cpu_pot<false, true>(partial_density, exc, corr, y2a);
      factor = point_weight * y2a;
      t.pause();

      trmm.start();
      for (uint i = 0; i < m; i++) {
        for (uint j = i; j < m; j++) {
          float Fi = function_values[m * point + i];
          float Fj = (i == j ? Fi : function_values[m * point + j]);
          rmm_output[COALESCED_DIMENSION(m) * j + i] += Fi * Fj * factor;
        }
      }
      trmm.pause();
    }
  }
}

#define POT_ALPHA 		-0.738558766382022447 // -(3/PI)^(1/3)
#define POT_GL 				0.620350490899400087

#define POT_VOSKO_A1 	0.03109205
#define POT_VOSKO_B1 	3.72744
#define POT_VOSKO_C1 	12.9352
#define POT_VOSKO_X0 	-0.10498

#define POT_VOSKO_Q 	6.15199066246304849
#define POT_VOSKO_A16 0.005182008333
#define POT_VOSKO_A2 	0.015546025
#define POT_VOSKO_B2 	7.06042
#define POT_VOSKO_C2 	18.0578
#define POT_VOSKO_X02 -0.32500
#define POT_VOSKO_Q2 	4.7309269
#define POT_VOSKO_A26 0.0025910042

#define POT_XX0 12.5549141492 // POT_VOSKO_X0 * POT_VOSKO_X0 + POT_VOSKO_B1 * POT_VOSKO_X0 + POT_VOSKO_C1
#define POT_T6 -0.00836166609762834 // POT_VOSKO_X0 / POT_XX0
#define POT_T4 -0.0311676086789438 // POT_VOSKO_B1 * POT_VOSKO_X0 / POT_XX0
#define POT_VOSKO_2C1 25.8704 // 2 * POT_VOSKO_C1

#define POT_VOSKO_2B1Q 1.21178337371132 // 2 * POT_VOSKO_B1 / POT_VOSKO_Q
#define POT_VOSKO_B2X0Q 1.14352579286644 // 2 * (POT_VOSKO_B1 + 2 * POT_VOSKO_X0) / POT_VOSKO_Q
#define POT_VOSKO_4B1 14.90976 // 4.0 * POT_VOSKO_B1
#define POT_VOSKO_QSQ 37.8469891110325 // POT_VOSKO_Q * POT_VOSKO_Q
#define POT_VOSKO_B1X0 1.0329232240928 // (1.0f - t6 * (POT_VOSKO_B1 - 2.0f * POT_VOSKO_X0))

template<bool compute_exc, bool compute_y2a> void cpu_pot(float dens, float& ex, float& ec, float& y2a)
{
	// data X alpha

	if (dens == 0) {
		if (compute_exc) { ex = 0.0f; ec = 0.0f; }
		if (compute_y2a) y2a = 0.0f;
		return;
	}

	float y = powf(dens, 0.333333333333333333f);  // rho^(1/3)
	float v0 = -0.984745021842697f * y; // -4/3 * (3/PI)^(1/3) * rho^(1/3)

	if (compute_exc) ex = POT_ALPHA * y; // -(3/PI)^(1/3) * rho^(1/3)

	switch(fortran_vars.iexch) {
		case 1:
		{
			if (compute_exc) ec = 0;
			if (compute_y2a) y2a = v0;
		}
		break;
		case 2:
		{
			float rs = POT_GL / y;
			float x1 = rs / 11.4f;
			float vc;

			if (x1 > 1.0f) {
				ec = -0.0333f * (0.5f * x1 - 0.33333333333333f);
				if (compute_y2a) vc = 0.0111f * x1 * 0.5f;
			}
			else {
				float t1 = (1.0f + x1 * x1 * x1);
				float t2 = logf(1.0f + 1.0f / x1);
				float t3 = x1 * x1;
        ec = -0.0333f * (t1 * t2 - t3 + 0.5f * x1 - 0.33333333333333f);
        if (compute_y2a) vc = 0.0111f * x1 * (3.0f * t3 * t2 - t1 / (x1 * (x1 + 1.0f)) - 2.0f * x1 + 0.5f);
			}
			if (compute_y2a) y2a = v0 + ec + vc;
		}
		break;
		case 3:
		{
			float rs = POT_GL / y;
			float x1 = sqrtf(rs);
			float Xx = rs + POT_VOSKO_B1 * x1 + POT_VOSKO_C1;
			float t1 = 2 * x1 + POT_VOSKO_B1;
			float t2 = logf(Xx);
			float t3 = atanf(POT_VOSKO_Q/t1);
      float t5 = (POT_VOSKO_B1 * x1 + POT_VOSKO_2C1) / x1;

      ec = POT_VOSKO_A1 * (2 * logf(x1) - t2 + POT_VOSKO_2B1Q * t3 - POT_T4 * (2 * logf(x1 - POT_VOSKO_X0) - t2 + POT_VOSKO_B2X0Q * t3));

			float vc;
      if (compute_y2a) {
				vc = ec - POT_VOSKO_A16 * x1 * (t5 / Xx - POT_VOSKO_4B1 / (t1 * t1 + POT_VOSKO_QSQ) * POT_VOSKO_B1X0 - POT_T4 * (2.0f / (x1 - POT_VOSKO_X0) - t1 / Xx));
				y2a = v0 + vc;
			}
		}
		break;
	}
}
