/****************************************************************************
 * Copyright (C) 2014 by Brendan Duncan.                                    *
 *                                                                          *
 * This file is part of DartRay.                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License");          *
 * you may not use this file except in compliance with the License.         *
 * You may obtain a copy of the License at                                  *
 *                                                                          *
 * http://www.apache.org/licenses/LICENSE-2.0                               *
 *                                                                          *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 *                                                                          *
 * This project is based on PBRT v2 ; see http://www.pbrt.org               *
 * pbrt2 source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.  *
 ****************************************************************************/
part of core;

double Fdr(double eta) {
  if (eta >= 1.0) {
    return -1.4399 / (eta * eta) + 0.7099 / eta + 0.6681 + 0.0636 * eta;
  } else {
    return -0.4399 + 0.7099 / eta - 0.3319 / (eta * eta) + 0.0636 / (eta * eta *
        eta);
  }
}

double _RdIntegral(double alphap, double A) {
  double sqrtTerm = Math.sqrt(3.0 * (1.0 - alphap));
  return alphap / 2.0 * (1.0 + Math.exp(-4.0 / 3.0 * A * sqrtTerm)) * Math.exp(
      -sqrtTerm);
}


double _RdToAlphap(double reflectance, double A) {
  double alphaLow = 0.0;
  double alphaHigh = 1.0;
  double kd0 = _RdIntegral(alphaLow, A);
  double kd1 = _RdIntegral(alphaHigh, A);
  for (int i = 0; i < 16; ++i) {
    //Assert(kd0 <= reflectance  kd1 >= reflectance);
    double alphaMid = (alphaLow + alphaHigh) * 0.5;
    double kd = _RdIntegral(alphaMid, A);
    if (kd < reflectance) {
      alphaLow = alphaMid;
      kd0 = kd;
    } else {
      alphaHigh = alphaMid;
      kd1 = kd;
    }
  }
  return (alphaLow + alphaHigh) * 0.5;
}


double PhaseIsotropic(Vector a, Vector b) {
  return 1.0 / (4.0 * Math.PI);
}


double PhaseRayleigh(Vector w, Vector wp) {
  double costheta = Vector.Dot(w, wp);
  return 3.0 / (16.0 * Math.PI) * (1 + costheta * costheta);
}


double PhaseMieHazy(Vector w, Vector wp) {
  double costheta = Vector.Dot(w, wp);
  return (0.5 + 4.5 * Math.pow(0.5 * (1.0 + costheta), 8.0)) / (4.0 * Math.PI);
}


double PhaseMieMurky(Vector w, Vector wp) {
  double costheta = Vector.Dot(w, wp);
  return (0.5 + 16.5 * Math.pow(0.5 * (1.0 + costheta), 32.0)) / (4.0 *
      Math.PI);
}


double PhaseHG(Vector w, Vector wp, double g) {
  double costheta = Vector.Dot(w, wp);
  return 1.0 / (4.0 * Math.PI) * (1.0 - g * g) / Math.pow(1.0 + g * g - 2.0 * g
      * costheta, 1.5);
}


double PhaseSchlick(Vector w, Vector wp, double g) {
  // improved g->k mapping derived by Thies Heidecke
  // see http://pbrt.org/bugtracker/view.php?id=102
  double alpha = 1.5;
  double k = alpha * g + (1.0 - alpha) * g * g * g;
  double kcostheta = k * Vector.Dot(w, wp);
  return 1.0 / (4.0 * Math.PI) * (1.0 - k * k) / ((1.0 - kcostheta) * (1.0 -
      kcostheta));
}


bool GetVolumeScatteringProperties(String name, Spectrum sigma_a, Spectrum
    sigma_prime_s) {
  if (!MEASURED_SS.containsKey(name)) {
    return false;
  }

  sigma_a.setRGB(MEASURED_SS[name][1][0], MEASURED_SS[name][1][1],
      MEASURED_SS[name][2][2]);

  sigma_prime_s.setRGB(MEASURED_SS[name][0][0], MEASURED_SS[name][0][1],
      MEASURED_SS[name][0][2]);
  return true;
}

void SubsurfaceFromDiffuse(Spectrum Kd, double meanPathLength, double
    eta, Spectrum sigma_a, Spectrum sigma_prime_s) {
  double A = (1.0 + Fdr(eta)) / (1.0 - Fdr(eta));
  List<double> rgb = [Kd.c[0], Kd.c[1], Kd.c[2]];

  List<double> sigma_prime_s_rgb = new List<double>(3);
  List<double> sigma_a_rgb = new List<double>(3);

  for (int i = 0; i < 3; ++i) {
    // Compute alpha' for RGB component, compute scattering properties
    double alphap = _RdToAlphap(rgb[i], A);
    double sigma_tr = 1.0 / meanPathLength;
    double sigma_prime_t = sigma_tr / Math.sqrt(3.0 * (1.0 - alphap));
    sigma_prime_s_rgb[i] = alphap * sigma_prime_t;
    sigma_a_rgb[i] = sigma_prime_t - sigma_prime_s_rgb[i];
  }

  sigma_a.setRGB(sigma_a_rgb[0], sigma_a_rgb[1], sigma_a_rgb[2]);
  sigma_prime_s.setRGB(sigma_prime_s_rgb[0], sigma_prime_s_rgb[1],
      sigma_prime_s_rgb[2]);
}

// name: [sigma_prime_s][sigma_a]
const Map<String, List<List<double>>> MEASURED_SS = const {
  // From "A Practical Model for Subsurface Light Transport"
  // Jensen, Marschner, Levoy, Hanrahan
  // Proc SIGGRAPH 2001
  'Apple': const [const [2.29, 2.39, 1.97], const [0.0030, 0.0034, 0.046]],
  'Chicken1': const [const [0.15, 0.21, 0.38], const [0.015, 0.077, 0.19]],
  'Chicken2': const [const [0.19, 0.25, 0.32], const [0.018, 0.088, 0.20]],
  'Cream': const [const [7.38, 5.47, 3.15], const [0.0002, 0.0028, 0.0163]],
  'Ketchup': const [const [0.18, 0.07, 0.03], const [0.061, 0.97, 1.45]],
  'Marble': const [const [2.19, 2.62, 3.00], const [0.0021, 0.0041, 0.0071]],
  'Potato': const [const [0.68, 0.70, 0.55], const [0.0024, 0.0090, 0.12]],
  'Skimmilk': const [const [0.70, 1.22, 1.90], const [0.0014, 0.0025, 0.0142]],
  'Skin1': const [const [0.74, 0.88, 1.01], const [0.032, 0.17, 0.48]],
  'Skin2': const [const [1.09, 1.59, 1.79], const [0.013, 0.070, 0.145]],
  'Spectralon': const [const [11.6, 20.4, 14.9], const [0.00, 0.00, 0.00]],
  'Wholemilk': const [const [2.55, 3.21, 3.77], const [0.0011, 0.0024, 0.014]],

  // From 'Acquiring Scattering Properties of Participating Media by Dilution',
  // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
  // Proc SIGGRAPH 2006
  'Lowfat Milk': const [const [0.912600, 1.074800, 1.250000],
                        const [0.000200, 0.000400, 0.000800]],
  'Reduced Milk': const [const [1.075000, 1.221300, 1.394100],
                         const [0.000200, 0.000400, 0.001000]],
  'Regular Milk': const [const [1.187400, 1.329600, 1.460200],
                         const [0.000100, 0.000300, 0.001300]],
  'Espresso': const [const [0.437600, 0.511500, 0.604800],
                     const [0.166900, 0.228700, 0.307800]],
  'Mint Mocha Coffee': const [const [0.190000, 0.260000, 0.350000],
                              const [0.098400, 0.151900, 0.204000]],
  'Lowfat Soy Milk': const [const [0.141900, 0.162500, 0.274000],
                            const [0.000100, 0.000500, 0.002500]],
  'Regular Soy Milk': const [const [0.243400, 0.271900, 0.459700],
                             const [0.000100, 0.000500, 0.003400]],
  'Lowfat Chocolate Milk': const [const [0.428200, 0.501400, 0.579100],
                                  const [0.000500, 0.001600, 0.006800]],
  'Regular Chocolate Milk': const [const [0.735900, 0.917200, 1.068800],
                                   const [0.000700, 0.003000, 0.010000]],
  'Coke': const [const [0.714300, 1.168800, 1.716900],
                 const [0.696600, 1.148000, 1.716900]],
  'Pepsi': const [const [0.643300, 0.999000, 1.442000],
                  const [0.637500, 0.984900, 1.442000]],
  'Sprite': const [const [0.129900, 0.128300, 0.139500],
                   const [0.123000, 0.119400, 0.130600]],
  'Gatorade': const [const [0.400900, 0.418500, 0.432400],
                     const [0.161700, 0.125800, 0.057900]],
  'Chardonnay': const [const [0.157700, 0.174800, 0.351200],
                       const [0.154700, 0.170100, 0.344300]],
  'White Zinfandel': const [const [0.176300, 0.237000, 0.291300],
                            const [0.173200, 0.232200, 0.284700]],
  'Merlot': const [const [0.763900, 1.642900, 1.919600],
                   const [0.758600, 1.642900, 1.919600]],
  'Budweiser Beer': const [const [0.148600, 0.321000, 0.736000],
                           const [0.144900, 0.314100, 0.728600]],
  'Coors Light Beer': const [const [0.029500, 0.066300, 0.152100],
                             const [0.026800, 0.060800, 0.152100]],
  'Clorox': const [const [0.160000, 0.250000, 0.330000],
                   const [0.017500, 0.077700, 0.137200]],
  'Apple Juice': const [const [0.121500, 0.210100, 0.440700],
                        const [0.101400, 0.185800, 0.408400]],
  'Cranberry Juice': const [const [0.270000, 0.630000, 0.830000],
                            const [0.257200, 0.614500, 0.810400]],
  'Grape Juice': const [const [0.550000, 1.250000, 1.530000],
                        const [0.542800, 1.250000, 1.530000]],
  'Ruby Grapefruit Juice': const [const [0.251300, 0.351700, 0.430500],
                                  const [0.089600, 0.191100, 0.263600]],
  'White Grapefruit Juice': const [const [0.360900, 0.380000, 0.563200],
                                   const [0.009600, 0.013100, 0.039500]],
  'Shampoo': const [const [0.028800, 0.071000, 0.095200],
                    const [0.018400, 0.059600, 0.080500]],
  'Strawberry Shampoo': const [const [0.021700, 0.078800, 0.102200],
                               const [0.018900, 0.075600, 0.098900]],
  'Head & Shoulders Shampoo': const [const [0.367400, 0.452700, 0.521100],
                                     const [0.088300, 0.163700, 0.212500]],
  'Lemon Tea': const [const [0.340000, 0.580000, 0.880000],
                      const [0.260200, 0.490200, 0.772700]],
  'Orange Juice Powder': const [const [0.337700, 0.557300, 1.012200],
                                const [0.144900, 0.344100, 0.786300]],
  'Pink Lemonade': const [const [0.240000, 0.370000, 0.450000],
                          const [0.116500, 0.236600, 0.319500]],
  'Cappuccino Powder': const [const [0.257400, 0.353600, 0.484000],
                              const [0.192000, 0.265400, 0.327200]],
  'Salt Powder': const [const [0.760000, 0.868500, 0.936300],
                        const [0.511500, 0.586300, 0.614700]],
  'Sugar Powder': const [const [0.079500, 0.175900, 0.278000],
                         const [0.065000, 0.159700, 0.257800]],
  'Suisse Mocha': const [const [0.509800, 0.647600, 0.794400],
                         const [0.187500, 0.289300, 0.379600]],
  'Pacific Ocean Surface Water': const [const [3.364500, 3.315800, 3.242800],
                                        const [3.184500, 3.132400, 3.014700]],
};
