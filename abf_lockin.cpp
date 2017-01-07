/*
This program takes temperature vs 
fluorescence power at APD vs laser power 
data and transforms it into temperature 
vs flux at the SADiCS viewport data

12/08/2016 Kristen approved, lastest version

1/6/2016 Kristen editted, SI units now besides printing flux as cm^-2*s^-1
*/


#include<iostream>
using std::cin; using std::cout; using std::endl;
#include<string>
using std::string;
using std::stod;
#include<algorithm>
using std::transform;
#include<cctype>
#include<cmath>
#include<fstream>
using std::ofstream;
using std::ifstream;



int main() {

  //-----CONSTANTS

  const double kPi = 3.1415926535897;
  //--- BoltzConst = Boltzmann's constant in m^2*kg*s^(-2)K^(-1)
  const double kBoltzConst = 1.38064852e-23;
  //--- SpeedLight = speed of light in m/s
  const double kSpeedLight = 2.9979e8;
  // --- Planck = planck's constant in m^2*kg/s
  const double kPlanck = 6.62607004e-34;
  //--- ClasElecRad = classical electron radius in m
  const double kClasElecRad = 2.8179403267e-15;

  //-----END CONSTANTS

 
  string writeFileName = "";
  
  // cout<<"What file would you like to read from?"<<endl;
  //cin>>readFileName;

  cout<<"What file would you like to write to?"<<endl;
  cin>>writeFileName;

  ofstream write_file(writeFileName);

  //--- nozzle_radius = radius of nozzle slit, nozzle_length = length of nozzle, nozzle2viewport = distance from the oven to the observation point, atomic_radius = radius of atomic beam at viewport (all in cm)
  double nozzle_radius = .001575, nozzle_length = .0127, nozzle2viewport = .12751, atomic_radius = 2*nozzle_radius*nozzle2viewport/nozzle_length;

  //--- w_a = semi-major radius of laser in m, w_b = semi-minor radius of laser in m...
  double w_a = .002220, w_b = .0013425;

  //--- solid angle of (1-cos(th))/2, th = tan^(-1)(.25mm/135.184mm) deg
  double solid_angle = 8.55004355e-7;

  cout<<endl;
  cout<<endl;

  /* The following arrays contain [0] mass in kg, [1] a, [2] b, [3] c, [4] excitation wavelength in m, [5] oscillator strength, where a, b and c correspond to the equation for vapor pressure*/

  //recalculate mg wavelength
  double yb[8] = {2.87339675e-25, 14.117, -8111.0, -1.0849, 398.9114186e-9, 1.37, 555.647e-9, .016};
  double mg[8] = {4.0359e-26, 13.495, -7813, -0.8253, 285.213e-9, 1.80, 457.11e-9, 2.38e-6};

  string source;
  double mass = 0;
  double a = 0, b = 0, c = 0;
  double lambda = 0;
  double oscillator_strength = 0;

  bool valid = true;

  cout<<"Will you be using Yb or Mg?"<<endl;
  while (valid == true) {

    cin>>source;
    transform(source.begin(), source.end(), source.begin(), tolower);
  
    //--- Assign source characteristics

    if (source == "yb") {
      mass = yb[0];
      a = yb[1];
      b = yb[2];
      c = yb[3];
      lambda = yb[4];
      oscillator_strength = yb[5];
      valid = false;
    }
    else if (source == "mg") {
      mass = mg[0];
      a = mg[1];
      b = mg[2];
      c = mg[3];
      lambda = mg[4];
      oscillator_strength = mg[5];
      valid = false;
    }
    else {
      cout<<"Source not valid, are you using Yb or Mg?"<<endl;
    }
  }
  
  double temp = 0;
  double apd_power = 0;
  double laser_power = 0;
  // double apd_power_sd = 0;

  

  string readFileName = "/home/kristen/SADiCS/C++/2016-12-05-183641-abf-sample-set-data-analysis.txt";
  //cin>>readFileName;
  ifstream read_file(readFileName);
  while(read_file >> temp >> apd_power >> laser_power) {
    //-- power in kg*m^2/s^3

    //--- temp to K
    temp += 273.15;

    //--- mp_velocity = most probable velocity
    double mp_velocity = pow((2.0*kBoltzConst*temp)/mass,0.5);
    double velocity = pow((3.0*kBoltzConst*temp)/mass,0.5);

    //--- vapor_pressure = vapor pressure kg/(m*s^2)
    double vapor_pressure = pow(10.0,(a + (b/temp) + (c*log10(temp))));\

    double num_density = vapor_pressure/(kBoltzConst*temp);

    //--- f = frequency, f_emit = shifted frequency
    double freq = kSpeedLight/lambda;
    double freq_emit = freq/(1+mp_velocity/kSpeedLight); 
    double lambda_emit = kSpeedLight/freq_emit;

    double laser_flux = (2.0*laser_power)/(kPi*w_a*w_b*kPlanck*freq_emit);

    double fullwidth_halfmax = freq*sqrt((8*kBoltzConst*temp*log(2))/(mass*pow(kSpeedLight,2.0)));
    double lorentzian = (fullwidth_halfmax/(2.0*kPi))/(pow((freq_emit-freq),2.0)+pow(fullwidth_halfmax/2.0,2.0));

    double cross_section = kPi*kClasElecRad*kSpeedLight*oscillator_strength*lorentzian;
    double net_rate = 1.92e8/2;//laser_flux*cross_section;

    double apd_photon_count = apd_power/(kPlanck*freq);
    //double apd_photon_count_sd = apd_power_sd/(kPlanck*f);

    double num_fiducial = apd_photon_count/(solid_angle*net_rate*.4);
    //double num_fid_sd = apd_photon_count_sd/(solid_angle*net_rate)
;
    double fiducial_volume = kPi*w_a*w_b*1/(cross_section*num_density);    
    
    double flux = num_fiducial*velocity/fiducial_volume;
    //double flux_sd = num_fid_sd*v/fid_vol;
    
    temp -= 273.15;

    write_file << temp << "\t" << flux/1e4 << "\n";
  }
  read_file.close();
  write_file.close();

  return 0;
}
