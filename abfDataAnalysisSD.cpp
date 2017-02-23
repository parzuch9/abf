/*
This program takes temperature vs 
fluorescence power at APD vs laser power 
data and transforms it into temperature 
vs flux at the SADiCS viewport data as
well as uncertainty in flux

outputs flux as kg*cm^-2*s^-1
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

  
  //  string readFileName = "2016-10-06-061904-abf-sample-set-data-analysis.txt";
  string writeFileName = "";
  
  // cout<<"What file would you like to read from?"<<endl;
  //cin>>readFileName;

  cout<<"What file would you like to write to?"<<endl;
  cin>>writeFileName;

  ofstream write_file(writeFileName);

  
  const double Pi = 3.1415926535897;
  //--- BoltzConst = Boltzmann's constant in cm^2*kg*s^(-2)K^(-1)
  const double BoltzConst = 1.38064852e-19;
  //--- SpeedLight = speed of light in cm/s
  const double SpeedLight = 2.9979e10;
  // --- H = planck's constant in cm^2*kg/s
  const double H = 6.62607004e-30;
  //--- ClasElecRad = classical electron radius in cm
  const double ClasElecRad = 2.8179403267e-13;

  //--- r = radius of source slit, l = length of source canal, d = distance from the oven to the observation point, atomic_rad = radius of atomic beam at viewport (all in cm)
  double r = 0.1575, l = 1.27, d = 12.751, atomic_rad = 2*r*d/l;

  //--- w_a = semi-major radius of laser, w_b = semi-minor...
  double w_a = .11805, w_b = .09605;

  //--- solid angle of (1-cos(th))/2, th = tan^(-1)(.25mm/135.184mm) deg
  double solid_ang = 8.55004355e-7;

  cout<<endl;
  cout<<endl;

  /* The following arrays contain [0] mass, [1] a, [2] b, [3] c, [4] excitation wavelength, [5] oscillator strength, where a, b and c correspond to the equation for vapor pressure*/

  //recalculate mg wavelength
  double yb[8] = {2.87339675e-25, 14.117, -8111.0, -1.0849, 398.9114186e-7, 1.37, 555.647e-7, .016};
  double mg[8] = {4.0359e-26, 13.495, -7813, -0.8253, 285.213e-7, 1.80, 457.11e-7, 2.38e-6};

  string source;
  double mass = 0;
  double a = 0, b = 0, c = 0;
  double lambda = 0;
  double osc_str = 0;

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
      osc_str = yb[5];
      valid = false;
    }
    else if (source == "mg") {
      mass = mg[0];
      a = mg[1];
      b = mg[2];
      c = mg[3];
      lambda = mg[4];
      osc_str = mg[5];
      valid = false;
    }
    else {
      cout<<"Source not valid, are you using Yb or Mg?"<<endl;
    }
  }
  
  double temp = 0;
  double apd_power = 0;
  double p = 0;
  double apd_power_sd=0;

  string readFileName = "/home/kristen/SADiCS/2016-10-06-061904-abf-sample-set-data-analysis.txt";
  ifstream read_file(readFileName);
  while(read_file >> temp >> apd_power >> p >> apd_power_sd) {
    
    //--- convert power to kg cm^2 per s^3
    apd_power *= 1e4;
    p *= 1e4;

    //--- temp to K
    temp += 273.15;

    //--- v_mp = most probable velocity
    double v_mp = pow((2.0*BoltzConst*temp)/mass,0.5);
    double v = pow((3.0*BoltzConst*temp)/mass,0.5);

    //--- vp = vapor pressure & conversion to kg / (cm*s^2)
    double vp = pow(10.0,(a + (b/temp) + (c*log10(temp))));
    vp /= 100;

    double n = vp/(BoltzConst*temp);

    //--- f = frequency, f_emit = shifted frequency
    double f = SpeedLight/lambda;
    double f_emit = f/(1+v_mp/SpeedLight); 
    double lambda_emit = SpeedLight/f_emit;

    double laser_flux = (2.0*p)/(Pi*w_a*w_b*H*f_emit);

    double FWHM = f*sqrt((8*BoltzConst*temp*log(2))/(mass*pow(SpeedLight,2.0)));
    double lorentz = (FWHM/(2.0*Pi))/(pow((f_emit-f),2.0)+pow(FWHM/2.0,2.0));

    double cross_section = Pi*ClasElecRad*SpeedLight*osc_str*lorentz;
    double net_rate = laser_flux*cross_section;

    double sig = apd_power/(H*f);
    double sig_sd = apd_power_sd/(H*f);

    double num_fid = sig/(solid_ang*net_rate);
    double num_fid_sd = sig_sd/(solid_ang*net_rate);
    double fid_vol = Pi*w_a*w_b*1/(cross_section*n);
    
    double flux = num_fid*v/fid_vol;
    double flux_sd = num_fid_sd*v/fid_vol;
    
    temp -= 273.15;

    write_file << temp << "\t" << flux << "\t" << flux_sd << "\n";
  }
  read_file.close();
  write_file.close();

  return 0;
}
