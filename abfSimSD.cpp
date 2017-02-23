/*
This program predicts the APD voltage reading and uncertainty
of fluorescence from Yb or Mg that is produced at the viewport 
in the SADiCS experiment over a range of temperatures and 
records the data in a file

*/


#include<iostream>
using std::cin; using std::cout; using std::endl;
#include<string>
using std::string;
#include<algorithm>
using std::transform;
#include<cctype>
#include<cmath>
#include<fstream>
using std::ofstream;
#include<math.h>



int main() {

  const double Pi = 3.1415926535897;
  //--- BoltzConst = Boltzmann's constant in cm^2*kg*s^(-2)K^(-1)
  const double BoltzConst = 1.38064852e-19;
  //--- SpeedLight = speed of light in cm/s
  const double SpeedLight = 2.9979e10;
  // --- H = planck's constant in cm^2*kg/s
  const double H = 6.62607004e-30;
  //--- ClasElecRad = classical electron radius in cm
  const double ClasElecRad = 2.8179403267e-13;

  //--- r = radius of source slit, l = length of source canal, d = distance from the oven to the observation point, atomic_rad = radius of atomic beam at viewport (all in cm
  double r = 0.1575, l = 1.27, d = 12.751, atomic_rad = 2*r*d/l;


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

  string filename = "";
  
  cout<<"What file would you like to write to?"<<endl;
  cin>>filename;

  ofstream write_file(filename);

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

    //---
  //--- Laser flux
  //---
  valid = true;

  cout<<"What power will you operate your laser at? (Max = 3.5W)"<<endl;
  double p = 0;
  while (valid == true) {
    cin>>p;
    if (p <= 3.5 && p > 0) {
      valid = false;
    }
    else {
      cout<<"Power not valid. Please enter a power between 0 and 3.5. "<<endl;
    }
  }
  // --- convert power to kg cm^2 per s^3
  p *= 1e4;


  //--- w_a = semi-major radius of laser, w_b = semi-minor...
  double w_a = 0, w_b = 0;

  cout<<"What was the semi-major and semi-minor radius of the laser spot in micrometers?"<<endl;
  cin>>w_a;
  cin>>w_b;

  //--- convert radii to cm
  w_a /= 1e4;
  w_b /= 1e4;

  for (double i = 0; i<= 350; i++) {

    //--- Convert temperature to Kelvin to agree with our equations
    double temp = i + 273.15;
  
    double flux = 0;

    //--- vp = vapor pressure & conversion to kg / (cm*s^2)
    double vp = pow(10.0,(a + (b/temp) + (c*log10(temp))));
    vp /= 100;

    //--- n = number density, v = rms velocity
    double n = vp/(BoltzConst*temp);
    double v = pow((3.0*BoltzConst*temp)/mass,0.5);
    double one_kappa = 8.0/3*(r/l);


    //cout<<"vapor pressure: "<<vp<<endl;
    //cout<<"number density: "<<n<<endl;
    //cout<<"rms velocity: "<<v<<endl;

    flux = 1.0/4*(one_kappa)*pow((r/d),2.0)*n*v;
  
    //
    //--- Find the doppler shifted wavelength
    //

    //--- v_mp = most probable velocity
    double v_mp = pow((2.0*BoltzConst*temp)/mass,0.5);

    // cout<<"Most probable velocity: "<<v_mp<<endl;
  
    //--- f = frequency, f_emit = shifted frequency
    double f = SpeedLight/lambda;
    double f_emit = f/(1+v_mp/SpeedLight); 
    double lambda_emit = SpeedLight/f_emit;


    double laser_flux = (2.0*p)/(Pi*w_a*w_b*H*f_emit);

    //cout<<"The laser has an initial output of "<<laser_flux<<" photons per cm^2*s"<<endl;

    //
    //--- cross section
    //  

    //--- Full width at half maximum and lorentzian
    double FWHM = f*sqrt((8*BoltzConst*temp*log(2))/(mass*pow(SpeedLight,2.0)));
    double lorentz = (FWHM/(2.0*Pi))/(pow((f_emit-f),2.0)+pow(FWHM/2.0,2.0));

    double cross_section = Pi*ClasElecRad*SpeedLight*osc_str*lorentz;

    //cout<<"Cross section: "<<cross_section<<" cm^2"<<endl;

    //
    //--- net rate
    //
  
    double net_rate = laser_flux*cross_section;

    //cout<<"Net rate: "<<net_rate<<" photons per s"<<endl;

    //
    //--- fiducial volume
    //

    double rad = pow(pow(w_a,2.0)+pow(w_b,2.0),1/2);
    double num_fid = 4*flux*pow(rad,2.0)*pow(atomic_rad,2.0)*(Pi/2*rad*sin(Pi*rad/(2*atomic_rad))+2*atomic_rad*cos(Pi*rad/(2*atomic_rad))-2*atomic_rad);


    // cout<<"Fiducial vol: "<<fid_vol<<" cm^3"<<endl;
    // cout<<"n3: "<<n3<<" particles per cm^3"<<endl;
    // cout<<"Number of atoms in the fiducial volume: "<<num_fid<<endl;

    //
    //--- light signal
    //


    //--- solid angle of (1-cos(th))/2, th =tan^(-1)(.25mm/135.184mm) deg
    double solid_ang = 8.55004355e-7;
    double sig =  solid_ang*num_fid*net_rate;
    double pow_sig = sig*H*f_emit;

    pow_sig /= 1e4;

    double APD_volt = pow_sig*11.5/50*50*500000;

    write_file << i << "\t" << pow_sig/(p/1e4) << "\n"; 

  }
  
  return 0;
}
