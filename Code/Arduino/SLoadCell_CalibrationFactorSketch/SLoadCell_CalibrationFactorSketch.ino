
/* sketch adapted by Lindie Burgess for use by the Agile and Adaptive Robotics Lab.

  This sketch is NOT meant to be used with Matlab.  Instead, use Arduino's built in serial monitor.

  This arduino sketch is adapted from the "SparkFun_HX711_Calibration" example, which
  can be found here
  https://github.com/sparkfun/HX711-Load-Cell-Amplifier/tree/master/firmware

  Source Code for the HX711 can be found here:
  https://github.com/bogde/HX711/tree/master/src

*/


#include "HX711.h" //This library can be obtained here http://librarymanager/All#Avia_HX711

#define LOADCELL_DOUT_PIN 3 //define the Serial Data Output Pin
#define LOADCELL_SCK_PIN 2  // define the Power Down and Serial Clock Input Pin

HX711 scale;

float calibration_factor = -17550; // initialize variable calibration_factor to some guess value


void setup() {

  Serial.begin(9600);  // initialize arduino serial communication
  scale.begin(LOADCELL_DOUT_PIN, LOADCELL_SCK_PIN); // initialize load cell
  scale.set_scale(); //set the scale value;
  //this value is used to convert the raw
  //data to "human readable" data (measure units)
  scale.set_offset(-14462);
  //This sets the offset value to a known zero.
  //There is no need for taring the scale once the zero point is known for a scale in a
  // set configuration.
  // The Zero Factor -17200 was found using a separate Arduino Sketch
  //called "SLoadCell_ZeroFactorSketch"
}

void loop() {
  //Protocol to read calibration_factor
  scale.set_scale(calibration_factor); //Adjust to this calibration factor

  Serial.print("Reading: ");
  Serial.print(scale.get_units(), 1);
  Serial.print(" lbs"); //Change this to kg and re-adjust the calibration factor if you follow SI units like a sane person
  Serial.print(" calibration_factor: ");
  Serial.print(calibration_factor);
  Serial.println();

  if (Serial.available())
    // if the following characters are sent over serial, the calibration factor
    // will adjust as commanded
  {
    char temp = Serial.read();
    if (temp == '+' || temp == 'a')
      calibration_factor += 10;
    if (temp == '+' || temp == 'a')
      calibration_factor += 10;
    else if (temp == '-' || temp == 's')
      calibration_factor -= 10;
  }
}
