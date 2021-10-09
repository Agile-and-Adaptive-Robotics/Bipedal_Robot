/* sketch adapted by Lindie Burgess for use by the Agile and Adaptive Robotics Lab at Portland State  
University.

This sketch is used in conjunction with Matlab function called, "KneeTest.m"

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


void setup() {
Serial.begin(115200);  // initialize arduino serial communication
}

void loop() 
{
  if (Serial.available()) // if information is sent over serial from matlab
  {  char choose_branch = Serial.read(); // read the serial data into variable "choose_branch"
      if (choose_branch == '1') // If choose_branch is equal to '1', iterate through the following for loop
      {  scale.begin(LOADCELL_DOUT_PIN, LOADCELL_SCK_PIN); // initialize load cell 
         scale.set_scale(); //set the scale value to 1; 
         scale.tare(); // sets OFFSET to 0
            for (int i = 0; i <150; i++)
            {
              long zero_factor = scale.read_average(); //read a baseline value into variable zero_factor
              Serial.println(zero_factor);   //print the variable zero_factor to serial
            }  
      }       

  }         
  else  //if there is no information to read over serial from matlab, wait...
  {  while(Serial.available() <1);
  }
}
