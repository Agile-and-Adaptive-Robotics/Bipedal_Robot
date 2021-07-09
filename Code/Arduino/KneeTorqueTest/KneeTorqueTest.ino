/* sketch adapted by Lindie Burgess for use by the Agile and Adaptive Robotics Lab at Portland State
  University.

  This sketch is used in conjunction with Matlab function called, "KneeTest.m"

  This arduino sketch is adapted from the "SparkFun_HX711_KnownZeroStartup" example, which
  can be found here
  https://github.com/sparkfun/HX711-Load-Cell-Amplifier/tree/master/firmware

  Source Code for the HX711 can be found here:
  https://github.com/bogde/HX711/tree/master/src


  Example using the SparkFun HX711 breakout board with a scale
  By: Nathan Seidle
  SparkFun Electronics
  Date: November 19th, 2014
  License: This code is public domain but you buy me a beer if you use this and we meet someday (Beerware license).

  Most scales require that there be no weight on the scale during power on. This sketch shows how to pre-load tare values
  so that you don't have to clear the scale between power cycles. This is good if you have something on the scale
  all the time and need to reset the Arduino and not need to tare the scale.

  This example code uses bogde's excellent library: https://github.com/bogde/HX711
  bogde's library is released under a GNU GENERAL PUBLIC LICENSE

  The HX711 does one thing well: read load cells. The breakout board is compatible with any wheat-stone bridge
  based load cell which should allow a user to measure everything from a few grams to tens of tons.
  Arduino pin 2 -> HX711 CLK
  3 -> DOUT
  5V -> VCC
  GND -> GND

  The HX711 board can be powered from 2.7V to 5V so the Arduino 5V power should be fine.

*/


#include "HX711.h" //This library can be obtained here http://librarymanager/All#Avia_HX711

#define LOADCELL_DOUT_PIN 3  //define the Serial Data Output Pin
#define LOADCELL_SCK_PIN 2  // define the Power Down and Serial Clock Input Pin

HX711 scale;

int choose_branch; // initialize the variable "choose_branch"
float calibration_factor = -15400;
// initialized variable calibration_factor.
// Calibration factor found using a separate arduino sketch
// called SLoadCell_CalibrationFactorSketch"


void setup() {
  pinMode(4, OUTPUT);
  Serial.begin(115200);  // initialize arduino serial communication
  Serial.setTimeout(1);
  scale.begin(LOADCELL_DOUT_PIN, LOADCELL_SCK_PIN); // initialize load cell
  scale.set_scale(calibration_factor);
  //set the scale value;
  //this value is used to convert the raw
  //data to "human readable" data.  Output units are in lbs
  scale.set_offset(-17200);
  //This sets the offset value to a known zero.
  //There is no need for taring the scale once the zero point is known for a scale in a
  // set configuration.
  // The Zero Factor -17200 was found using a separate Arduino Sketch
  //called "SLoadCell_ZeroFactorSketch"
}

void loop() {
  char choose_branch = '0';
  int total = 0;
  if (Serial.available() > 0) {     // if information is sent over serial from matlab]
    choose_branch = Serial.read();  // read the serial data into variable "choose_branch"
    if (choose_branch == '2') {     // If choose_branch is equal to '2', iterate through the following for loop

      Serial.println("running"); // tell matlab that the arduino recieved the protocol id
      String reading = "2";
      // initiate variable to store the serial data when it is being read


      while (true) {
        reading = Serial.readString();                    //read from serial
        //if the reading is anything but 2 then break and start collecting data (we need to ignore 2 because matlab sends many instances
        if (reading.length() > 1 and reading != "2\n") {  //of the protocol id, so we want to make sure that it checks for the next unique value, which would be how many data points to collect)
          break;
        }
      }

      total = reading.toInt();     //convert the read string to an integer
      double start = millis();     //start timer
      for (int i = 0; i < total; i++) {
        digitalWrite(4,HIGH);
        Serial.println(scale.get_units(), 1); //scale.get_units() returns a float representing the force on the load cell
        Serial.println(analogRead(A0));       //reads raw pressure sensor data
        Serial.println(millis() - start);     //record time stamp of data collection
      }
      digitalWrite(4, LOW);
    
  } else {  //if there is no information to read over serial from matlab, wait...
  digitalWrite(4,LOW );
  }
}
}
