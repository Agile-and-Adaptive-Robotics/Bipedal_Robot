#include "HX711.h"

HX711 scale;

// Valve Pin Assignments
const int fillPin = 11; // Input 1 of TBD62381APG (Slice 1, Solenoid 14)
const int exhaustPin = 6; // Input 2 of TBD62381APG (Slice 2, Solenoid 12)
//HX711 Pins
const uint8_t DATA_PIN  = 3;
const uint8_t CLOCK_PIN = 2;
//Analog reading
const int sensorPin = A0; // Analog Input from MPX5700GP Pin 
const int encoderPin = A2; // Analog Input from AS5600 Magnetic Encoder 
 unsigned long  currentTimer = 0;
  float previousTimer = 0;
// Replace these with your printed calibration values
const int32_t OFFSET = -12228;
const float SCALE_FACTOR = -4489.861328;

void setup(){
  pinMode(fillPin, OUTPUT);
  pinMode(exhaustPin, OUTPUT);
  
  // Enforce safe default state on boot
  digitalWrite(fillPin, LOW);
  digitalWrite(exhaustPin, LOW);
  
  Serial.begin(115200);

  scale.begin(DATA_PIN, CLOCK_PIN);

  scale.set_offset(OFFSET);
  scale.set_scale(SCALE_FACTOR);
}

void loop()
{
  currentTimer = millis();
  //if (currentTimer - previousTimer >=1){
  float forceN = scale.get_units(1);
  float pressureRaw = analogRead(sensorPin);
  float pressureCalibrated = 0.7572 * pressureRaw - 13.955;
  float angleRaw = analogRead(encoderPin);
  float angleCalibrated = -(angleRaw*0.359)+215.25;
  // Serial.print("Time = ");
  Serial.print(currentTimer);
  // Serial.print(" ms ");
  Serial.print(",");
  // Serial.print("Position = ");
  Serial.print(angleCalibrated, 3);
  // Serial.print(" degrees ");
  Serial.print(",");
  // Serial.print("Force = ");
  Serial.print(forceN, 3);
  // Serial.print(" N ");
  Serial.print(",");
  // Serial.print("Pressure = ");
  Serial.println(pressureCalibrated, 3);
  // Serial.println(" kPa ");
  //delay(10);
  //previousTimer = currentTimer; 

}
