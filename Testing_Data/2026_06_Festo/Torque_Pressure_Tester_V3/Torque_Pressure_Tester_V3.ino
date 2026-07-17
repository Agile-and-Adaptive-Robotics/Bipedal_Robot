#include "HX711.h"

HX711 scale;

const uint8_t DATA_PIN  = 3;
const uint8_t CLOCK_PIN = 2;
const int sensorPin = A0; // Analog Input from MPX5700GP Pin 
const int encoderPin = A2; // Analog Input from AS5600 Magnetic Encoder 
 unsigned long  currentTimer = 0;
  float previousTimer = 0;
// Replace these with your printed calibration values
const int32_t OFFSET = -12228;
const float SCALE_FACTOR = -4489.861328;

void setup()
{
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
  float pressureCalibrated = 0.7633 * pressureRaw - 13.744;
  float angleRaw = analogRead(encoderPin);
  float angleCalibrated = -(angleRaw*0.3629)+216.74;
  //Serial.print("Time = ");
  //Serial.print(currentTimer);
  //Serial.print(" ms ");
  Serial.print("Position (raw) = ");
  Serial.print(angleRaw, 3);
  Serial.print("Position (deg) = ");
  Serial.print(angleCalibrated, 3);
  //Serial.print(" degrees ");
  Serial.print("Force (N) = ");
  Serial.print(forceN, 3);
  //Serial.print(" N ");
  Serial.print("Pressure (kPa) = ");
  Serial.println(pressureCalibrated, 3);
  //Serial.println(" kPa ");
  //delay(10);
  //previousTimer = currentTimer; 

}
