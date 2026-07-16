#include "HX711.h"

HX711 scale;

const uint8_t DATA_PIN  = 3;
const uint8_t CLOCK_PIN = 2;
const int sensorPin = A0; // Analog Input from MPX5700GP Pin 
const int encoderPin = A1; // Analog Input from AS5600 Magnetic Encoder 

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
  float forceN = scale.get_units(1);
  float pressureRaw = analogRead(sensorPin);
  float pressureCalibrated = 0.7633 * pressureRaw - 13.744;
  float angleRaw = analogRead(encoderPin);
  float angleCalibrated = (angleRaw/1023)*360;
  Serial.print("Position = ");
  Serial.print(angleCalibrated, 3);
  Serial.print(" degrees ");
  Serial.print("Force = ");
  Serial.print(forceN, 3);
  Serial.print(" N ");
  Serial.print("Pressure = ");
  Serial.print(pressureCalibrated, 3);
  Serial.println(" kPa ");
  //delay(10);
}