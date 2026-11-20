#include "HX711.h"

// Valve Pin Assignments
const int fillPin = 8; // Input 1 of TBD62381APG (Slice 1, Solenoid 14)
const int exhaustPin = 9; // Input 2 of TBD62381APG (Slice 2, Solenoid 12)
const int sensorPin = A0; // Analog Input from MPX5700GP Pin 1

// HX711 Pin Assignments
const int hx711DataPin = 3; // DT Pin on HX711 Module
const int hx711ClockPin = 2; // SCK Pin on HX711 Module

HX711 scale;

// Target Window Constants (620 kPa +/- 10 kPa)
const float targetPressure = 620.0;
const float tolerance = 10.0;
const float lowerBound = targetPressure - tolerance; // 610.0 kPa
const float upperBound = targetPressure + tolerance; // 630.0 kPa

// MPX5700GP Calibration Parameters
const int adcOffset = 41;      
const int adcFullScale = 963;  

// HX711 Calibration Factor
// Change this value to match your load cell rating (grams, kg, or lbs calibration)
const float calibrationFactor = 420.0; 

void setup() {
  pinMode(fillPin, OUTPUT);
  pinMode(exhaustPin, OUTPUT);
  
  // Enforce safe default state on boot
  digitalWrite(fillPin, LOW);
  digitalWrite(exhaustPin, LOW);
  
  Serial.begin(9600);

  // Initialize the HX711 Scale
  scale.begin(hx711DataPin, hx711ClockPin);
  scale.set_scale(calibrationFactor); 
  scale.tare(); // Zero out the scale assuming it is empty on startup
}

void loop() {
  // 1. Read Pressure Sensor (MPX5700GP)
  long adcSum = 0;
  for (int i = 0; i < 20; i++) {
    adcSum += analogRead(sensorPin);
    delayMicroseconds(50); 
  }
  float averageAdc = (float)adcSum / 20.0;
  float currentPressure = (averageAdc - adcOffset) * 700.0 / (adcFullScale - adcOffset);
  if (currentPressure < 0.0) currentPressure = 0.0;

  // 2. Read Weight Scale (HX711)
  // use non-blocking raw check to see if HX711 data stream is ready
  float currentWeight = 0.0;
  if (scale.is_ready()) {
    currentWeight = scale.get_units(1); // Read a single data sample quickly
  }

  // 3. Bang-Bang Pressure Control Logic Loop
  if (currentPressure < lowerBound) {
    digitalWrite(fillPin, HIGH);    
    digitalWrite(exhaustPin, LOW);   
    Serial.print("[FILL] ");
  } 
  else if (currentPressure > upperBound) {
    digitalWrite(fillPin, LOW);     
    digitalWrite(exhaustPin, HIGH);  
    Serial.print("[EXHAUST] ");
  } 
  else {
    digitalWrite(fillPin, LOW);     
    digitalWrite(exhaustPin, LOW);   
    Serial.print("[HOLD] ");
  }

  // 4. Combined Metrics Output
  Serial.print("Pressure: ");
  Serial.print(currentPressure, 1);
  Serial.print(" kPa | ");
  Serial.print("Weight: ");
  Serial.println(currentWeight, 2);

  // Maintain sample frequency limit without choking the active scale data stream
  delay(100); 
}
