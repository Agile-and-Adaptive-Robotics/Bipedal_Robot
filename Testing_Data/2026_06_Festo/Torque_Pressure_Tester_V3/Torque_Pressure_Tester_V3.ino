#include "HX711.h"

HX711 scale;

// Valve Pin Assignments
const int fillPin = 11;       // Input 1 of TBD62381APG
const int exhaustPin = 6;     // Input 2 of TBD62381APG

// HX711 Pins
const uint8_t DATA_PIN  = 3;
const uint8_t CLOCK_PIN = 2;

// Analog inputs
const int sensorPin = A0;     // MPX5700GP pressure sensor
const int encoderPin = A2;    // AS5600 magnetic encoder

// Replace these with your printed calibration values
const int32_t OFFSET = -12228;
const float SCALE_FACTOR = -4489.861328;

// Maximum output rate.
// The actual rate may be lower if the HX711 is operating at 10 samples/s.
const unsigned long sampleIntervalMs = 50;

unsigned long currentTimer = 0;
unsigned long previousTimer = 0;


void sendValveState()
{
  Serial.print("STATE,");
  Serial.print(digitalRead(fillPin));
  Serial.print(",");
  Serial.println(digitalRead(exhaustPin));
}


void setValves(bool valveState)
{
  digitalWrite(fillPin, valveState ? HIGH : LOW);
  digitalWrite(exhaustPin, valveState ? HIGH : LOW);

  // Immediately report the new output states
  sendValveState();
}


void handleSerialCommands()
{
  while (Serial.available() > 0)
  {
    char command = Serial.read();

    // Accept either upper- or lower-case commands
    if (command >= 'a' && command <= 'z')
    {
      command = command - ('a' - 'A');
    }

    switch (command)
    {
      case 'V':
        // Both valve-control outputs HIGH
        setValves(true);
        break;

      case 'O':
        // Both valve-control outputs LOW
        setValves(false);
        break;

      case '?':
        // MATLAB can request the present valve state
        sendValveState();
        break;

      default:
        // Ignore line endings and unrecognized characters
        break;
    }
  }
}


void setup()
{
  pinMode(fillPin, OUTPUT);
  pinMode(exhaustPin, OUTPUT);

  // Enforce safe default state on boot
  digitalWrite(fillPin, LOW);
  digitalWrite(exhaustPin, LOW);

  Serial.begin(115200);

  scale.begin(DATA_PIN, CLOCK_PIN);
  scale.set_offset(OFFSET);
  scale.set_scale(SCALE_FACTOR);

  sendValveState();
}


void loop()
{
  // Check commands on every pass through loop()
  handleSerialCommands();

  currentTimer = millis();

  // is_ready() prevents get_units() from waiting for a new HX711 sample
  if ((unsigned long)(currentTimer - previousTimer) >= sampleIntervalMs &&
      scale.is_ready())
  {
    previousTimer = currentTimer;

    float forceN = scale.get_units(1);

    float pressureRaw = analogRead(sensorPin);
    float pressureCalibrated = 0.7633 * pressureRaw - 13.744;

    float angleRaw = analogRead(encoderPin);
    float angleCalibrated = -(angleRaw * 0.3629) + 216.74;

    int fillState = digitalRead(fillPin);
    int exhaustState = digitalRead(exhaustPin);
    // Time
    Serial.print("Time = ");
    Serial.print(currentTimer);
    Serial.print(",");
    // Angle raw
    Serial.print("Angle (raw) = ");
    Serial.print(angleRaw, 3);
    Serial.print(",");
    // Angle degree
    Serial.print("Angle (deg) = ");
    Serial.print(angleCalibrated, 3);
    Serial.print(",");
    // Force (N)
    Serial.print("Force (N) = ");
    Serial.print(forceN, 3);
    Serial.print(",");
    // Pressure (kPa)
    Serial.print("Pressure (kPa) = ");
    Serial.print(pressureCalibrated, 3);
    Serial.print(",");
    // Fill Valve Pin State
    Serial.print("Fill Valve State = ");
    Serial.print(fillState);
    Serial.print(",");
    // Hold Valve Pin State
    Serial.print("Hold Valve State = ");
    Serial.println(exhaustState);
  }
}
