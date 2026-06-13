#include "LibraryBase.h"

class basic_HX711 : public LibraryBase
{
    public:
        basic_HX711(MWArduinoClass& a)
        {
            // Must match the MATLAB LibraryName exactly.
            libName = "basicHX711/basic_HX711";
            a.registerLibrary(this);
        }
    public:
    void commandHandler(byte cmdID, byte* dataIn, unsigned int payloadSize)
    {
        switch (cmdID){
            case 0x01:{
                byte DOUT = dataIn[0];
                byte SCK = dataIn[1];
                byte data[3] = { 0 };

                digitalWrite(SCK, LOW);

                while (digitalRead(DOUT))
                {
                    // Wait
                }

                data[2] = shiftIn(DOUT, SCK, MSBFIRST);
                data[1] = shiftIn(DOUT, SCK, MSBFIRST);
                data[0] = shiftIn(DOUT, SCK, MSBFIRST);

                digitalWrite(SCK, HIGH);
                digitalWrite(SCK, LOW);

                data[2] ^= 0x80;

                sendResponseMsg(cmdID, data, 3);
            break;
            }
            default:{
            }
        }
    }
};
