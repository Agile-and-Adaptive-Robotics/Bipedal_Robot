/* Copyright 2018, Nicholas Giacoboni

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/

#include "LibraryBase.h"

class basic_HX711 : public LibraryBase
{
    public:
        basic_HX711(MWArduinoClass& a)
        {
            libName = "basicHX711/HX711";
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

