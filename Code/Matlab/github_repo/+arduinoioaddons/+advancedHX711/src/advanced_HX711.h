/**************************************************************************
 *   Copyright (C) 2019  Nicholas Giacoboni
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **************************************************************************
*/

#include "LibraryBase.h"

#define TEST 0x00
#define READ 0x01
#define POWER_DOWN 0x02
#define POWER_UP 0x03

class advanced_HX711 : public LibraryBase
{
    public:
        advanced_HX711(MWArduinoClass& a)
        {
            libName = "advancedHX711/advanced_HX711";
            a.registerLibrary(this);
        }
	public:
    void commandHandler(byte cmdID, byte* dataIn, unsigned int payloadSize)
    {
        switch (cmdID){
            case TEST:{
                /* Send back the "Hello World" string to check the 
                 *comunication, if everything works fine the string is
                 *shown in the command window;
                 */
                byte val [13] = "Hello World";  
                sendResponseMsg(cmdID, val, 13);
            break;
            } 
            case READ:{
                byte DOUT = dataIn[0];
                byte SCK = dataIn[1];
                byte Gain = dataIn[2];
                bool disable_interrupt = dataIn[3];
                byte data[3] = { 0 };

                digitalWrite(SCK, LOW);

                while (digitalRead(DOUT))
                {
                    /*Wait until the next output data is available
                     *according to the HX711 data sheet;
                     */
                }
                if(disable_interrupt) {
                    /*Serial interface between Arduino and the ADC 
                     *is time sensitive: diseable interrupt;
                     */
                    noInterrupts();
                }
                
                //Data retrieval phase
                data[2] = shiftIn(DOUT, SCK, MSBFIRST);
                data[1] = shiftIn(DOUT, SCK, MSBFIRST);
                data[0] = shiftIn(DOUT, SCK, MSBFIRST);
                
                for (int i = 0; i < Gain; ++i) {
                    digitalWrite(SCK, HIGH);
                    digitalWrite(SCK, LOW);
                }
                
                if(disable_interrupt) {
                    //Interrupt are now available;
                    interrupts();
                }

                data[2] ^= 0x80;

                sendResponseMsg(cmdID, data, 3);
            break;
            } 
            case POWER_DOWN:{
                byte DOUT = dataIn[0];
                byte SCK = dataIn[1];
                bool disable_interrupt = dataIn[2];
                
                if(disable_interrupt) {
                    noInterrupts();
                }
                
                digitalWrite(SCK, LOW);
                digitalWrite(SCK, HIGH);
                /*According to the datasheet: when SCK pin changes 
                 * from low to high and stays at high for longer 
                 * than 60µs, HX711 enters power down mode;
                 */
                
                if(disable_interrupt) {
                    //Interrupt are now available;
                    interrupts();
                }
                
                byte val [14] = "Powered down";  
                sendResponseMsg(cmdID, val, 14);
            break;
            }
            case POWER_UP:{
                byte DOUT = dataIn[0];
                byte SCK = dataIn[1];
                
                /*According to the datasheet: when SCK pin returns
                 * to low, chip will reset and enter normal 
                 * operation mode;
                 */ 
                digitalWrite(SCK, LOW);
                
                byte val [12] = "Powered up";  
                sendResponseMsg(cmdID, val, 12);
            break;
            }
            default:{
            }
        }
	}
};

