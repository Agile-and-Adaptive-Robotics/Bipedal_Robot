#include "LibraryBase.h"

const char MSG_EXAMPLE_DEBUG[]  PROGMEM = "Example debug message: cmdID %d\n";

class HelloWorld : public LibraryBase 
{  
    public:
        HelloWorld(MWArduinoClass& a)
        {            
            libName = "ExampleAddon/HelloWorld";            
            a.registerLibrary(this);
        }
    
    public:
   void commandHandler(byte cmdID, byte* dataIn, unsigned int payloadSize)
   {            
     
      debugPrint(MSG_EXAMPLE_DEBUG, cmdID);
      switch (cmdID){
               
         case 0x01:{  
            byte val [13] = "Hello World!";                  
            sendResponseMsg(cmdID, val, 13);
            break;
         }

	  // Other cases with appropriate cmdIDs  

         default:{
            // Do nothing
         }
      }
   }


};