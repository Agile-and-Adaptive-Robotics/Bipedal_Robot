[![View Advanced Custom Arduino Library for HX711 on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/73692-advanced-custom-arduino-library-for-hx711)
# Adv_UNO_HX711_Matlab
Get data from the HX711 load cell amplifier and import into MATLAB workspace 
with Arduino Uno and Mega2560.

This manual explains how to register a custom library and get started with
the HX711 load cell ampliﬁer in MATLAB environment. You need to have
already installed the MATLAB Support Package for ARDUINO Hardware.
It does require basic knowledge about MATLAB and his functionality, furthermore 
it assumes your load cell works properly and the connection is
correct. This particular library has been tested with Arduino UNO and
MEGA2560, if you are using Arduino DUE you can use this one:
https://github.com/GiacoboniNicholas/Adv_DUE_HX711_Matlab.

***IMPORTANT: if you are not familiar with custom Arduino Addons I suggest to download 
my basic Addons for HX711 here:
https://it.mathworks.com/matlabcentral/fileexchange/66641-basic-custom-arduino-library-for-hx711
which contains a more complete User's Manual. Once you get confident with the basic one 
you can use this library which is more complete and offers several functionality such as 
PowerUp and PowerDown function, a more stable serial comunication, an improved calibration
class and so on...***

## Overview
An add-on library is a collection of MATLAB and C++ code that provides a user
easy access to features on the Arduino hardware or attached shields. 
The HX711 add-on library just develops the two wire comunication protocol
between Matlab workspace and HX711 itself through ARDUINO, thus you can built
your own commands to calibrate your load cell and acquire weight. 

The HX711 is a precision 24-bit analog to digital converter (ADC) designed for 
weigh scales and industrial control applications to interface directly with a 
bridge sensor. All you have to do is to specify the name of the digital pins in
according to phisical connection. 

## Limitations
• This library has been tested with ARDUINO UNO and MEGA2560 under windows environment. 
  It does not work with ARDUINO DUE and it might not work with other board. If you are
  using Arduino DUE use this library https://github.com/GiacoboniNicholas/Adv_DUE_HX711_Matlab
  
• You can make static or dynamic measurements with high precision even with low cost 
  components, but in the second case do not exceed the sampling rate of 10 [Hz] because 
  you may have some comunication error (0 or NaN).

## Caution
Here I sum up some tips you should follow in order to avoid wasting time and solve every
problem quickly. Please don’t ask for help in comment section if you are not going to 
check the following recommendations.

• I suggest to test your load cell with an usual library in Arduino IDE environment. 
  You can ﬁnd a very good one at https://github.com/bogde/HX711. 
  After the calibration phase take note of the scale factor so you can compare it 
  with the one you get from MATLAB.
  
• Make sure that the connection between ARDUINO and HX711 is correct and particularly 
  the connection between the load cell and HX711 must be stable and strong 
  (I personally recommend a soldered connection, it also reduces noise eﬀect if it’s 
  well executed, instead of other temporary connection).
  
• The correct functioning of the load cell depends on the installation of the strain 
  gauges, I warn you that on the market there are load cell with a low quality installation 
  technique and even a wrong bridge connection between strain gauges. If the calibration 
  phase is done properly you should be able to make measurements with at least 0.001 [kg] of
  precision with a 0÷10 [kg] load cell. In the best case the error is below 0.0001 [kg]. 

## How to Start
In this part you will learn how to use a Custom Arduino Libraries and then how to get data 
from HX711 with just one command. 

### Register Custom Library
The current folder contains two ﬁles named HX711.m and HX711.h. The ﬁrst one is the MATLAB 
Add-On Class that inherits from arduino LibraryBase class a variety of properties and methods. 
The C++ Header File is a class that includes librarybase.h and the code segments executed on the
Arduino device.

First you should extract content from .zip ﬁle and paste it in your own personal folder (C:\work). 
In order to register the HX711 custom library you have to add the path of the folder to Matlab search path. 
You can do that with the following command:
```c++
1 addpath ( ’C : \ work ’ );
```
The structure in the work folder is the following: +arduinoioaddons which contains several subfolders, 
one of these is +advancedHX711 folder. You can have other +NameLibrary subfolder in +arduinoioaddons folder. 
The folder structure is really important so if you have some questions look here:
https://it.mathworks.com/help/supportpkg/arduinoio/ug/create-custom-folder-structure.html. 
Make sure the basicHX711/basic HX711 library is available with listArduinoLibraries command.
```c++  
2 listArduinoLibraries  
```  

### Upload Arduino Server
Now you have to initialize the connection between Matlab and Arduino. You can do that in two ways. 
The ﬁrst one consists of using the following command, which requires as argument the Arduino Board 
you’re currently using and the name of serial port (UNO and COM6 in the example).
```c++  
3 a= arduino (’com6’,’Uno’,’libraries’,’advancedHX711/advanced_HX711’);
```  
The second way consists of using the hardware setup procedure which is initialized by the following command:
```c++  
3 arduinosetup
```
Now you can start creating the objects of the classes: 
```c++  
4 a= arduino (’com6’,’Uno’) // If you ’ve used Hardware Setup ;
5 LoadCell = addon (a,’advancedHX711/advanced_HX711’,’Pins’,{’D2’,’D3’})
``` 
The default gain is 128, if you wish you can use 64 for channel A and 32 for channel B with an additional parameter:
```c++  
5 LoadCell = addon (a,’advancedHX711/advanced_HX711’,’Pins’,{’D2’,’D3’},’Gain’,64)
``` 
The serial comunication protocol is time sensitive so if you need to use Interrupt set "true" 
with an additional parameter: 
```c++  
5 LoadCell = addon (a,’advancedHX711/advanced_HX711’,’Pins’,{’D2’,’D3’},’Interrupt’,true)
``` 
Finally you can get data from HX711 with the following command:
```c++  
6 read_HX711(LoadCell)
``` 
Now you can write your own function to calibrate the loadcell or use the calibration class.

### Troubleshooting
If Custom Arduino Library class is not detected visit Custom Arduino Library Issues page: https://it.mathworks.com/help/supportpkg/arduinoio/ug/custom-arduino-library-issues.html.

## Calibration Class
The calibration class allows you to calibrate your loadcell with a few simple builtin command.

•***1.*** Create the object of the calibration class:
  ```c++  
    cal = calibration(number_of_readings,known_weight); 
  ```
  The first parameter is the number of readings during calibration phase, in order to have a more
  precise calibration every partial results is alwais an averege of "number_of_readings" values.
  ***I suggest at least 100 readings.***
  The second parameter is the weight of your sample which is necessary to determine the scale factor.
  ***The unit of measure is nested in the scale factor, so if you set the known weight in [kg] then 
  the result you'll get after calibration phase will also be in [kg]***.
  
  ***Example (100 readings, weight 250 grams):***
  ```c++
    1 cal = calibration(100,250)
  ``` 
  
•***2.*** Call the tare function:
  ```c++  
    tare(name_calibration_object,name_HX711_object); 
  ```
  the first parameter is the name of calibration's object (in the point 1 I've chosen "cal"), the 
  second parameter is the name of HX711's object (in the previous paragraph I've chosen "LoadCell").
  
  ***Example:***
  ```c++  
    2 tare(cal,LoadCell); 
  ```
  
•***3.*** Call the scale function:
 ```c++  
    scale(name_calibration_object,name_HX711_object); 
  ```
  ***Example:***
  ```c++  
    3 scale(cal,LoadCell); 
  ```
  
•***4.*** The calibration phase is done. Now you can call the get_weight function which give you a single raw reading:
  ```c++  
     get_weight(name_calibration_object,name_HX711_object); 
  ```
   ***Example:***
  ```c++  
    4 get_weight(cal,LoadCell);
  ```  
  If you want an average value of multiple readings add an additional parameter: 
  ```c++  
     get_weight(name_calibration_object,name_HX711_object,number_of_readings); 
  ```
   ***Example:***
  ```c++  
    4 get_weight(cal,LoadCell,10);
  ```  
  
•***OPT***. The number of readings and the known weight are public properties, so if you want to change them 
  after the object is created you can do it in the following way:
  ```c++  
  5 cal.n = new_number_of_readings; 
  6 cal.known_weight = new_known_weight;
  ```
  now you have to repeat steps 2 and 3.
  
 •***OPT.*** If you are interested in how your loadcell behaves with a speciﬁc static load in terms of precision 
   and accuracy you can use the stat function in order to get the average and the standard deviation of 
   multiple readings: 
   ```c++  
      stat(name_calibration_object,name_HX711_object,number_of_readings);
   ```    
   ***Example:***    
   ```c++  
     7 stat(cal,LoadCell,100);
   ```   
    
•***OPT.*** Finally if you are interested in a visual representation of the stat function you can use the 
  plot data function:
  ```c++  
     plot_data(name_calibration_object,name_HX711_object,number_of_readings,known_weight);
  ```
  ***Example:***
  ```c++  
      8 plot_data(cal,LoadCell,100,250);
  ```   
***Important***: according to the Avia semiconductor's HX711 datasheet the output 24 bits of data is in 
2’s complement format. When input differential signal goes out of the 24 bit range, the output data will 
be saturated at 800000h (MIN) or 7FFFFFh (MAX), until the input signal comes back to the input range. 
In order to garantee a more accurate calibration of your loadcell when the Calibration class starts 
the serial comunication with HX711 it automatically checks the output, if the data is 800000h or 7FFFFFh the 
Calibration class assigns NaN. During calibration phase NaN values are automatically discarded. If you are
acquiring multiple readings with get_weight function it will be easier to filter the data. 








  
  
