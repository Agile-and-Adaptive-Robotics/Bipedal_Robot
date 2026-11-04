const {chromium}=require(process.env.AH001_NODE_MODULES+'/playwright');
const path=require('path');
const fs=require('fs');
(async()=>{
 const output=path.resolve(__dirname,'../packaging/output');
 const browser=await chromium.launchPersistentContext(path.resolve(__dirname,'../packaging/build/browser_profile'),{executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,viewport:{width:1600,height:1100},deviceScaleFactor:1,args:['--no-sandbox','--disable-crash-reporter','--allow-file-access-from-files','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 const page=await browser.newPage();
 const errors=[];page.on('pageerror',e=>errors.push(String(e)));
 await page.goto('file:///'+output.replaceAll('\\','/')+'/AH001_Packaging.html');
 try {await page.waitForFunction(()=>window.packagingReady,null,{timeout:60000});}
 catch(e){console.log('PAGE',await page.locator('#loading').textContent(),errors);await page.screenshot({path:path.join(output,'failed_load.png')});await browser.close();throw e}
 await page.screenshot({path:path.join(output,'preview.png')});
 const neutral=await page.evaluate(()=>window.packagingState);
 await page.selectOption('#route','ercspn_r');await page.check('#isolate');await page.click('#focus');
 await page.screenshot({path:path.join(output,'erector_spinae.png')});
 await page.evaluate(()=>{document.getElementById('i').value=0;document.getElementById('j').value=0;update()});
 const A=await page.evaluate(()=>window.packagingState);
 await page.evaluate(()=>{document.getElementById('i').value=4;document.getElementById('j').value=4;update()});
 const B=await page.evaluate(()=>window.packagingState);
 await page.selectOption('#route','vas_int_r');await page.evaluate(()=>{document.getElementById('knee').value=18;update()});
 await page.click('#reset');await page.screenshot({path:path.join(output,'knee_flexion.png')});
 const knee=await page.evaluate(()=>window.packagingState);
 if(errors.length)throw Error(errors.join('\n'));
 if(A.diameter_mm===B.diameter_mm)throw Error('Saved state did not change diameter');
 if(knee.knee!==18)throw Error('Knee slider failed');
 const result={neutral,erector_states:[A,B],knee,errors};
 fs.writeFileSync(path.join(output,'browser_verification.json'),JSON.stringify(result,null,2));
 console.log(JSON.stringify(result,null,2));await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
