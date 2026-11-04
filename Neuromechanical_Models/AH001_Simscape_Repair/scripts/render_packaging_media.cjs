const {chromium}=require(process.env.AH001_NODE_MODULES+'/playwright');
const path=require('path'),fs=require('fs');
(async()=>{
 const output=path.resolve(__dirname,'../packaging/output'),frames=path.resolve(__dirname,'../packaging/build/frames');fs.mkdirSync(frames,{recursive:true});
 const browser=await chromium.launchPersistentContext(path.resolve(__dirname,'../packaging/build/media_profile'),{executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,viewport:{width:1280,height:960},args:['--no-sandbox','--allow-file-access-from-files','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 const page=await browser.newPage();await page.goto('file:///'+output.replaceAll('\\','/')+'/AH001_Packaging.html');await page.waitForFunction(()=>window.packagingReady,null,{timeout:60000});
 await page.evaluate(()=>{document.querySelectorAll('aside,#status,#help').forEach(x=>x.style.display='none');distance=2.85;target=[-.03,.65,0];yaw=.9;pitch=.12;render()});
 await page.screenshot({path:path.join(output,'Defense_Packaging.png')});
 await page.evaluate(()=>{document.getElementById('route').value='vas_int_r';document.getElementById('isolate').checked=true;document.getElementById('sync').checked=true;target=[-.06,.48,.02];distance=2.0;update()});
 for(let k=0;k<19;k++){
  await page.evaluate(k=>{document.getElementById('knee').value=k;update()},k);
  await page.screenshot({path:path.join(frames,`knee_${String(k).padStart(2,'0')}.png`)});
 }
 await page.evaluate(()=>{document.getElementById('knee').value=0;document.getElementById('route').value='ercspn_r';document.getElementById('isolate').checked=true;update();let r=selected(),P=routePoints(r);target=P.reduce((s,p)=>add(s,mul(p,1/P.length)),[0,0,0]);distance=.95;yaw=1.1;pitch=.10;render()});
 await page.screenshot({path:path.join(output,'Defense_Erector_Spinae.png')});
 await browser.close();console.log('Rendered clean stills and 19 checked knee frames.');
})().catch(e=>{console.error(e);process.exit(1)});
