try { await import('./revise.mjs'); }
catch(e) {console.error(e.message);console.error((e.stack??'').split('\n').filter(x=>!x.includes('artifact_tool.mjs')).slice(0,10).join('\n'));process.exitCode=1;}
