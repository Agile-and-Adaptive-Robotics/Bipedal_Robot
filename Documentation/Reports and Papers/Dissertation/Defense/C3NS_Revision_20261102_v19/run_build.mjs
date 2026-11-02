try { await import('./rebuild.mjs'); }
catch(e) { console.error(e.message); console.error((e.stack||'').split('\n').filter(x=>!x.includes('artifact_tool.mjs')).slice(0,8).join('\n')); process.exitCode=1; }
