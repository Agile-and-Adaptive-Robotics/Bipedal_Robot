import fs from 'node:fs/promises';
import path from 'node:path';

export async function apply(h){
  const {d,begin,newSlide,im,tx,table,replaceText,addAfter,a,old,GREEN,INK}=h;
  const f=name=>a(path.join('force',name));
  const cite='Modified from Bolen (2026)';
  const reference='Bolen, Elzein, Pang, and Hunt (2026), Actuators 15(5):230. Published source and submitted Chapter 3 are unchanged.';

  // Retain the original native quantitative chart on these two slides.
  replaceText(15,/112 mm initial length/,'⌀10 mm BPA   ·   l₀ = 112 mm   ·   P = 620 kPa');
  replaceText(25,/112 mm resting length/,'⌀10 mm BPA   ·   l₀ = 112 mm');

  begin(16,'Existing models mispredict force for ⌀10 mm and ⌀20 mm BPAs',
    'The measured force values remain on the horizontal axis. Only the proposed-model prediction points were removed. The plotted Sarosi and Martens–Boblan predictions are grouped as Others. This panel is the published 20 mm, 450 mm case; the 10 mm manufacturer example is on the previous slide, and the complete four-case numerical comparison follows later. The equality reference remains. '+reference,cite);
  im(16,old('slide17_comparison.png'),35,106,890,327);
  tx(16,'l₀ = 450 mm, ⌀20 mm BPA',255,441,450,28,23,false,INK,'center');

  for(const [n,diameter,title,b1,b2,b3] of [
    [19,10,'Shorter ⌀10 mm BPAs produce less maximum force','303.5 N','19.03 m⁻¹','0.0075 m'],
    [20,20,'The ⌀20 mm BPAs show the same loss of force at short initial lengths','922.4 N','15.37 m⁻¹','0.013 m']
  ]){
    begin(n,title,`Published maximum force at P = 620 kPa. The arctangent form is F620(l0) = b1 atan[b2(l0-b3)], with b1 = ${b1}, b2 = ${b2}, and b3 = ${b3}. The chart data, curves, and line colors are retained from the saved publication figure. The engineering diameter symbol replaces the old phi glyph. `+reference,cite);
    im(n,f(`maximum_force_${diameter}.png`),18,106,548,367);
    im(n,f('max_force_general_stacked.png'),570,151,364,116);
    tx(n,`b₁ = ${b1}\nb₂ = ${b2}\nb₃ = ${b3}`,603,294,326,123,26,false,INK);
  }

  begin(21,'An arbitrary pressure limit Pₘₐₓ scales the magnitude and shape of the maximum-force response',
    'This is the pressure-dependent maximum-force result for the 10 mm dataset. Pmax appears both in the prefactor and in the arctangent argument, so changing the supply limit changes magnitude and length-dependent shape. The plotted values and colors are unchanged. The later normalized force model retains F/F620, P/(620 kPa), and epsilon/epsilon620; do not normalize pressure a second time. The paper supports pressures below 620 kPa and describes higher-pressure use as unconfirmed extrapolation. The fixed 7.5 mm offset is verbatim from the published pressure equation. '+reference,cite);
  im(21,f('maximum_force_pressure_10.png'),18,111,542,360);
  im(21,f('max_force_pressure_stacked.png'),568,150,370,171);
  tx(21,'a₁ = 0.4895 N/kPa\na₂ = 0.03068 (kPa·m)⁻¹',583,354,350,78,24,false,INK);
  // The free-contraction measurement result on baseline 22 is retained.

  const original23=d.slides.getItem(22).speakerNotes.text;
  begin(23,'Normalization combines measurements across initial lengths',
    'The separate pressure key preserves the original blue, purple, pink, coral, and orange colors. Dashed lines are Festo predictions, solid lines are the new model, and dots are measured data. Definitions use the original 620 kPa force, pressure, and free-contraction references. Saved plot objects are re-exported without refitting or removing points. '+original23+' '+reference,cite);
  im(23,f('normalized_force_10_clean.png'),20,112,520,357);
  im(23,f('pressure_legend.png'),565,107,367,86);
  im(23,f('model_style_legend.png'),565,199,367,59);
  im(23,f('normalized_force_definition.png'),595,276,300,48);
  im(23,f('normalized_pressure_definition.png'),585,338,320,48);
  im(23,f('normalized_contraction_definition.png'),595,400,300,48);

  begin(24,'The normalized model captures force for both tested BPA diameters',
    'The generalized normalized equation is the published positive-force expression written as max(0, expression). The 10 mm coefficients are c0 = 0.5682, c1 = 4.254, and c2 = 0.5597; the 20 mm coefficients are c0 = 0.2579, c1 = 6.477, and c2 = 1.321. All source plot points, curves, and pressure colors are retained; repeated legends are removed as requested. Coefficients are checked against the final paper and the frozen submitted chapter. '+reference,cite);
  im(24,f('normalized_force.png'),37,109,886,43);
  im(24,f('normalized_force_10_clean.png'),25,164,429,270);
  im(24,f('normalized_force_20_clean.png'),500,164,429,270);
  tx(24,'⌀10 mm: c₀ = 0.5682, c₁ = 4.254, c₂ = 0.5597',28,443,437,29,19,false,INK);
  tx(24,'⌀20 mm: c₀ = 0.2579, c₁ = 6.477, c₂ = 1.321',499,443,437,29,19,false,INK);

  const comparison=JSON.parse(await fs.readFile(f('model_comparison_table.json'),'utf8'));
  const added=[];
  for(const diameter of [10,20]){
    const n=newSlide(`force_gof_${diameter}`,`The new model lowers prediction error for ⌀${diameter} mm BPAs`,
      'All values are transcribed from the published model-comparison table and independently matched to the frozen submitted Chapter 3. RMSE and maximum absolute error are in N; FVU is dimensionless. This is the reported comparison on these four actuator configurations, not a new optimization or a new held-out-actuator experiment. Smaller values indicate better agreement. '+reference,cite);
    const rows=comparison.long_rows.filter(row=>row.diameter_mm===diameter);
    const values=[['l₀ (mm)','Model','RMSE (N)','FVU','Max. error (N)'],...rows.map(row=>[String(row.l0_mm),row.model==='New model (Bolen)'?'New model':row.model,...row.values_display])];
    const t=table(n,values,35,128,890,289,[111,291,145,133,210],22);
    for(const r of [1,4])for(let c=0;c<5;c++)t.getCell(r,c).text.style={fontFamily:'Arial',fontSize:22,color:GREEN,alignment:c===0?'left':'center',bold:true};
    tx(n,'FVU: fraction of variance unexplained. Lower values indicate better agreement.',38,441,885,30,21,false,INK);
    added.push(n);
  }
  addAfter(27,...added);
  return {added};
}
