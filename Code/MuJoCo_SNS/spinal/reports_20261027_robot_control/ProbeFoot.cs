using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Web.Script.Serialization;
using SolidWorks.Interop.sldworks;
using SolidWorks.Interop.swconst;
class ProbeFoot {
 [STAThread] static void Main(string[] args) {
  AppDomain.CurrentDomain.AssemblyResolve += delegate(object sender,ResolveEventArgs a) {
   string p=Path.Combine(@"D:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\api\redist",new AssemblyName(a.Name).Name+".dll");
   return File.Exists(p)?Assembly.LoadFrom(p):null;
  }; Inspect(args);
 }
 static void Inspect(string[] args) {
  SldWorks sw=(SldWorks)Marshal.GetActiveObject("SldWorks.Application");
  int errors=0,warnings=0;
  ModelDoc2 doc=(ModelDoc2)sw.OpenDoc6(args[0],(int)swDocumentTypes_e.swDocASSEMBLY,
   (int)swOpenDocOptions_e.swOpenDocOptions_Silent|(int)swOpenDocOptions_e.swOpenDocOptions_ReadOnly,"",ref errors,ref warnings);
  if(doc==null) throw new Exception("Foot open failed: "+errors+" warnings "+warnings);
  var records=new List<object>();
  foreach(Component2 comp in (object[])((AssemblyDoc)doc).GetComponents(false)) {
   var record=new Dictionary<string,object>{{"name",comp.Name2},{"path",comp.GetPathName()},{"suppressed",comp.IsSuppressed()},{"transform",comp.Transform2==null?null:comp.Transform2.ArrayData}};
   object info; object[] bodies=comp.GetBodies3((int)swBodyType_e.swSolidBody,out info) as object[];
   var masses=new List<object>();
   if(bodies!=null) foreach(Body2 body in bodies) masses.Add(body.GetMassProperties(1.0));
   record["mass_properties_density_one"]=masses;
   records.Add(record);
  }
  var output=new Dictionary<string,object>{{"source",doc.GetPathName()},{"open_errors",errors},{"open_warnings",warnings},{"box",((AssemblyDoc)doc).GetBox((int)swBoundingBoxOptions_e.swBoundingBoxIncludeRefPlanes)},{"components",records}};
  var serializer=new JavaScriptSerializer(); serializer.MaxJsonLength=Int32.MaxValue;
  File.WriteAllText(Path.Combine(args[1],"foot_inventory.json"),serializer.Serialize(output));
  doc.ShowNamedView2("*Isometric",(int)swStandardViews_e.swIsometricView); doc.ViewZoomtofit2();
  doc.SaveAs3(Path.Combine(args[1],"foot_cad.png"),0,(int)swSaveAsOptions_e.swSaveAsOptions_Silent|(int)swSaveAsOptions_e.swSaveAsOptions_Copy);
  bool binary=sw.GetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLBinaryFormat);
  bool translate=sw.GetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLDontTranslateToPositive);
  bool merged=sw.GetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLComponentsIntoOneFile);
  int units=sw.GetUserPreferenceIntegerValue((int)swUserPreferenceIntegerValue_e.swExportStlUnits);
  try {
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLBinaryFormat,true);
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLDontTranslateToPositive,true);
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLComponentsIntoOneFile,true);
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLShowInfoOnSave,false);
   sw.SetUserPreferenceIntegerValue((int)swUserPreferenceIntegerValue_e.swExportStlUnits,2);
   int status=doc.SaveAs3(Path.Combine(args[1],"foot_cad.stl"),0,(int)swSaveAsOptions_e.swSaveAsOptions_Silent|(int)swSaveAsOptions_e.swSaveAsOptions_Copy);
   Console.WriteLine("Foot components {0}; open errors {1}; warnings {2}; STL status {3}",records.Count,errors,warnings,status);
  } finally {
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLBinaryFormat,binary);
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLDontTranslateToPositive,translate);
   sw.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swSTLComponentsIntoOneFile,merged);
   sw.SetUserPreferenceIntegerValue((int)swUserPreferenceIntegerValue_e.swExportStlUnits,units);
  }
 }
}
