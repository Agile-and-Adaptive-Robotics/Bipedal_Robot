using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Web.Script.Serialization;
using SolidWorks.Interop.sldworks;
using SolidWorks.Interop.swconst;

class InspectAH001
{
    [STAThread]
    static void Main(string[] args)
    {
        AppDomain.CurrentDomain.AssemblyResolve += delegate(object sender, ResolveEventArgs a) {
            string name = new AssemblyName(a.Name).Name + ".dll";
            string root = @"D:\Program Files\SOLIDWORKS Corp\SOLIDWORKS";
            foreach (string dir in new string[] { Path.Combine(root,"api","redist"), Path.Combine(root,"URDFExporter") }) {
                string file=Path.Combine(dir,name);
                if(File.Exists(file)) return Assembly.LoadFrom(file);
            }
            return null;
        };
        Inspect(args[0]);
    }

    static void Inspect(string outDir)
    {
        SldWorks sw=(SldWorks)Marshal.GetActiveObject("SldWorks.Application");
        ModelDoc2 doc=(ModelDoc2)sw.ActiveDoc;
        if(doc==null || doc.GetTitle().IndexOf("10_AH_001",StringComparison.OrdinalIgnoreCase)<0)
            throw new Exception("The intended assembly is not active.");
        AssemblyDoc assembly=(AssemblyDoc)doc;
        var records=new List<object>();
        foreach(Component2 comp in (object[])assembly.GetComponents(false)) {
            ModelDoc2 part=comp.GetModelDoc2() as ModelDoc2;
            if(part==null || part.GetType()!=(int)swDocumentTypes_e.swDocPART) continue;
            object info;
            object[] bodies=comp.GetBodies3((int)swBodyType_e.swSolidBody,out info) as object[];
            var record=new Dictionary<string,object>();
            record["name"]=comp.Name2;
            record["path"]=comp.GetPathName();
            record["transform"]=comp.Transform2==null?null:comp.Transform2.ArrayData;
            record["component_box"]=comp.GetBox(false,false);
            record["body_count"]=bodies==null?0:bodies.Length;
            Console.WriteLine("Reading "+comp.Name2);
            var surfaces=new List<object>();
            if(bodies!=null && bodies.Length>0) {
                // Density-one body properties are purely geometric and avoid
                // assembly overrides. The model builder applies material densities.
                record["mass_properties_density_one"]=new List<object>();
                foreach(Body2 body in bodies) {
                    ((List<object>)record["mass_properties_density_one"]).Add(body.GetMassProperties(1.0));
                }
            }
            record["surfaces"]=surfaces;
            records.Add(record);
        }
        var output=new Dictionary<string,object>{{"source",doc.GetPathName()},{"components",records}};
        var serializer=new JavaScriptSerializer(); serializer.MaxJsonLength=Int32.MaxValue;
        File.WriteAllText(Path.Combine(outDir,"robot_mass_geometry.json"),serializer.Serialize(output));
        doc.ShowNamedView2("*Isometric",(int)swStandardViews_e.swIsometricView);
        doc.ViewZoomtofit2();
        int saveResult=doc.SaveAs3(Path.Combine(outDir,"robot_cad.png"),(int)swSaveAsVersion_e.swSaveAsCurrentVersion,
            (int)swSaveAsOptions_e.swSaveAsOptions_Silent | (int)swSaveAsOptions_e.swSaveAsOptions_Copy);
        Console.WriteLine("Inspected {0} component instances; image export status {1}",records.Count,saveResult);
    }
}
