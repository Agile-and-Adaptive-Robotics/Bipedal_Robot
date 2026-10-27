using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Web.Script.Serialization;
using SolidWorks.Interop.sldworks;
using SW2URDF.URDF;
using SW2URDF.URDFExport;

class ExportAH001
{
    class PackageNotice : SW2URDF.UI.IMessageBox {
        public System.Windows.MessageBoxResult Show(string message) {
            if(!message.StartsWith("Creating URDF Package")) throw new Exception(message);
            Console.WriteLine(message); return System.Windows.MessageBoxResult.OK;
        }
        public System.Windows.MessageBoxResult Show(string message,string caption,System.Windows.MessageBoxButton buttons) {
            throw new Exception("Unexpected exporter question: "+message);
        }
    }
    public class Group { public string name; public string[] components; }
    [STAThread]
    static void Main(string[] args)
    {
        AppDomain.CurrentDomain.AssemblyResolve += delegate(object sender, ResolveEventArgs a) {
            string name=new AssemblyName(a.Name).Name+".dll";
            string root=@"D:\Program Files\SOLIDWORKS Corp\SOLIDWORKS";
            foreach(string dir in new string[]{Path.Combine(root,"api","redist"),Path.Combine(root,"URDFExporter")}) {
                string file=Path.Combine(dir,name);
                if(File.Exists(file)) return Assembly.LoadFrom(file);
            }
            return null;
        };
        Export(args[0],args[1]);
    }
    static void Export(string configuration,string destination)
    {
        SldWorks sw=(SldWorks)Marshal.GetActiveObject("SldWorks.Application");
        ModelDoc2 doc=(ModelDoc2)sw.ActiveDoc;
        if(doc==null || doc.GetTitle().IndexOf("10_AH_001",StringComparison.OrdinalIgnoreCase)<0)
            throw new Exception("The intended source assembly must be active.");
        AssemblyDoc assembly=(AssemblyDoc)doc;
        var components=new Dictionary<string,Component2>();
        foreach(Component2 component in (object[])assembly.GetComponents(false)) components[component.Name2]=component;
        var groups=new JavaScriptSerializer().Deserialize<Group[]>(File.ReadAllText(configuration));
        Robot robot=new Robot();
        Link root=null;
        foreach(Group group in groups) {
            Link link=new Link(); link.Name=group.name; link.STLQualityFine=true;
            link.Joint.CoordinateSystemName=""; // Original assembly frame for each mesh.
            // Installed exporter embeds interop types in its generic list;
            // IList avoids the C# compiler's cross-assembly No-PIA restriction.
            System.Collections.IList assigned=(System.Collections.IList)typeof(Link).GetField("SWComponents").GetValue(link);
            foreach(string name in group.components) assigned.Add(components[name]);
            link.SWMainComponent=components[group.components[0]];
            // This intermediate package only groups meshes. The separate builder
            // assigns audited manual axes, inertia, and four-bar loop closures.
            if(root==null) { root=link; robot.SetBaseLink(root); }
            else {
                link.Parent=root; root.Children.Add(link);
                link.Joint.Name="mesh_group_"+group.name;
                link.Joint.Type="fixed";
                link.Joint.Parent.Name=root.Name; link.Joint.Child.Name=link.Name;
                link.Joint.Origin.SetXYZ(new double[]{0,0,0});
                link.Joint.Origin.SetRPY(new double[]{0,0,0});
            }
        }
        ExportHelper helper=new ExportHelper(sw);
        helper.SetComputeInertial(false); helper.SetComputeJointKinematics(false);
        helper.SetComputeJointLimits(false); helper.SetComputeVisualCollision(false);
        helper.URDFRobot=robot; helper.PackageName="ah001_cad_meshes"; helper.SavePath=destination;
        URDFPackage.MessageBox=new PackageNotice();
        helper.ExportRobot(true);
        Console.WriteLine("Exported {0} manually assigned rigid mesh groups without saving source CAD.",groups.Length);
    }
}
