// DengExport.cs - headless AnimatLab standalone-sim exporter for the Deng 2019
// original project. Must live in the AnimatLab bin dir (the app resolves its
// modules against the entry-exe directory). Loads the project via the app's
// own -PROJECT command-line handling inside StartApplication, then calls
// ExportStandAloneSim (inline if StartApplication returns, else via the UI
// thread). Auto-dismisses Win32 "Error" dialogs. Never saves the project.
using System;
using System.IO;
using System.Reflection;
using System.Threading;
using System.Runtime.InteropServices;
using System.Text;
using System.Collections.Generic;

class DengExport {
    [DllImport("user32.dll")] static extern bool EnumWindows(EnumWindowsProc cb, IntPtr lp);
    [DllImport("user32.dll")] static extern bool EnumChildWindows(IntPtr h, EnumChildProc cb, IntPtr lp);
    [DllImport("user32.dll")] static extern int GetClassName(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll")] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
    [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] static extern IntPtr SendMessage(IntPtr h, uint msg, IntPtr wp, IntPtr lp);
    delegate bool EnumWindowsProc(IntPtr h, IntPtr lp);
    delegate bool EnumChildProc(IntPtr h, IntPtr lp);

    static uint _pid;
    static volatile bool _stop;
    static object _app;
    static string _aproj;
    static string _bin;
    [ThreadStatic] static bool _inResolve;

    const BindingFlags FPUB = BindingFlags.InvokeMethod | BindingFlags.GetProperty | BindingFlags.Public | BindingFlags.Instance;
    const BindingFlags FPRIV = BindingFlags.InvokeMethod | BindingFlags.GetProperty | BindingFlags.NonPublic | BindingFlags.Instance;

    static object CallAny(object target, string member, object[] par) {
        var t = target.GetType();
        const BindingFlags ALL = BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance;
        var mi = t.GetMethod(member, ALL);
        if (mi == null) {
            foreach (var i in t.GetInterfaces()) {
                var map = t.GetInterfaceMap(i);
                for (int k = 0; k < map.InterfaceMethods.Length; k++) {
                    if (map.InterfaceMethods[k].Name == member) { mi = map.TargetMethods[k]; break; }
                }
                if (mi != null) break;
            }
        }
        if (mi == null) throw new MissingMethodException(t.FullName, member);
        return mi.Invoke(target, par ?? new object[0]);
    }

    static void Log(string s) { Console.Out.WriteLine(s); Console.Out.Flush(); }

    [STAThread]
    static void Main(string[] args) {
        // accept: DengExport <path>  OR  DengExport -PROJECT <path> (pass-through form)
        string path = null;
        for (int i = 0; i < args.Length; i++) {
            if (args[i].ToUpper() == "-PROJECT" && i + 1 < args.Length) path = args[++i];
            else if (path == null) path = args[i];
        }
        if (path == null || !File.Exists(path)) { Log("usage: DengExport [-PROJECT] <project.aproj>"); Environment.Exit(1); }
        _aproj = Path.GetFullPath(path);
        _bin = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
        Directory.SetCurrentDirectory(_bin);
        Log("bin=" + _bin);
        Log("aproj=" + _aproj);

        AppDomain.CurrentDomain.AssemblyResolve += (s, e) => {
            if (_inResolve) return null;
            _inResolve = true;
            try {
                string p = Path.Combine(_bin, e.Name.Split(',')[0] + ".dll");
                return File.Exists(p) ? Assembly.LoadFrom(p) : null;
            } finally { _inResolve = false; }
        };

        _pid = (uint)System.Diagnostics.Process.GetCurrentProcess().Id;
        var killer = new Thread(KillErrorDialogs); killer.IsBackground = true; killer.Start();

        // all app work happens on a dedicated STA UI thread that pumps messages
        var ui = new Thread(UiMain);
        ui.SetApartmentState(System.Threading.ApartmentState.STA);
        ui.Start();
        // wait for the app form's handle (created inside StartApplication's modal loop)
        var sw = System.Diagnostics.Stopwatch.StartNew();
        System.Windows.Forms.Form form = null;
        while (sw.Elapsed.TotalSeconds < 120) {
            form = _app as System.Windows.Forms.Form;
            if (form != null && form.IsHandleCreated) break;
            Thread.Sleep(500);
        }
        if (form == null || !form.IsHandleCreated) { Log("TIMEOUT: no form handle"); _stop = true; Thread.Sleep(500); Environment.Exit(6); }
        Log("form handle up after " + sw.Elapsed.TotalSeconds.ToString("F1") + "s; settling 20s...");
        Thread.Sleep(20000);
        Log("queuing work");
        form.BeginInvoke((Action)delegate {
            Log("work: LoadProject...");
            try {
                CallAny(_app, "LoadProject", new object[]{ _aproj });
                Thread.Sleep(3000);
                Log("ProjectIsOpen=" + SafeProp("ProjectIsOpen"));
                Log("ProjectName=" + SafeProp("ProjectName"));
                Log("ProjectPath=" + SafeProp("ProjectPath"));
                Log("work: ExportStandAloneSim regardless...");
                ExportNow();
            } catch (Exception ex) { Log("work error: " + ex); }
            _done = true;
            Log("work: closing form to end modal loop...");
            try { form.Close(); } catch {}
        });
        // supervise: wait up to 180 s for the UI flow to finish
        sw.Restart();
        while (!_done && sw.Elapsed.TotalSeconds < 180) Thread.Sleep(500);
        if (!_done) { Log("TIMEOUT: UI thread did not finish"); _stop = true; Thread.Sleep(500); Environment.Exit(5); }
        _stop = true;
        Thread.Sleep(600);
        _stop = true;
        Thread.Sleep(600);
        ReportAsims();
        Log("DONE");
        Environment.Exit(0);
    }

    static string SafeProp(string name) {
        try {
            object v = CallAny(_app, name, null);
            return v == null ? "<null>" : v.ToString();
        } catch (Exception ex) { return "<ex: " + ex.Message + ">"; }
    }

    static bool BoolProp(string name) {
        try {
            object v = CallAny(_app, name, null);
            return v is bool && (bool)v;
        } catch { return false; }
    }

    static void ExportNow() {
        try {
            CallAny(_app, "ExportStandAloneSim", new object[]{ _aproj });
            Log("export returned");
        } catch (Exception ex) { Log("EXPORT-FAILED: " + ex.InnerException + " / " + ex.Message); }
    }

    static volatile bool _done;
    static volatile bool _workQueued;

    static void UiMain() {
        try {
            var asm = Assembly.LoadFrom(Path.Combine(_bin, "AnimatGUI.dll"));
            var appType = asm.GetType("AnimatGUI.Forms.AnimatApplication");
            Log("creating AnimatApplication...");
            _app = Activator.CreateInstance(appType);
            Log("app OK");

            Log("AnimatServer.Server.Initialize...");
            try {
                var srvAsm = Assembly.LoadFrom(Path.Combine(_bin, "AnimatServer.dll"));
                var srvType = srvAsm.GetType("AnimatServer.Server");
                const BindingFlags ALL2 = BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static;
                var init = srvType.GetMethod("Initialize", ALL2);
                object srv = Activator.CreateInstance(srvType);
                init.Invoke(srv, new object[]{ _app, 0 });
                Log("Server.Initialize OK");
            } catch (Exception ex) { Log("Server.Initialize failed: " + ex.Message); }

            Log("StartApplication (MODAL - runs its own loop)...");
            _workQueued = false;
            CallAny(_app, "StartApplication", new object[]{ true, false, null });
            Log("StartApplication returned (modal loop ended)");
        } catch (Exception ex) {
            Log("FATAL(UI): " + ex);
            _done = true;
            try { System.Windows.Forms.Application.ExitThread(); } catch {}
        }
    }

    static void ReportAsims() {
        foreach (string d in new[]{ Path.GetDirectoryName(_aproj), _bin, Directory.GetCurrentDirectory(),
                                    Path.Combine(Path.GetDirectoryName(_aproj), "..") }) {
            try {
                if (!Directory.Exists(d)) continue;
                foreach (string f in Directory.GetFiles(d, "*.asim"))
                    Log("ASIM: " + f + " " + new FileInfo(f).Length + " " + File.GetLastWriteTime(f));
            } catch {}
        }
    }

    static void KillErrorDialogs() {
        var seen = new HashSet<IntPtr>();
        while (!_stop) {
            try {
                EnumWindows((h, lp) => {
                    uint wpid; GetWindowThreadProcessId(h, out wpid);
                    if (wpid == _pid && IsWindowVisible(h) && !seen.Contains(h)) {
                        var cn = new StringBuilder(256); GetClassName(h, cn, 256);
                        var ti = new StringBuilder(512); GetWindowText(h, ti, 512);
                        if (cn.ToString() == "#32770" && ti.ToString() == "Error") {
                            seen.Add(h);
                            Thread.Sleep(1500);  // give outside observers time to read it
                            EnumChildWindows(h, (c, lp2) => {
                                var ccn = new StringBuilder(256); GetClassName(c, ccn, 256);
                                var cti = new StringBuilder(1024); GetWindowText(c, cti, 1024);
                                Log("   ERRDLG child " + ccn + ": '" + cti + "'");
                                if (ccn.ToString() == "Button" && cti.ToString() == "OK") {
                                    SendMessage(c, 0x00F5, IntPtr.Zero, IntPtr.Zero); // BM_CLICK
                                    Log("   ERRDLG: clicked OK");
                                }
                                return true;
                            }, IntPtr.Zero);
                        }
                    }
                    return true;
                }, IntPtr.Zero);
            } catch {}
            Thread.Sleep(300);
        }
    }
}
