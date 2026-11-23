import java.nio.file.*;
import org.opensim.modeling.*;
import org.opensim.threejs.ModelVisualizationJson;
import org.opensim.view.pub.OpenSimDB;
import org.json.simple.JSONObject;
class NativeScene {
  public static void main(String[] args) throws Exception {
    System.loadLibrary("osimJavaJNI");
    ModelVisualizer.addDirToGeometrySearchPaths("C:/OpenSim 4.6/Geometry");
    Path path=Path.of(args[0]);
    ModelVisualizer.addDirToGeometrySearchPaths(path.getParent().resolve("Geometry").toString());
    Model model=new Model(path.toString());
    for (int i=0; i<model.getMuscles().getSize(); i++) { model.getMuscles().get(i).getGeometryPath().upd_Appearance().set_visible(true); }
    State state=model.initSystem();
    model.realizePosition(state);
    OpenSimContext context=new OpenSimContext(state, model);
    OpenSimDB.getInstance().setContext(model,context);
    ModelVisualizationJson view=new ModelVisualizationJson(new JSONObject(),model);
    Files.writeString(Path.of(args[1]),view.toJSONString());
    System.out.println("OpenSim Desktop model export: "+model.getName()+" ("+model.getMuscles().getSize()+" muscles)");
    System.exit(0);
  }
}
