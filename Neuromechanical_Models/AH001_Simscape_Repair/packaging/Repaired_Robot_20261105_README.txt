REPAIRED ROUTING CANDIDATES — NOT A QUALIFIED FINAL DESIGN

Blender: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\packaging\Repaired_Robot_20261105.blend
Open this .blend directly. The saved scene is REPAIRED COMPLIANT; the scene
selector also offers REPAIRED RIGID. Numpad 0 selects camera view; Space plays
frames 1–289. Knees cycle from 0 to -120 degrees and back; camera orbits 180
degrees about Blender Z. World/render background is white.

Collections 01/02/03 contain native OpenSim bones, gray original right-side
routes and blue repaired routes. Collection 04 contains full solved BPA
channels/sections. Collection 05 contains fitting/tackle envelopes, with
unresolved tackle wireframes excluded from rendering. Select a BPA and read
its custom properties for nominal diameter, parallel count, gain, status,
contraction, length and rendered diameter by knee pose.

OpenSim compliant model: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\packaging\gait2392_robotbody_repaired_compliant_blue.osim
OpenSim rigid model: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\packaging\gait2392_robotbody_repaired_rigid_blue.osim
Load the model, then motion: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\packaging\knee_comparison_20261105.mot
Gray = original paths. Blue = disabled, display-only optimized PathActuators.
Original model files were not overwritten.

Torque, clearance and requirement audit: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\audit_20261105\ROUTING_AUDIT_20261105.md
Torque figure: C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\audit_20261105\Repaired_Torque_Comparison.pdf
Exact tables, source hashes and native API checks are in that audit folder.

Some candidates fail torque or clearance. The cam, complete physical reeving,
independent series-section mechanics and full robot-CAD qualification remain
unresolved. These files expose failures for design review; they do not hide
them or certify packaging for a dissertation defense.
