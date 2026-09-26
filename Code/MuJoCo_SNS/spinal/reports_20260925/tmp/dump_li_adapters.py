import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
root = ET.parse(path).getroot()

def txt(el, tag):
    x = el.find(tag)
    return x.text if x is not None else None

nm = root.find('.//NervousSystem/Node')
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Adapter' not in cls and not cls.endswith('.RigidBody') and not cls.endswith('.Joint'):
        continue
    name = txt(nd, 'Text') or txt(nd, 'Name') or ''
    if 'Adapter' in cls:
        print("=" * 80)
        print(f"ADAPTER {name!r} cls={cls.split('.')[-1]} id={txt(nd,'ID')}")
        # dump all direct children compactly, skipping cosmetic ones
        skip = {'AssemblyFile', 'ClassName', 'ID', 'Alignment', 'AutoSize', 'BackMode',
                'DashStyle', 'DrawColor', 'DrawWidth', 'Font', 'Hidden', 'Jump',
                'LineStyle', 'OrthogonalDynamic', 'OrientedText', 'Selectable',
                'Stretchable', 'Text', 'ToolTip', 'Url', 'ZOrder', 'FillColor',
                'Gradient', 'GradientColor', 'GradientMode', 'DiagramImageName',
                'ImageName', 'ImageLocation', 'ImagePosition', 'InLinkable',
                'LabelEdit', 'Location', 'OutLinkable', 'ShadowStyle', 'ShadowColor',
                'ShadowSize', 'Shape', 'ShapeOrientation', 'TextColor', 'TextMargin',
                'Transparent', 'XMoveable', 'XSizeable', 'YMoveable', 'YSizeable',
                'InLinks', 'OutLinks', 'TemplateNode', 'TemplateNodeCount',
                'TemplateChangeScript', 'Enabled', 'SynchWithRobot',
                'SynchUpdateInterval', 'InitIODisableDuration', 'RobotIOScale',
                'DelayBufferMode', 'DelayBufferInterval'}
        for ch in nd:
            if ch.tag in skip:
                continue
            print("  <" + ch.tag + ">",
                  {k: v for k, v in ch.attrib.items()},
                  (ch.text or '').strip()[:60])
            for gc in ch:
                print("      <" + gc.tag + ">", {k: v for k, v in gc.attrib.items()},
                      (gc.text or '').strip()[:60])
                for ggc in gc:
                    print("          <" + ggc.tag + ">", dict(ggc.attrib),
                          (ggc.text or '').strip()[:60])
