within ;
model ViralDecay

  parameter Real kdec(min=0) = 0.48
    "Decay rate of virus";

  parameter Real V(min=0) = 100
    "Room volume";
  Modelica.Blocks.Math.Gain gain(k=-kdec*1.2*V/3600)
    annotation (Placement(transformation(extent={{-6,-10},{14,10}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
equation
  connect(u, gain.u)
    annotation (Line(points={{-120,0},{-8,0}}, color={0,0,127}));
  connect(gain.y, y)
    annotation (Line(points={{15,0},{110,0}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")));
end ViralDecay;
