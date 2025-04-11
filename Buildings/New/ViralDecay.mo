within ;
model ViralDecay

  replaceable package Medium =
    Modelica.Media.Interfaces.PartialMedium "Medium in the component";

  parameter Real kdec[Medium.nC](min=0) = {0,0.48}
    "Decay rate of virus";

  parameter Real V(min=0) = 100
    "Room volume";
  Modelica.Blocks.Interfaces.RealInput[Medium.nC] u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealOutput[Medium.nC] y
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
equation
  for i in 1:Medium.nC loop
    y[i] = -1.2*V*kdec[i]*u[i]/3600;
  end for;
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Modelica(version="4.0.0")));
end ViralDecay;
