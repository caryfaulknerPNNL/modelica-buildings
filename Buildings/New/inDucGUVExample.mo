within ;
model inDucGUVExample
  replaceable package MediumA = Buildings.Media.Air(extraPropertiesNames={"CO2", "COVID"}) "Medium model for air";
  parameter Modelica.Units.SI.MassFlowRate mAir_flow_nominal = 6.71
    "Design mass flow rate";
  Buildings.Fluid.FixedResistances.InDuctGUV
                                   inDucGUV(
    dp_nominal=0,
    eff=1,
    kpow=200,
    m_flow_nominal=mAir_flow_nominal,
    redeclare package Medium = MediumA)
    annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
  Buildings.Fluid.Sources.MassFlowSource_T boundary(
    use_C_in=false,
    C={0,1},
    use_m_flow_in=true,
    redeclare package Medium = MediumA,
    nPorts=1) annotation (Placement(transformation(extent={{-58,-12},{-38,8}})));
  Buildings.Fluid.Sources.Boundary_pT bou(nPorts=1, redeclare package Medium = MediumA)
    annotation (Placement(transformation(extent={{72,-10},{52,10}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=0.75*mAir_flow_nominal,
    duration=3600,
    offset=0.25*mAir_flow_nominal,
    startTime=120)
    annotation (Placement(transformation(extent={{-88,26},{-68,46}})));
  Modelica.Blocks.Sources.BooleanConstant booleanConstant
    annotation (Placement(transformation(extent={{-66,-46},{-46,-26}})));
equation
  connect(boundary.ports[1], inDucGUV.port_a)
    annotation (Line(points={{-38,-2},{-10,-2}}, color={0,127,255}));
  connect(inDucGUV.port_b, bou.ports[1]) annotation (Line(points={{10,-2},{46,-2},
          {46,0},{52,0}}, color={0,127,255}));
  connect(booleanConstant.y, inDucGUV.u) annotation (Line(points={{-45,-36},{-20,
          -36},{-20,-10},{-12,-10}}, color={255,0,255}));
  connect(ramp.y, boundary.m_flow_in) annotation (Line(points={{-67,36},{-62,36},
          {-62,12},{-66,12},{-66,6},{-60,6}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Buildings(version="10.0.0"), Modelica(version="4.0.0")));
end inDucGUVExample;
