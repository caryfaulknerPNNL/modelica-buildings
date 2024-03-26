within ;
model GUVExample
  replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "COVID"}) "Medium model for air";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6.71
    "Design mass flow rate";
  Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                           vol(
    redeclare package Medium = Medium,
    final energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
    final massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
    final V=500,
    final C_start={0,1},
    final mSenFac=1,
    final m_flow_nominal=m_flow_nominal,
    final prescribedHeatFlowRate=true,
    final nPorts=2,
    m_flow_small=1E-4*abs(m_flow_nominal),
    allowFlowReversal=true,
    final use_C_flow=true)        "Room air volume"
    annotation (Placement(transformation(extent={{8,-8},{28,12}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-60,22},{-40,42}})));
  Modelica.Blocks.Sources.RealExpression patConc(y=vol.C[2])
    "virus concentration"
    annotation (Placement(transformation(extent={{-90,-40},{-70,-20}})));
  Buildings.Fluid.Sources.MassFlowSource_T sup(
    use_C_in=false,
    C={0,1},
    use_m_flow_in=false,
    redeclare package Medium = Medium,
    m_flow=0*500*1.2/3600,
    nPorts=1)
    annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
  Buildings.Fluid.Sources.Boundary_pT bou(nPorts=1, redeclare package Medium =
        Medium)
    annotation (Placement(transformation(extent={{58,-38},{38,-18}})));
  Modelica.Blocks.Sources.BooleanConstant booleanConstant
    annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
  ViralDecay viralDecay
    annotation (Placement(transformation(extent={{-36,-90},{-16,-70}})));
  Modelica.Blocks.Sources.RealExpression patConc1(y=vol.C[1])
    "virus concentration"
    annotation (Placement(transformation(extent={{-62,-30},{-42,-10}})));
  GUV gUV(redeclare package Medium = Medium, krad={0,2.9e3})
    annotation (Placement(transformation(extent={{-30,-56},{-10,-36}})));
equation
  connect(const.y, vol.mWat_flow) annotation (Line(points={{-39,32},{-22,32},{
          -22,10},{6,10}},
                       color={0,0,127}));
  connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{-4,0},{
          -4,-12},{16,-12},{16,-8}}, color={0,127,255}));
  connect(vol.ports[2], bou.ports[1])
    annotation (Line(points={{20,-8},{20,-28},{38,-28}},  color={0,127,255}));
  connect(patConc.y, viralDecay.u) annotation (Line(points={{-69,-30},{-44,-30},
          {-44,-80},{-38,-80}}, color={0,0,127}));
  connect(booleanConstant.y, gUV.u) annotation (Line(points={{-69,-62},{-38,-62},
          {-38,-48},{-32,-48}}, color={255,0,255}));
  connect(patConc1.y, gUV.C[1]) annotation (Line(points={{-41,-20},{-34,-20},{
          -34,-36},{-38,-36},{-38,-42},{-32,-42}}, color={0,0,127}));
  connect(patConc.y, gUV.C[2]) annotation (Line(points={{-69,-30},{-44,-30},{
          -44,-42},{-32,-42}}, color={0,0,127}));
  connect(gUV.yC_flow[1], vol.C_flow[1]) annotation (Line(points={{-9,-42},{0,
          -42},{0,-4},{6,-4}}, color={0,0,127}));
  connect(gUV.yC_flow[2], vol.C_flow[2]) annotation (Line(points={{-9,-42},{0,
          -42},{0,-4},{6,-4}}, color={0,0,127}));
  connect(gUV.diss, vol.heatPort) annotation (Line(points={{-9.6,-54},{-2,-54},
          {-2,2},{8,2}}, color={191,0,0}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Buildings(version="10.0.0"), Modelica(version="4.0.0")));
end GUVExample;
