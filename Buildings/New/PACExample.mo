within ;
model PACExample
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
  PAC pAC(
    redeclare package Medium = Medium,
    eff={0,0.9997},
    nPACs=2,
    flowPAC=400*1.2/3600,
    kpow=100)
    annotation (Placement(transformation(extent={{-36,-62},{-16,-42}})));
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
equation
  connect(const.y, vol.mWat_flow) annotation (Line(points={{-39,32},{-22,32},{
          -22,10},{6,10}},
                       color={0,0,127}));
  connect(pAC.diss, vol.heatPort) annotation (Line(points={{-15.6,-60},{-6,-60},
          {-6,2},{8,2}},
                color={191,0,0}));
  connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{-4,0},{
          -4,-12},{16,-12},{16,-8}}, color={0,127,255}));
  connect(vol.ports[2], bou.ports[1])
    annotation (Line(points={{20,-8},{20,-28},{38,-28}},  color={0,127,255}));
  connect(booleanConstant.y, pAC.u) annotation (Line(points={{-69,-62},{-46,-62},
          {-46,-54},{-38,-54}}, color={255,0,255}));
  connect(patConc.y, viralDecay.u) annotation (Line(points={{-69,-30},{-44,-30},
          {-44,-80},{-38,-80}}, color={0,0,127}));
  connect(pAC.yC_flow[1], vol.C_flow[1]) annotation (Line(points={{-15,-48},{0,
          -48},{0,-10},{2,-10},{2,-4},{6,-4}}, color={0,0,127}));
  connect(pAC.yC_flow[2], vol.C_flow[2]) annotation (Line(points={{-15,-48},{0,
          -48},{0,-10},{2,-10},{2,-4},{6,-4}}, color={0,0,127}));
  connect(patConc1.y, pAC.C[1]) annotation (Line(points={{-41,-20},{-38,-20},{
          -38,-38},{-46,-38},{-46,-48},{-38,-48}}, color={0,0,127}));
  connect(patConc.y, pAC.C[2]) annotation (Line(points={{-69,-30},{-44,-30},{
          -44,-48},{-38,-48}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Buildings(version="10.0.0"), Modelica(version="4.0.0")));
end PACExample;
