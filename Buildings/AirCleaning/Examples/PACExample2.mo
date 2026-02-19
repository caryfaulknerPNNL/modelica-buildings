within Buildings.AirCleaning.Examples;
model PACExample2
  extends Modelica.Icons.Example;
  replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "SARS-CoV-2"}) "Medium model for air";
  parameter Modelica.Units.SI.MassFlowRate m_flow_nominal=1.2*vol.V/3600
    "Design mass flow rate";
  Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir
                                           vol(
    redeclare package Medium = Medium,
    final energyDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
    final massDynamics=Modelica.Fluid.Types.Dynamics.DynamicFreeInitial,
    final V=500,
    final C_start={0,0},
    final mSenFac=1,
    final m_flow_nominal=m_flow_nominal,
    final prescribedHeatFlowRate=true,
    final nPorts=3,
    m_flow_small=1E-4*abs(m_flow_nominal),
    allowFlowReversal=true,
    final use_C_flow=true)        "Room air volume"
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  Modelica.Blocks.Sources.Constant conZero(k=0)
    annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  Modelica.Blocks.Sources.RealExpression speConc2(y=vol.C[2])
    "SARS-CoV-2 concentration"
    annotation (Placement(transformation(extent={{-90,-50},{-70,-30}})));
  Buildings.Fluid.Sources.MassFlowSource_T sup(
    use_C_in=false,
    C={0.1,0.1},
    use_m_flow_in=false,
    redeclare package Medium = Medium,
    m_flow=m_flow_nominal,
    nPorts=1) "Supply airflow"
    annotation (Placement(transformation(extent={{-68,-10},{-48,10}})));
  Modelica.Blocks.Sources.BooleanConstant conPACEna "True when PAC is on"
    annotation (Placement(transformation(extent={{-90,-72},{-70,-52}})));
  Modelica.Blocks.Sources.RealExpression speConc1(y=vol.C[1])
    "CO2 concentration"
    annotation (Placement(transformation(extent={{-90,-30},{-70,-10}})));
  PAC_new pAC_new(redeclare package Medium = Medium, m_flow_in=1)
    annotation (Placement(transformation(extent={{20,-50},{40,-30}})));
equation
  connect(sup.ports[1], vol.ports[1]) annotation (Line(points={{-48,0},{-8,0},{
          -8,-12},{48.6667,-12},{48.6667,-10}},
                                     color={0,127,255}));
  connect(conZero.y, vol.mWat_flow)
    annotation (Line(points={{-19,50},{38,50},{38,8}}, color={0,0,127}));
  connect(conZero.y, vol.C_flow[1]) annotation (Line(points={{-19,50},{38,50},{
          38,16},{30,16},{30,-6},{38,-6}}, color={0,0,127}));
  connect(conZero.y, vol.C_flow[2]) annotation (Line(points={{-19,50},{38,50},{
          38,16},{30,16},{30,-6},{38,-6}}, color={0,0,127}));
  connect(vol.ports[2], pAC_new.port_a) annotation (Line(points={{50,-10},{50,
          -16},{2,-16},{2,-40},{20,-40}}, color={0,127,255}));
  connect(pAC_new.port_b, vol.ports[3]) annotation (Line(points={{40,-40},{
          51.3333,-40},{51.3333,-10}}, color={0,127,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>This example simulates <a href=\"modelica://Buildings.AirCleaning.PAC\">AirCleaning.PAC</a> for a scenario with CO<sub>2. </sub>and SARS-CoV-2 as trace species. Airflow with a constant rate and concentration of trace species is supplied to a mixing volume. A PAC is connected to the mixing volume to simulate removal of trace species from the volume, where the PAC cannot remove the gaseous CO<sub>2</sub> but filters SARS-CoV-2 with an efficiency of 99&percnt;.</p>
</html>"),
    experiment(StopTime=3600, __Dymola_Algorithm="Dassl"));
end PACExample2;
