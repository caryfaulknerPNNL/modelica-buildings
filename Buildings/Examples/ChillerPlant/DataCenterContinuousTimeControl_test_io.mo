within Buildings.Examples.ChillerPlant;
model DataCenterContinuousTimeControl_test_io
  "Model of data center that approximates the trim and respond logic"
  extends Buildings.Examples.ChillerPlant.BaseClasses.DataCenter(
    break connect(cooCoi.port_a1, val6.port_b),
    break connect(val4.port_b, cooTow.port_a),
    break feedback,
    break TAirSet,
    break linPieTwo,
    break gain,
    break chiSwi,
    break chiCon,
    break KMinusU,
    break connect(cooTowFanCon.y, cooTow.y),
    break connect(val5.port_b, cooTow.port_a));
  extends Modelica.Icons.Example;

  Modelica.Blocks.Sources.Constant dpSet(k=45000) annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, origin={38,80})));
  Modelica.Blocks.Sources.Constant CHWSet(k=273.15 + 6) annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, origin={80,100})));
  Fluid.Sensors.TemperatureTwoPort TCHWSup(redeclare package Medium = MediumW,
      m_flow_nominal=mCHW_flow_nominal)
    "Temperature of chilled water leaving the cooling coil" annotation (
      Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={348,-90})));
  Fluid.Sensors.TemperatureTwoPort TCWReturn(redeclare package Medium = MediumW,
      m_flow_nominal=mCW_flow_nominal)
    "Temperature of condenser water leaving the cooling tower" annotation (
      Placement(transformation(extent={{10,-10},{-10,10}}, origin={184,211})));
  Modelica.Blocks.Sources.Constant zero(k=0)
    "Set temperature for air supply to the room" annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, origin={-150,158})));
  Modelica.Blocks.Sources.Constant one(k=1)
    "Set temperature for air supply to the room" annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, origin={-152,190})));
  Modelica.Blocks.Sources.BooleanConstant booleanConstant
    annotation (Placement(transformation(extent={{-62,74},{-42,94}})));
  Controls.Continuous.LimPID conPID(Ti=40)
    annotation (Placement(transformation(extent={{158,256},{178,276}})));
  Modelica.Blocks.Sources.Constant CWSet(k=273.15 + 20) annotation (Placement(
        transformation(extent={{-10,-10},{10,10}}, origin={112,268})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
    table=[0,0; 1,2; 2,4; 3,6; 4,8; 5,10; 6,8; 7,6; 8,4; 9,2],
    smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
    extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic,
    timeScale=3600*4,
    offset={273.15 + 20})
    annotation (Placement(transformation(extent={{-190,228},{-170,248}})));
  Modelica.Blocks.Interfaces.RealInput u_CWST
    annotation (Placement(transformation(extent={{-440,120},{-400,160}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{400,94},{420,114}})));
  Utilities.IO.SignalExchange.Read read
    annotation (Placement(transformation(extent={{372,72},{392,92}})));
  Modelica.Blocks.Interfaces.RealInput u_Qflow
    annotation (Placement(transformation(extent={{-442,40},{-402,80}})));
equation

  connect(dpSet.y, pumCHW.dp_in) annotation (Line(points={{49,80},{56,80},{56,-24},
          {64,-24},{64,-120},{206,-120}}, color={0,0,127}));
  connect(CHWSet.y, chi.TSet) annotation (Line(points={{91,100},{100,100},{100,84},
          {96,84},{96,76},{284,76},{284,90},{276,90}}, color={0,0,127}));
  connect(val6.port_b, TCHWSup.port_b) annotation (Line(points={{358,30},{360,30},
          {360,20},{348,20},{348,-80}}, color={0,127,255}));
  connect(TCHWSup.port_a, cooCoi.port_a1) annotation (Line(points={{348,-100},{348,
          -164},{300,-164}}, color={0,127,255}));
  connect(zero.y, valByp.y) annotation (Line(points={{-139,158},{224,158},{224,44},
          {288,44},{288,32}}, color={0,0,127}));
  connect(one.y, val6.y) annotation (Line(points={{-141,190},{-24,190},{-24,164},
          {228,164},{228,48},{336,48},{336,40},{346,40}}, color={0,0,127}));
  connect(booleanConstant.y, chi.on) annotation (Line(points={{-41,84},{20,84},
          {20,64},{96,64},{96,68},{136,68},{136,52},{216,52},{216,40},{292,40},
          {292,96},{276,96}}, color={255,0,255}));
  connect(one.y, val5.y) annotation (Line(points={{-141,190},{-24,190},{-24,164},
          {200,164},{200,180},{206,180}}, color={0,0,127}));
  connect(booleanConstant.y, or1.u2)
    annotation (Line(points={{-41,84},{18,84},{18,192}}, color={255,0,255}));
  connect(conPID.y, cooTow.y) annotation (Line(points={{179,266},{208,266},{208,
          247},{257,247}}, color={0,0,127}));
  connect(TCWLeaTow.T, conPID.u_s) annotation (Line(points={{330,130},{330,136},
          {380,136},{380,292},{156,292},{156,266}}, color={0,0,127}));
  connect(TCHWSup.T, read.u) annotation (Line(points={{337,-90},{337,-92},{328,
          -92},{328,92},{364,92},{364,82},{370,82}}, color={0,0,127}));
  connect(read.y, y) annotation (Line(points={{393,82},{428,82},{428,104},{410,
          104}}, color={0,0,127}));
  connect(val4.port_b, cooTow.port_a) annotation (Line(points={{98,190},{98,232},
          {220,232},{220,239},{259,239}}, color={0,0,127}));
  connect(val5.port_b, TCWReturn.port_b) annotation (Line(points={{218,190},{
          218,228},{192,228},{192,236},{164,236},{164,211},{174,211}}, color={0,
          127,255}));
  connect(TCWReturn.port_a, cooTow.port_a) annotation (Line(points={{194,211},{
          240,211},{240,239},{259,239}}, color={0,127,255}));
  connect(u_Qflow, roo.uQflow) annotation (Line(points={{-422,60},{-368,60},{
          -368,-230},{228,-230},{228,-240},{236,-240}}, color={0,0,127}));
  connect(u_CWST, conPID.u_m) annotation (Line(points={{-420,140},{-204,140},{
          -204,256},{96,256},{96,252},{152,252},{152,248},{168,248},{168,254}},
        color={0,0,127}));
  annotation (
    __Dymola_Commands(file=
          "modelica://Buildings/Resources/Scripts/Dymola/Examples/ChillerPlant/DataCenterContinuousTimeControl.mos"
        "Simulate and plot"), Documentation(info="<html>
<p>
This model is the chilled water plant with continuous time control.
The trim and respond logic is approximated by a PI controller which
significantly reduces computing time. The model is described at
<a href=\"Buildings.Examples.ChillerPlant\">
Buildings.Examples.ChillerPlant</a>.
</p>
<p>
See
<a href=\"Buildings.Examples.ChillerPlant.DataCenterContinuousTimeControl\">
Buildings.Examples.ChillerPlant.DataCenterContinuousTimeControl</a>
for an implementation with the discrete time trim and respond logic.
</p>
</html>", revisions="<html>
<ul>
<li>
January 13, 2015, by Michael Wetter:<br/>
Moved base model to
<a href=\"Buildings.Examples.ChillerPlant.BaseClasses.DataCenter\">
Buildings.Examples.ChillerPlant.BaseClasses.DataCenter</a>.
</li>
<li>
December 5, 2012, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-400,-300},{400,
            300}})),
    experiment(
      StopTime=31536000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end DataCenterContinuousTimeControl_test_io;
