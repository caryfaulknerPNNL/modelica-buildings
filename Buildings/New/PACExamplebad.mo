within ;
model PACExamplebad
  replaceable package Medium = Buildings.Media.Air(extraPropertiesNames={"CO2", "COVID"}) "Medium model for air";
  parameter Modelica.Units.SI.Length wExtSou=49.91
    "South zone exterior wall width";
  parameter Modelica.Units.SI.Length wExtNor=49.91
    "North zone exterior wall width";
  parameter Modelica.Units.SI.Length wExtEas=33.27
    "East zone exterior wall width";
  parameter Modelica.Units.SI.Length wExtWes=33.27
    "West zone exterior wall width";

  parameter Buildings.HeatTransfer.Types.InteriorConvection intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature
    "Convective heat transfer model for room-facing surfaces of opaque constructions";
  parameter Real winWalRat(
    min=0.01,
    max=0.99) = 0.33 "Window to wall ratio for exterior walls";
  parameter Modelica.Units.SI.Length hWin=1.5 "Height of windows";
  constant Modelica.Units.SI.Height hRoo=2.74 "Room height";
  parameter Boolean sampleModel = false
   "Set to true to time-sample the model, which can give shorter simulation time if there is already time sampling in the system model";
  parameter Modelica.Units.SI.Area AFloSou=568.77/hRoo;

  Buildings.ThermalZones.Detailed.MixedAir_virus sou(
    redeclare package Medium = Medium,
    AFlo=AFloSou,
    hRoo=hRoo,
    nConExt=0,
    nConExtWin=1,
    datConExtWin(
      layers={conExtWal},
      A={wExtSou*hRoo},
      glaSys={glaSys},
      wWin={winWalRat/hWin*wExtSou*hRoo},
      each hWin=hWin,
      fFra={0.1},
      til={Buildings.Types.Tilt.Wall},
      azi={Buildings.Types.Azimuth.S}),
    nConPar=2,
    datConPar(
      layers={conFlo,conFur},
      A={AFloSou,414.68},
      til={Buildings.Types.Tilt.Floor,Buildings.Types.Tilt.Wall}),
    nConBou=3,
    datConBou(
      layers={conIntWal,conIntWal,conIntWal},
      A={6.47,40.76,6.47}*hRoo,
      til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall}),
    nSurBou=0,
    use_C_flow=true,
    kdec=0,
    frad=0,
    Eavg=0,
    krad=0,
    kpow_GUV=0,
    nPACs=1,
    flowPAC=0.133,
    kpow_PAC=0,
    intConMod=intConMod,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    final sampleModel=sampleModel,
    nPorts=5)                      "South zone"
    annotation (Placement(transformation(extent={{24,-24},{64,16}})));

  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-36,-14},{-16,6}})));
  parameter Buildings.HeatTransfer.Data.Solids.Plywood
                                             matFur(x=0.15, nStaRef=5)
    "Material for furniture"
    annotation (Placement(transformation(extent={{-42,66},{-22,86}})));
  parameter Buildings.HeatTransfer.Data.Resistances.Carpet
                                                 matCar "Carpet"
    annotation (Placement(transformation(extent={{-2,66},{18,86}})));
  parameter Buildings.HeatTransfer.Data.Solids.Concrete
                                              matCon(
    x=0.1,
    k=1.311,
    c=836,
    nStaRef=5) "Concrete"
    annotation (Placement(transformation(extent={{-36,44},{-16,64}})));
  parameter Buildings.HeatTransfer.Data.Solids.Plywood
                                             matWoo(
    x=0.01,
    k=0.11,
    d=544,
    nStaRef=1) "Wood for exterior construction"
    annotation (Placement(transformation(extent={{-80,8},{-60,28}})));
  parameter Buildings.HeatTransfer.Data.Solids.Generic
                                             matIns(
    x=0.087,
    k=0.049,
    c=836.8,
    d=265,
    nStaRef=5) "Steelframe construction with insulation"
    annotation (Placement(transformation(extent={{12,42},{32,62}})));
  parameter Buildings.HeatTransfer.Data.Solids.GypsumBoard
                                                 matGyp(
    x=0.0127,
    k=0.16,
    c=830,
    d=784,
    nStaRef=2) "Gypsum board"
    annotation (Placement(transformation(extent={{-82,-20},{-62,0}})));
  parameter Buildings.HeatTransfer.Data.Solids.GypsumBoard
                                                 matGyp2(
    x=0.025,
    k=0.16,
    c=830,
    d=784,
    nStaRef=2) "Gypsum board"
    annotation (Placement(transformation(extent={{-80,36},{-60,56}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic
                                                          conExtWal(final nLay=3,
      material={matWoo,matIns,matGyp}) "Exterior construction"
    annotation (Placement(transformation(extent={{74,62},{94,82}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic
                                                          conIntWal(final nLay=1,
      material={matGyp2}) "Interior wall construction"
    annotation (Placement(transformation(extent={{48,44},{68,64}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic
                                                          conFlo(final nLay=1,
      material={matCon})
                 "Floor construction (opa_a is carpet)"
    annotation (Placement(transformation(extent={{74,26},{94,46}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic
                                                          conFur(final nLay=1,
      material={matFur})
                 "Construction for internal mass of furniture"
    annotation (Placement(transformation(extent={{68,-4},{88,16}})));
  parameter Buildings.HeatTransfer.Data.Solids.Plywood
                                             matCarTra(
    k=0.11,
    d=544,
    nStaRef=1,
    x=0.215/0.11) "Wood for floor"
    annotation (Placement(transformation(extent={{-80,66},{-60,86}})));
  parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear
                                                                   glaSys(
    UFra=2,
    shade=Buildings.HeatTransfer.Data.Shades.Gray(),
    haveInteriorShade=false,
    haveExteriorShade=false) "Data record for the glazing system"
    annotation (Placement(transformation(extent={{36,66},{56,86}})));
  Modelica.Blocks.Sources.BooleanConstant booleanConstant
    annotation (Placement(transformation(extent={{-30,-54},{-10,-34}})));
  Buildings.Fluid.Sources.Boundary_pT bou(
    p(displayUnit="kPa") = 101000,                  redeclare package Medium =
        Medium,
    nPorts=1)
    annotation (Placement(transformation(extent={{46,-90},{26,-70}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam=
        Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos"),
      computeWetBulbTemperature=false) "Weather data reader"
    annotation (Placement(transformation(extent={{-62,-58},{-42,-38}})));
  Buildings.Fluid.Sources.Boundary_pT bou1(
    p(displayUnit="kPa") = 102000,
    redeclare package Medium = Medium,
    nPorts=1)
    annotation (Placement(transformation(extent={{78,-72},{58,-52}})));
  Modelica.Blocks.Sources.Constant const1(k=1)
    annotation (Placement(transformation(extent={{-16,-28},{4,-8}})));
  Buildings.Fluid.Sources.Boundary_pT bou2(
    p(displayUnit="kPa") = 101000,
    redeclare package Medium = Medium,
    nPorts=1)
    annotation (Placement(transformation(extent={{16,-92},{-4,-72}})));
  Buildings.Fluid.Sources.Boundary_pT bou3(
    p(displayUnit="kPa") = 101000,
    redeclare package Medium = Medium,
    nPorts=1)
    annotation (Placement(transformation(extent={{-10,-92},{-30,-72}})));
  Buildings.Fluid.Sources.Boundary_pT bou4(
    p(displayUnit="kPa") = 102000,
    redeclare package Medium = Medium,
    nPorts=1)
    annotation (Placement(transformation(extent={{94,-40},{74,-20}})));
  Buildings.Fluid.FixedResistances.PressureDrop res
    annotation (Placement(transformation(extent={{36,-42},{56,-22}})));
  Buildings.Fluid.FixedResistances.PressureDrop res1
    annotation (Placement(transformation(extent={{28,-60},{48,-40}})));
  Buildings.Fluid.FixedResistances.PressureDrop res2
    annotation (Placement(transformation(extent={{10,-70},{30,-50}})));
  Buildings.Fluid.FixedResistances.PressureDrop res3
    annotation (Placement(transformation(extent={{-18,-74},{2,-54}})));
  Buildings.Fluid.FixedResistances.PressureDrop res4
    annotation (Placement(transformation(extent={{-62,-84},{-42,-64}})));
equation
  connect(const.y, sou.uSha[1]) annotation (Line(points={{-15,-4},{14,-4},{14,14},
          {22.4,14}}, color={0,0,127}));
  connect(const.y, sou.uWin[1]) annotation (Line(points={{-15,-4},{14,-4},{14,9},
          {22.4,9}}, color={0,0,127}));
  connect(const.y, sou.qGai_flow[1]) annotation (Line(points={{-15,-4},{14,-4},{
          14,2.93333},{22.4,2.93333}}, color={0,0,127}));
  connect(const.y, sou.qGai_flow[2]) annotation (Line(points={{-15,-4},{14,-4},{
          14,4},{22.4,4}}, color={0,0,127}));
  connect(const.y, sou.qGai_flow[3]) annotation (Line(points={{-15,-4},{14,-4},{
          14,5.06667},{22.4,5.06667}}, color={0,0,127}));
  connect(booleanConstant.y, sou.u_on_off) annotation (Line(points={{-9,-44},{16,
          -44},{16,-18.2},{22.2,-18.2}}, color={255,0,255}));
  connect(weaDat.weaBus, sou.weaBus) annotation (Line(
      points={{-42,-48},{-34,-48},{-34,-18},{-40,-18},{-40,22},{60,22},{60,18},{
          61.9,18},{61.9,13.9}},
      color={255,204,51},
      thickness=0.5));
  connect(const.y, sou.C_flow[1]) annotation (Line(points={{-15,-4},{14,-4},{14,
          -1.2},{22.4,-1.2}}, color={0,0,127}));
  connect(const1.y, sou.C_flow[2]) annotation (Line(points={{5,-18},{14,-18},{14,
          -1.2},{22.4,-1.2}}, color={0,0,127}));
  connect(res.port_b, bou4.ports[1]) annotation (Line(points={{56,-32},{68,-32},
          {68,-30},{74,-30}}, color={0,127,255}));
  connect(res.port_a, sou.ports[1]) annotation (Line(points={{36,-32},{14,-32},{
          14,-16},{16,-16},{16,-17.2},{29,-17.2}}, color={0,127,255}));
  connect(res1.port_b, bou1.ports[1]) annotation (Line(points={{48,-50},{52,-50},
          {52,-62},{58,-62}}, color={0,127,255}));
  connect(res1.port_a, sou.ports[2]) annotation (Line(points={{28,-50},{22,-50},
          {22,-32},{14,-32},{14,-15.6},{29,-15.6}}, color={0,127,255}));
  connect(bou.ports[1], res2.port_b) annotation (Line(points={{26,-80},{20,-80},
          {20,-96},{-34,-96},{-34,-58},{4,-58},{4,-46},{24,-46},{24,-36},{30,-36},
          {30,-60}}, color={0,127,255}));
  connect(res2.port_a, sou.ports[3]) annotation (Line(points={{10,-60},{6,-60},{
          6,-32},{14,-32},{14,-14},{29,-14}}, color={0,127,255}));
  connect(bou2.ports[1], res3.port_b) annotation (Line(points={{-4,-82},{-4,-98},
          {-36,-98},{-36,-38},{-32,-38},{-32,-30},{-6,-30},{-6,-38},{8,-38},{8,-56},
          {4,-56},{4,-60},{2,-60},{2,-64}}, color={0,127,255}));
  connect(res3.port_a, sou.ports[4]) annotation (Line(points={{-18,-64},{-30,-64},
          {-30,-56},{-34,-56},{-34,-46},{-32,-46},{-32,-32},{14,-32},{14,-12.4},
          {29,-12.4}}, color={0,127,255}));
  connect(bou3.ports[1], res4.port_b) annotation (Line(points={{-30,-82},{-38,-82},
          {-38,-88},{-42,-88},{-42,-74}}, color={0,127,255}));
  connect(res4.port_a, sou.ports[5]) annotation (Line(points={{-62,-74},{-66,-74},
          {-66,-32},{14,-32},{14,-10.8},{29,-10.8}}, color={0,127,255}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    uses(Buildings(version="10.0.0"), Modelica(version="4.0.0")));
end PACExamplebad;
