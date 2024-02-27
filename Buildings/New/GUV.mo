within ;
model GUV "In-room GUV"

  parameter Real frad(min=0, max=1) = 0.2
    "Fraction of irradiated space";

  parameter Real Eavg(min=0, max=1) = 50e-6
    "Effluence rate";

  parameter Real krad(min=0) = 2.93e3
    "Inactivation constant";

  parameter Real kpow(min=0) = 120
    "Rated power";

  parameter Real V(min=0) = 120
    "Zone volume";

  Modelica.Blocks.Interfaces.RealInput C "Zone concentration"
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealOutput yC_flow "Concentration outflow"
    annotation (Placement(transformation(extent={{100,30},{120,50}})));
  Modelica.Blocks.Interfaces.BooleanInput u "on/off"
    annotation (Placement(transformation(extent={{-140,-40},{-100,0}})));
  Modelica.Blocks.Interfaces.RealOutput yP_GUV "Power output"
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Interfaces.RealOutput yE_GUV "Energy output"
    annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal
    annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
  Modelica.Blocks.Continuous.Integrator E_GUV(use_reset=false)
    annotation (Placement(transformation(extent={{60,-50},{80,-30}})));
  Modelica.Blocks.Math.MultiProduct multiProduct(nu=2)
    annotation (Placement(transformation(extent={{-80,34},{-68,46}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a diss "Heat dissipation"
    annotation (Placement(transformation(extent={{94,-90},{114,-70}})));

protected
  Modelica.Blocks.Math.Gain pow(final k=kpow)
                                             "power of GUV"
    annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
  Modelica.Blocks.Math.Gain UV_calc(final k=-1.2*V*frad*Eavg*krad)
    "death of virus from UV"
    annotation (Placement(transformation(extent={{-16,30},{4,50}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prePow(final alpha=0)
    "Prescribed power (=heat and flow work) flow for dynamic model"
    annotation (Placement(transformation(extent={{10,-70},{30,-50}})));

equation
  connect(u, booleanToReal.u)
    annotation (Line(points={{-120,-20},{-82,-20}}, color={255,0,255}));
  connect(C, multiProduct.u[1]) annotation (Line(points={{-120,40},{-100,40},{
          -100,42.1},{-80,42.1}}, color={0,0,127}));
  connect(booleanToReal.y, multiProduct.u[2]) annotation (Line(points={{-59,-20},
          {-54,-20},{-54,28},{-84,28},{-84,37.9},{-80,37.9}}, color={0,0,127}));
  connect(multiProduct.y, UV_calc.u)
    annotation (Line(points={{-66.98,40},{-18,40}}, color={0,0,127}));
  connect(booleanToReal.y, pow.u)
    annotation (Line(points={{-59,-20},{-22,-20}}, color={0,0,127}));
  connect(pow.y, yP_GUV) annotation (Line(points={{1,-20},{96,-20},{96,0},{110,
          0}}, color={0,0,127}));
  connect(pow.y, prePow.Q_flow) annotation (Line(points={{1,-20},{8,-20},{8,-46},
          {2,-46},{2,-60},{10,-60}}, color={0,0,127}));
  connect(prePow.port, diss) annotation (Line(points={{30,-60},{90,-60},{90,-80},
          {104,-80}}, color={191,0,0}));
  connect(pow.y, E_GUV.u) annotation (Line(points={{1,-20},{8,-20},{8,-40},{58,
          -40}}, color={0,0,127}));
  connect(E_GUV.y, yE_GUV)
    annotation (Line(points={{81,-40},{110,-40}}, color={0,0,127}));
  connect(UV_calc.y, yC_flow)
    annotation (Line(points={{5,40},{110,40}}, color={0,0,127}));
  annotation (Icon(graphics={
        Polygon(
          points={{-60,58},{-60,-64},{-40,-64},{-40,58},{-60,58}},
          lineColor={28,108,200},
          fillColor={28,108,200},
          fillPattern=FillPattern.Solid),
        Line(points={{-36,48},{42,48}}, color={28,108,200}),
        Line(points={{-36,26},{42,26}}, color={28,108,200}),
        Line(points={{-36,8},{42,8}},   color={28,108,200}),
        Line(points={{-36,-32},{42,-32}}, color={28,108,200}),
        Line(points={{-36,-52},{42,-52}}, color={28,108,200}),
        Line(points={{-36,-12},{42,-12}}, color={28,108,200})}), uses(Buildings(
          version="10.0.0")));
end GUV;
